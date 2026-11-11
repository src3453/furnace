/**
 * Furnace Tracker - multi-system chiptune tracker
 * Copyright (C) 2021-2026 tildearrow and contributors
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

// 3HS88PWN4 visualizer: live view of the 8 harmonic (8-operator) channels
// and the 4 PCM/WT/noise channels, drawn as one canvas.

#include "gui.h"
#include <imgui.h>
#include "imgui_internal.h"
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include "../engine/platform/s3hs.h"

namespace {

const ImU32 visBg=IM_COL32(10,13,18,255);
const ImU32 visPanel=IM_COL32(15,19,26,255);
const ImU32 visGrid=IM_COL32(44,52,64,255);
const ImU32 visDim=IM_COL32(120,130,145,255);
const ImU32 visText=IM_COL32(222,228,236,255);
const ImU32 visBarBg=IM_COL32(34,40,50,255);

const ImU32 visChanColors[12]={
  IM_COL32(64,128,255,255),  // HS1 blue
  IM_COL32(40,180,150,255),  // HS2 teal
  IM_COL32(150,100,235,255), // HS3 purple
  IM_COL32(240,130,50,255),  // HS4 orange
  IM_COL32(235,90,150,255),  // HS5 pink
  IM_COL32(220,200,70,255),  // HS6 yellow
  IM_COL32(90,200,100,255),  // HS7 green
  IM_COL32(70,200,230,255),  // HS8 cyan
  IM_COL32(255,120,90,255),  // PCM1 coral
  IM_COL32(245,170,60,255),  // PCM2 amber
  IM_COL32(170,215,70,255),  // PCM3 lime
  IM_COL32(175,120,240,255)  // PCM4 violet
};

const ImU32 visStageColors[5]={
  IM_COL32(240,205,80,255),  // attack
  IM_COL32(240,150,70,255),  // decay
  IM_COL32(80,205,120,255),  // sustain
  IM_COL32(230,90,90,255),   // release
  IM_COL32(90,98,110,255)    // off
};
const char* visStageNames[5]={"ATK","DEC","SUS","REL","OFF"};

// Algorithms as routed in hs816::Channel: name and the operators that reach
// the output (bit i = OPi+1).
const char* visAlgNames[13]={
  "ADD 8", "FM 4>4", "RING 4x4", "FM PAIRS", "FM 8 CHAIN", "FM 4x2", "FM 2x4",
  "FM 4x2 RING", "RING PAIRS", "FM 2x4 RING", "DIRECT 4>4", "DIRECT PAIRS", "DIRECT 8"
};
const unsigned char visAlgCarriers[13]={
  0xff, 0x0f, 0xff, 0x03, 0x01, 0x11, 0x55, 0x11, 0xff, 0x55, 0x0f, 0x03, 0x01
};

const char* visPcmModes[6]={"PCM","WAVE~","NOISE","NOISE64","WAVE","FIFO"};
const char* visFilterModes[4]={"LP","HP","BP","NOTCH"};

ImU32 withAlpha(ImU32 c, float a) {
  if (a<0.0f) a=0.0f;
  if (a>1.0f) a=1.0f;
  return (c&~IM_COL32_A_MASK)|((ImU32)(((c>>IM_COL32_A_SHIFT)&0xff)*a)<<IM_COL32_A_SHIFT);
}

ImU32 mixCol(ImU32 a, ImU32 b, float t) {
  ImVec4 va=ImGui::ColorConvertU32ToFloat4(a);
  ImVec4 vb=ImGui::ColorConvertU32ToFloat4(b);
  return ImGui::ColorConvertFloat4ToU32(ImLerp(va,vb,t));
}

void visNoteName(char* out, size_t len, double hz) {
  static const char* names[12]={"C-","C#","D-","D#","E-","F-","F#","G-","G#","A-","A#","B-"};
  if (hz<=0.0) {
    snprintf(out,len,"---");
    return;
  }
  int note=(int)floor(12.0*log2(hz/440.0)+57.5);
  if (note<0 || note>=12*10) {
    snprintf(out,len,"???");
    return;
  }
  snprintf(out,len,"%s%d",names[note%12],note/12);
}

// Peak-to-peak of the last 30 ms of a channel scope, 0..1.
float scopeLevel(DivDispatchOscBuffer* buf) {
  if (buf==NULL) return 0.0f;
  const int n=65536*0.03;
  unsigned short needle=buf->needle>>16;
  short lo=32767, hi=-32768;
  for (unsigned short i=needle-n; i!=needle; i++) {
    short y=buf->data[i];
    if (y==-1) continue;
    if (y<lo) lo=y;
    if (y>hi) hi=y;
  }
  if (hi<lo) return 0.0f;
  return MIN(1.0f,(float)(hi-lo)/65536.0f);
}

// Min/max overview of a PCM RAM range, rebuilt when the range or width changes.
struct PcmOverview {
  unsigned int start=0, end=0;
  int cols=0;
  double builtAt=-100.0;
  std::vector<unsigned char> lo, hi;
};

struct VisCache {
  const void* chip=NULL;
  PcmOverview pcm[4];
  unsigned int lastTrigger[12]={};
  double triggerAt[12]={};
  float level[12]={};
};

std::vector<VisCache> visCaches;

}  // namespace

static void visTextV(ImDrawList* dl, ImFont* font, float size, ImVec2 pos, ImU32 col, const char* fmt, va_list args) {
  char buf[256];
  vsnprintf(buf,sizeof(buf),fmt,args);
  dl->AddText(font,size,pos,col,buf);
}

static void visPrint(ImDrawList* dl, ImFont* font, float size, ImVec2 pos, ImU32 col, const char* fmt, ...) {
  va_list args;
  va_start(args,fmt);
  visTextV(dl,font,size,pos,col,fmt,args);
  va_end(args);
}

static float visTextWidth(ImFont* font, float size, const char* text) {
  return font->CalcTextSizeA(size,FLT_MAX,0.0f,text).x;
}

static void visBar(ImDrawList* dl, ImVec2 a, ImVec2 b, float frac, ImU32 col, float back=-1.0f) {
  dl->AddRectFilled(a,b,visBarBg);
  if (back>0.0f) {
    dl->AddRectFilled(a,ImVec2(a.x+(b.x-a.x)*MIN(back,1.0f),b.y),withAlpha(col,0.25f));
  }
  if (frac>0.0f) {
    dl->AddRectFilled(a,ImVec2(a.x+(b.x-a.x)*MIN(frac,1.0f),b.y),col);
  }
}

// Stereo pan meter: two bars, left and right gain in fifteenths.
static void visPan(ImDrawList* dl, ImFont* font, float fs, ImVec2 a, ImVec2 b, int left, int right, ImU32 col) {
  float mid=(a.x+b.x)*0.5f;
  float gap=2.0f;
  dl->AddRectFilled(a,b,visBarBg);
  dl->AddRectFilled(ImVec2(mid-gap-(mid-gap-a.x)*left/15.0f,a.y),ImVec2(mid-gap,b.y),col);
  dl->AddRectFilled(ImVec2(mid+gap,a.y),ImVec2(mid+gap+(b.x-mid-gap)*right/15.0f,b.y),col);
  char label[3]={left?'L':'-',right?'R':'-',0};
  dl->AddText(font,fs,ImVec2(b.x+fs*0.4f,a.y+(b.y-a.y-fs)*0.5f),visDim,label);
}

// Edge-triggered scope of the last ~25 ms of a channel's output.
static void visScope(ImDrawList* dl, ImVec2 a, ImVec2 b, DivDispatchOscBuffer* buf, ImU32 col, float gain) {
  float midY=(a.y+b.y)*0.5f;
  dl->AddLine(ImVec2(a.x,midY),ImVec2(b.x,midY),visGrid);
  if (buf==NULL) return;
  const int n=65536*0.025;
  unsigned short needle=buf->needle>>16;
  unsigned short start=needle-n;
  // Look one window further back for a rising zero crossing.
  short prev=0;
  for (int i=0; i<n; i++) {
    unsigned short at=needle-2*n+i;
    short y=buf->data[at];
    if (y==-1) continue;
    if (prev<0 && y>=0) {
      start=at;
      break;
    }
    prev=y;
  }
  int cols=MAX(2,(int)(b.x-a.x));
  static std::vector<ImVec2> points;
  points.resize(cols);
  short last=0;
  for (int c=0; c<cols; c++) {
    unsigned short at=start+(unsigned short)((long long)c*n/cols);
    short y=buf->data[at];
    if (y!=-1) last=y;
    float v=MAX(-1.0f,MIN(1.0f,(float)last/32768.0f*gain));
    points[c]=ImVec2(a.x+c,midY-v*(b.y-a.y)*0.5f);
  }
  dl->PushClipRect(a,b,true);
  dl->AddPolyline(points.data(),cols,col,0,1.2f);
  dl->PopClipRect();
}

// Horizontal A0..C8 keyboard scale with a marker at `hz`.
static void visPitchBar(ImDrawList* dl, ImFont* font, float fs, ImVec2 a, ImVec2 b, double hz, ImU32 col, float alpha) {
  dl->AddRectFilled(a,b,visPanel);
  dl->AddRect(a,b,visGrid);
  const double lo=log2(27.5), hi=log2(4186.009);
  // Octave ticks at every C.
  for (int oct=1; oct<=8; oct++) {
    double c=log2(440.0*pow(2.0,(oct*12-57)/12.0));
    float x=a.x+(float)((c-lo)/(hi-lo))*(b.x-a.x);
    dl->AddLine(ImVec2(x,a.y),ImVec2(x,b.y),visGrid);
  }
  float labelY=b.y-fs-1.0f;
  dl->AddText(font,fs*0.85f,ImVec2(a.x+3.0f,labelY),visDim,"A0");
  float a4=a.x+(float)((log2(440.0)-lo)/(hi-lo))*(b.x-a.x);
  dl->AddText(font,fs*0.85f,ImVec2(a4+2.0f,labelY),visDim,"A4");
  dl->AddText(font,fs*0.85f,ImVec2(b.x-visTextWidth(font,fs*0.85f,"C8")-3.0f,labelY),visDim,"C8");
  if (hz>0.0) {
    float t=(float)((log2(hz)-lo)/(hi-lo));
    if (t>=0.0f && t<=1.0f) {
      float x=a.x+t*(b.x-a.x);
      float w=MAX(3.0f,fs*0.35f);
      dl->AddRectFilled(ImVec2(x-w*0.5f,a.y+1.0f),ImVec2(x+w*0.5f,b.y-1.0f),withAlpha(col,0.35f+0.65f*alpha));
    }
  }
}

static void visPanelFrame(ImDrawList* dl, ImVec2 a, ImVec2 b, ImU32 col, float active) {
  dl->AddRectFilled(a,b,visPanel);
  dl->AddRect(a,b,mixCol(visGrid,col,0.35f+0.65f*active),0.0f,0,1.0f+active);
}

void FurnaceGUI::drawS3HSVis() {
  if (nextWindow==GUI_WINDOW_S3HS_VIS) {
    s3hsVisOpen=true;
    ImGui::SetNextWindowFocus();
    nextWindow=GUI_WINDOW_NOTHING;
  }
  if (!s3hsVisOpen) {
    visCaches.clear();
    return;
  }
  ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding,ImVec2(0,0));
  if (ImGui::Begin("3HS88PWN4 Visualizer",&s3hsVisOpen,globalWinFlags|ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoScrollWithMouse,_("3HS88PWN4 Visualizer"))) {
    std::vector<int> chips;
    for (int i=0; i<e->song.systemLen; i++) {
      if (e->song.system[i]==DIV_SYSTEM_S3HS) chips.push_back(i);
    }
    if (chips.empty()) {
      ImGui::SetCursorPos(ImVec2(8.0f*dpiScale,8.0f*dpiScale));
      ImGui::TextUnformatted(_("this song has no 3HS88PWN4."));
    } else {
      if (visCaches.size()!=chips.size()) visCaches.resize(chips.size());
      ImVec2 avail=ImGui::GetContentRegionAvail();
      // Several chips share the window vertically.
      float chipH=MAX(200.0f*dpiScale,avail.y/chips.size());
      for (size_t k=0; k<chips.size(); k++) {
        drawS3HSVisChip(chips[k],ImVec2(avail.x,chipH),k);
      }
      if (e->isPlaying()) {
        WAKE_UP;
      }
    }
  }
  if (ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows)) curWindow=GUI_WINDOW_S3HS_VIS;
  ImGui::End();
  ImGui::PopStyleVar();
}

void FurnaceGUI::drawS3HSVisChip(int sysIndex, ImVec2 size, size_t cacheIndex) {
  DivPlatformS3HS* dispatch=(DivPlatformS3HS*)e->getDispatch(sysIndex);
  ImVec2 origin=ImGui::GetCursorScreenPos();
  ImGui::Dummy(size);
  if (dispatch==NULL || dispatch->getChip()==NULL) return;
  const hs88::Engine* chip=dispatch->getChip();
  VisCache& cache=visCaches[cacheIndex];
  if (cache.chip!=chip) {
    cache=VisCache();
    cache.chip=chip;
    for (int i=0; i<12; i++) cache.lastTrigger[i]=dispatch->getVisTrigger(i);
  }
  const double now=ImGui::GetTime();
  for (int i=0; i<12; i++) {
    unsigned int t=dispatch->getVisTrigger(i);
    if (t!=cache.lastTrigger[i]) {
      cache.lastTrigger[i]=t;
      cache.triggerAt[i]=now;
    }
    float lv=scopeLevel(dispatch->getOscBuffer(i));
    cache.level[i]=MAX(cache.level[i]*0.85f,MIN(1.0f,sqrtf(lv)*1.4f));
  }

  ImDrawList* dl=ImGui::GetWindowDrawList();
  ImFont* font=patFont;
  ImFont* bigF=headFont?headFont:patFont;
  const float W=size.x, H=size.y;
  const float fs=MAX(9.0f*dpiScale,MIN(W/150.0f,H/82.0f));
  const float pad=MAX(3.0f,fs*0.45f);
  const ImVec2 end=origin+size;
  dl->AddRectFilled(origin,end,visBg);
  dl->PushClipRect(origin,end,true);

  float y=origin.y+pad;
  // Title: author - "song"
  {
    String title;
    if (!e->song.author.empty()) title+=e->song.author+" - ";
    title+="\""+(e->song.name.empty()?String("untitled"):e->song.name)+"\"";
    dl->AddRect(ImVec2(origin.x+pad,y),ImVec2(end.x-pad,y+fs*1.6f),visGrid);
    dl->AddText(font,fs,ImVec2(origin.x+pad*2.0f,y+fs*0.3f),visText,title.c_str());
    y+=fs*1.6f+pad;
  }
  // Chip header
  {
    dl->AddRect(ImVec2(origin.x+pad,y),ImVec2(end.x-pad,end.y-pad),IM_COL32(150,90,40,255));
    float x=origin.x+pad*2.0f;
    dl->AddText(bigF,fs*1.8f,ImVec2(x,y+fs*0.2f),visText,"3HS88PWN4");
    x+=visTextWidth(bigF,fs*1.8f,"3HS88PWN4")+fs*3.0f;
    visPrint(dl,font,fs,ImVec2(x,y+fs*0.7f),visDim,"48000 Hz  :  8 harmonic x 8-op + 4 PCM/WT/noise");
    y+=fs*2.4f;
  }

  const float innerX0=origin.x+pad*2.0f, innerX1=end.x-pad*2.0f;
  const float innerW=innerX1-innerX0;
  auto sectionLabel=[&](const char* name, const char* sub) {
    dl->AddRect(ImVec2(innerX0,y),ImVec2(innerX1,y+fs*1.3f),visGrid);
    dl->AddText(font,fs,ImVec2(innerX0+pad,y+fs*0.15f),visText,name);
    dl->AddText(font,fs*0.9f,ImVec2(innerX0+pad+fs*9.0f,y+fs*0.2f),visDim,sub);
    y+=fs*1.3f+pad;
  };

  const float bottom=end.y-pad*2.0f;
  const float remaining=bottom-y-2.0f*(fs*1.3f+pad);
  const float hsH=remaining*0.62f;
  const float pcmH=remaining-hsH;

  // ---- harmonic channels: 2 columns x 4 rows ----
  sectionLabel("HARMONIC","8 channels, 8 operators each");
  {
    const float cellW=(innerW-pad)*0.5f;
    const float cellH=(hsH-pad*3.0f)/4.0f;
    for (int ch=0; ch<8; ch++) {
      const hs816::ChannelParams& p=chip->harmonic(ch).params();
      const bool gate=chip->harmonic(ch).gateOn();
      const ImU32 col=visChanColors[ch];
      ImVec2 a(innerX0+(ch&1)*(cellW+pad),y+(ch>>1)*(cellH+pad));
      ImVec2 b=a+ImVec2(cellW,cellH);
      const float act=p.muted?0.0f:cache.level[ch];
      visPanelFrame(dl,a,b,col,act);
      dl->PushClipRect(a,b,true);

      // Left column: name and scope.
      const float lw=MIN(fs*8.0f,cellW*0.18f);
      char name[8];
      snprintf(name,sizeof(name),"HS%d",ch+1);
      dl->AddText(bigF,fs*1.5f,ImVec2(a.x+pad,a.y+pad*0.6f),p.muted?visDim:visText,name);
      dl->AddText(font,fs*0.85f,ImVec2(a.x+pad,a.y+pad*0.6f+fs*1.6f),visDim,p.muted?"muted":(gate?"gate on":"gate off"));
      visScope(dl,ImVec2(a.x+pad,a.y+pad+fs*2.8f),ImVec2(a.x+lw,b.y-pad),dispatch->getOscBuffer(ch),col,1.5f);

      // Header row: pitch, algorithm, feedback, pan.
      const float rx0=a.x+lw+pad, rx1=b.x-pad;
      const float rowH=fs*2.2f;
      const float ry=a.y+pad*0.6f;
      const float pitchW=(rx1-rx0)*0.58f;
      visPitchBar(dl,font,fs,ImVec2(rx0,ry),ImVec2(rx0+pitchW,ry+rowH),p.frequency,col,gate?1.0f:act);
      char note[8];
      visNoteName(note,sizeof(note),p.frequency);
      visPrint(dl,font,fs,ImVec2(rx0+3.0f,ry+1.0f),visText,"%5d Hz %s",p.frequency,note);
      const float infoX=rx0+pitchW+pad;
      const int alg=p.algorithm;
      const unsigned char carriers=(alg<13)?visAlgCarriers[alg]:0;
      visPrint(dl,font,fs,ImVec2(infoX,ry),visText,"AL %d %s",alg,(alg<13)?visAlgNames[alg]:"(none)");
      if (p.feedback!=128 && p.feedbackTarget<8) {
        visPrint(dl,font,fs,ImVec2(infoX,ry+fs*1.1f),visDim,"FB %+d > OP%d",(int)p.feedback-128,p.feedbackTarget+1);
      } else {
        visPrint(dl,font,fs,ImVec2(infoX,ry+fs*1.1f),visDim,"FB 0");
      }
      const float panW=MIN(fs*6.0f,(rx1-infoX)*0.3f);
      const int panL=(p.panLeft==0 && p.panRight==0)?15:p.panLeft;
      const int panR=(p.panLeft==0 && p.panRight==0)?15:p.panRight;
      visPan(dl,font,fs,ImVec2(rx1-panW-fs*1.6f,ry+fs*0.2f),ImVec2(rx1-fs*1.6f,ry+fs*0.9f),panL,panR,col);

      // Operator cells: one row of 8 when there is room, otherwise 4x2.
      const float opY0=ry+rowH+pad;
      const float opAreaW=rx1-rx0, opAreaH=b.y-pad-opY0;
      const int perRow=(opAreaW/8.0f>=fs*7.5f || opAreaH<fs*7.0f)?8:4;
      const int rows=8/perRow;
      const float opGap=MAX(2.0f,pad*0.6f);
      const float opW=(opAreaW-opGap*(perRow-1))/perRow;
      const float opH=(opAreaH-opGap*(rows-1))/rows;
      for (int op=0; op<8; op++) {
        const hs816::OperatorParams& o=p.operators[op];
        ImVec2 oa(rx0+(op%perRow)*(opW+opGap),opY0+(op/perRow)*(opH+opGap));
        ImVec2 ob=oa+ImVec2(opW,opH);
        const bool carrier=carriers&(1<<op);
        int stage=dispatch->getVisStage(ch,op);
        const int level=dispatch->getVisLevel(ch,op);
        if (stage==3 && level==0) stage=4;
        const float amp=(float)((level*o.volume)>>16)/255.0f;
        dl->AddRect(oa,ob,carrier?withAlpha(col,0.85f):visGrid);
        dl->PushClipRect(oa,ob,true);
        const float tx=oa.x+3.0f;
        float ty=oa.y+2.0f;
        visPrint(dl,font,fs*0.9f,ImVec2(tx,ty),carrier?visText:visDim,"OP%d",op+1);
        dl->AddText(font,fs*0.9f,ImVec2(tx+visTextWidth(font,fs*0.9f,"OP0 "),ty),visStageColors[stage],visStageNames[stage]);
        if (op==p.feedbackTarget && p.feedback!=128) {
          dl->AddText(font,fs*0.8f,ImVec2(ob.x-visTextWidth(font,fs*0.8f,"FB")-3.0f,ty),visDim,"FB");
        }
        ty+=fs*1.05f;
        const float barH=MAX(3.0f,fs*0.45f);
        visBar(dl,ImVec2(tx,ty),ImVec2(ob.x-3.0f,ty+barH),amp,visStageColors[stage==4?2:stage],o.volume/255.0f);
        ty+=barH+2.0f;
        const double ratio=(op==0)?1.0:o.ratio/4096.0;
        if (opH>=fs*4.2f) {
          visPrint(dl,font,fs*0.85f,ImVec2(tx,ty),visDim,"x%.3f V%d",ratio,o.volume);
          ty+=fs*0.95f;
          visPrint(dl,font,fs*0.85f,ImVec2(tx,ty),visDim,"%.1fHz",p.frequency*ratio);
          ty+=fs*0.95f;
        } else {
          visPrint(dl,font,fs*0.85f,ImVec2(tx,ty),visDim,"x%.2f V%d",ratio,o.volume);
          ty+=fs*0.95f;
        }
        // Waveform thumbnail.
        if (ob.y-2.0f-ty>=fs*0.8f) {
          ImVec2 wa(tx,ty+1.0f), wb(ob.x-3.0f,ob.y-2.0f);
          float midY=(wa.y+wb.y)*0.5f;
          ImVec2 pts[65];
          for (int i=0; i<=64; i++) {
            int s=chip->waves().sample(o.wave,(i*4)&255);
            pts[i]=ImVec2(wa.x+(wb.x-wa.x)*i/64.0f,midY-(s/128.0f)*(wb.y-wa.y)*0.5f);
          }
          dl->AddPolyline(pts,65,withAlpha(carrier?col:visDim,0.4f+0.6f*MAX(amp,carrier?0.25f:0.0f)),0,1.0f);
          char w[4];
          snprintf(w,sizeof(w),"%X",o.wave);
          dl->AddText(font,fs*0.75f,ImVec2(wb.x-visTextWidth(font,fs*0.75f,w),wb.y-fs*0.8f),visDim,w);
        }
        dl->PopClipRect();
      }
      dl->PopClipRect();
    }
    y+=hsH;
  }

  // ---- PCM/WT channels: 2 x 2 ----
  sectionLabel("PCM / WT","4 lanes: sample, wavetable, noise, FIFO");
  {
    const float cellW=(innerW-pad)*0.5f;
    const float cellH=(pcmH-pad)*0.5f;
    for (int i=0; i<4; i++) {
      const int ch=8+i;
      const hs88::PcmChannelParams& p=chip->pcm(i).params();
      const ImU32 col=visChanColors[ch];
      ImVec2 a(innerX0+(i&1)*(cellW+pad),y+(i>>1)*(cellH+pad));
      ImVec2 b=a+ImVec2(cellW,cellH);
      const float act=p.muted?0.0f:cache.level[ch];
      visPanelFrame(dl,a,b,col,act);
      dl->PushClipRect(a,b,true);
      const float tx=a.x+pad;
      float ty=a.y+pad*0.6f;
      const int mode=p.mode;
      const char* modeName=(mode<6)?visPcmModes[mode]:"(off)";
      visPrint(dl,font,fs,ImVec2(tx,ty),visText,"PCM %d",i+1);
      dl->AddText(font,fs,ImVec2(tx+fs*5.0f,ty),col,modeName);

      // State badge.
      const char* badge;
      ImU32 badgeCol;
      const float sinceTrigger=(float)(now-cache.triggerAt[ch]);
      if (p.muted) {
        badge="MUTE"; badgeCol=visDim;
      } else if (cache.triggerAt[ch]>0.0 && sinceTrigger<0.25f) {
        badge="TRIGGER"; badgeCol=col;
      } else if (p.volume>0 && p.frequency>0 && act>0.02f) {
        badge="PLAY"; badgeCol=mixCol(col,visText,0.3f);
      } else {
        badge="IDLE"; badgeCol=visDim;
      }
      {
        float bw=visTextWidth(font,fs*0.9f,badge)+fs*0.8f;
        ImVec2 ba(b.x-pad-bw,ty), bb(b.x-pad,ty+fs*1.1f);
        if (badgeCol==visDim) {
          dl->AddRect(ba,bb,visDim);
          dl->AddText(font,fs*0.9f,ba+ImVec2(fs*0.4f,fs*0.08f),visDim,badge);
        } else {
          dl->AddRectFilled(ba,bb,badgeCol);
          dl->AddText(font,fs*0.9f,ba+ImVec2(fs*0.4f,fs*0.08f),visBg,badge);
        }
      }
      ty+=fs*1.3f;
      // Volume and pan.
      const float panW=MIN(fs*6.0f,cellW*0.2f);
      visBar(dl,ImVec2(tx,ty),ImVec2(b.x-pad-panW-fs*2.5f,ty+fs*0.55f),p.volume/255.0f,col);
      const int panL=(p.panLeft==0 && p.panRight==0)?15:p.panLeft;
      const int panR=(p.panLeft==0 && p.panRight==0)?15:p.panRight;
      visPan(dl,font,fs,ImVec2(b.x-pad-panW-fs*1.6f,ty),ImVec2(b.x-pad-fs*1.6f,ty+fs*0.55f),panL,panR,col);
      ty+=fs*0.55f+pad;

      // Main display.
      ImVec2 da(tx,ty), db(b.x-pad,b.y-pad-fs*1.3f);
      dl->AddRectFilled(da,db,visBg);
      dl->AddRect(da,db,visGrid);
      const float midY=(da.y+db.y)*0.5f;
      const unsigned long long pos=dispatch->getVisPcmPos(i);
      char footer[160];
      const char* filt="";
      char filtBuf[40];
      if (p.cutoff!=0) {
        snprintf(filtBuf,sizeof(filtBuf),"  %s c%d r%d",(p.filterMode<4)?visFilterModes[p.filterMode]:"LP*",p.cutoff,p.resonance);
        filt=filtBuf;
      }
      if (mode==0) {
        // Sample in PCM RAM: start, end (inclusive) and loop as 24-bit BE.
        const unsigned int start=(p.table[0]<<16)|(p.table[1]<<8)|p.table[2];
        const unsigned int endAddr=(p.table[3]<<16)|(p.table[4]<<8)|p.table[5];
        const unsigned int loop=(p.table[6]<<16)|(p.table[7]<<8)|p.table[8];
        const bool valid=endAddr>start && endAddr<hs88::kPcmMemorySize;
        if (valid) {
          PcmOverview& ov=cache.pcm[i];
          const int cols=MAX(2,(int)(db.x-da.x)-2);
          if (ov.start!=start || ov.end!=endAddr || ov.cols!=cols || now-ov.builtAt>2.0) {
            ov.start=start; ov.end=endAddr; ov.cols=cols; ov.builtAt=now;
            ov.lo.assign(cols,128); ov.hi.assign(cols,128);
            const unsigned long long len=endAddr-start+1;
            for (int c=0; c<cols; c++) {
              unsigned long long s0=start+len*c/cols, s1=start+len*(c+1)/cols;
              if (s1<=s0) s1=s0+1;
              // Sample at most 256 points per column.
              unsigned long long step=MAX(1ULL,(s1-s0)/256);
              unsigned char lo=255, hi=0;
              for (unsigned long long s=s0; s<s1; s+=step) {
                unsigned char v=chip->readPcmMemory((unsigned int)s);
                if (v<lo) lo=v;
                if (v>hi) hi=v;
              }
              ov.lo[c]=lo; ov.hi[c]=hi;
            }
          }
          const float hh=(db.y-da.y)*0.5f-2.0f;
          const bool looped=loop<endAddr && loop!=0xffffff;
          // Cursor address, following hs88::PcmChannel::readPcm.
          unsigned long long addr;
          if (start+pos>endAddr && looped) {
            addr=start+pos%(endAddr-loop);
          } else {
            addr=MIN((unsigned long long)start+pos,(unsigned long long)endAddr);
          }
          const float span=(float)(endAddr-start);
          const float cx=da.x+1.0f+(float)(addr-start)/span*(cols-1);
          if (looped && loop>=start) {
            float lx=da.x+1.0f+(float)(loop-start)/span*(cols-1);
            dl->AddRectFilled(ImVec2(lx,da.y+1.0f),ImVec2(db.x-1.0f,db.y-1.0f),withAlpha(col,0.07f));
            dl->AddLine(ImVec2(lx,da.y+1.0f),ImVec2(lx,db.y-1.0f),withAlpha(col,0.5f));
          }
          for (int c=0; c<cols; c++) {
            float x=da.x+1.0f+c;
            const bool played=x<=cx;
            ImU32 lc=played?withAlpha(col,0.55f):col;
            dl->AddLine(ImVec2(x,midY-(ov.hi[c]-128)/128.0f*hh),ImVec2(x,midY-(ov.lo[c]-128)/128.0f*hh+1.0f),lc);
          }
          const bool ended=!looped && start+pos>=endAddr;
          if (!ended) {
            float cw=MAX(2.0f,fs*0.5f);
            dl->AddRectFilled(ImVec2(cx-cw*0.5f,da.y+1.0f),ImVec2(cx+cw*0.5f,db.y-1.0f),withAlpha(visText,0.18f));
            dl->AddLine(ImVec2(cx,da.y+1.0f),ImVec2(cx,db.y-1.0f),withAlpha(visText,0.7f));
          }
          snprintf(footer,sizeof(footer),"@%06X +%X%s  F %d Hz  V %d%s  %llu/%u",start,endAddr-start,looped?" LOOP":"",p.frequency,p.volume,filt,
            (unsigned long long)(MIN((unsigned long long)addr,(unsigned long long)endAddr)-start),endAddr-start);
        } else {
          dl->AddText(font,fs,ImVec2(da.x+pad,midY-fs*0.5f),visDim,"no sample range");
          snprintf(footer,sizeof(footer),"F %d Hz  V %d%s",p.frequency,p.volume,filt);
        }
      } else if (mode==1 || mode==4) {
        // 32-step wavetable with the cursor.
        const int cur=(int)(pos%32);
        const float sw=(db.x-da.x-2.0f)/32.0f;
        const float hh=(db.y-da.y)*0.5f-2.0f;
        ImVec2 pts[33];
        for (int s=0; s<32; s++) {
          float x0=da.x+1.0f+s*sw;
          float v=(p.table[s]-128)/128.0f;
          ImU32 bc=(s==cur)?visText:withAlpha(col,0.3f+0.4f*act);
          dl->AddRectFilled(ImVec2(x0+1.0f,MIN(midY,midY-v*hh)),ImVec2(x0+sw-1.0f,MAX(midY,midY-v*hh)+1.0f),bc);
          pts[s]=ImVec2(x0+sw*0.5f,midY-v*hh);
        }
        pts[32]=ImVec2(db.x-1.0f,midY-(p.table[0]-128)/128.0f*hh);
        dl->AddPolyline(pts,33,col,0,1.5f);
        dl->AddLine(ImVec2(da.x,midY),ImVec2(db.x,midY),visGrid);
        char note[8];
        visNoteName(note,sizeof(note),p.frequency);
        snprintf(footer,sizeof(footer),"F %d Hz %s  V %d%s  STEP %2d/32",p.frequency,note,p.volume,filt,cur);
      } else {
        // Noise and FIFO: live scope.
        visScope(dl,ImVec2(da.x+1.0f,da.y+1.0f),ImVec2(db.x-1.0f,db.y-1.0f),dispatch->getOscBuffer(ch),col,1.5f);
        if (mode==5) {
          snprintf(footer,sizeof(footer),"F %d Hz  V %d%s  FIFO %d",p.frequency,p.volume,filt,(int)chip->pcm(i).fifoLength());
        } else {
          snprintf(footer,sizeof(footer),"F %d Hz  V %d%s  PERIOD %s",p.frequency,p.volume,filt,mode==3?"64":"65536");
        }
      }
      dl->AddText(font,fs*0.9f,ImVec2(tx,b.y-pad-fs*1.05f),visDim,footer);
      dl->PopClipRect();
    }
  }
  dl->PopClipRect();
}
