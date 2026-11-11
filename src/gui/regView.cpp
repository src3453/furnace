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

#include "../engine/bsr.h"
#include "gui.h"
#include "IconsFontAwesome4.h"
#include <imgui.h>

void FurnaceGUI::drawRegView() {
  if (nextWindow==GUI_WINDOW_REGISTER_VIEW) {
    regViewOpen=true;
    ImGui::SetNextWindowFocus();
    nextWindow=GUI_WINDOW_NOTHING;
  }
  if (!regViewOpen) {
    regViewState.clear();
    return;
  }
  if (ImGui::Begin("Register View",&regViewOpen,globalWinFlags,_("Register View"))) {
    const double now=ImGui::GetTime();
    const double fadeDuration=0.5;
    regViewState.resize(e->song.systemLen);
    ImVec2 prevPos=ImGui::GetCursorPos();
    ImVec2 topPos=ImGui::GetCursorPos();
    topPos.x+=ImGui::GetContentRegionAvail().x-ImGui::CalcTextSize(ICON_FA_BARS).x;
    topPos.y+=ImGui::GetScrollY();
    if (ImGui::IsWindowHovered()) {
      ImGui::SetCursorPos(topPos);
      ImGui::TextUnformatted(ICON_FA_BARS "##regViewSettings");
      if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
        ImGui::SetTooltip(_("Register View settings"));
      }
      if (ImGui::IsItemClicked()) {
        ImGui::OpenPopup("regViewSettingsPopup");
      }
      ImGui::SetCursorPos(prevPos);
    }
    if (ImGui::BeginPopup("regViewSettingsPopup")) {
      if (ImGui::InputInt(_("Bytes per column##RegViewColumns"),&regViewColumns,1,4)) {
        if (regViewColumns<1) regViewColumns=1;
        if (regViewColumns>64) regViewColumns=64;
      }
      ImGui::Checkbox(_("Highlight register changes"),&regViewHighlight);
      ImGui::BeginDisabled(!regViewHighlight);
      ImGui::ColorEdit3(_("Highlight color"),(float*)&regViewHighlightColor,ImGuiColorEditFlags_DisplayHex);
      ImGui::EndDisabled();
      ImGui::EndPopup();
    }
    for (int i=0; i<e->song.systemLen; i++) {
      ImGui::Text("%d. %s",i+1,getSystemName(e->song.system[i]));
      int size=0;
      int depth=8;
      unsigned char* regPool=e->getRegisterPool(i,size,depth);
      unsigned short* regPoolW=(unsigned short*)regPool;
      RegViewState& state=regViewState[i];
      if (regPool==NULL) {
        state=RegViewState();
        ImGui::Text(_("- no register pool available"));
      } else {
        // A new pool establishes a baseline without flashing every register.
        const bool reset=state.pool!=regPool || state.system!=e->song.system[i] ||
          state.depth!=depth || state.values.size()!=(size_t)size;
        if (reset) {
          state.pool=regPool;
          state.system=e->song.system[i];
          state.depth=depth;
          state.values.resize(size);
          state.changedAt.assign(size,-fadeDuration);
        }
        // Update all entries even when their table cells are outside the viewport.
        for (int j=0; j<size; j++) {
          unsigned short value=depth==16?regPoolW[j]:regPool[j];
          if (!reset && state.values[j]!=value) state.changedAt[j]=now;
          state.values[j]=value;
          if (regViewHighlight && now-state.changedAt[j]<fadeDuration) {
            WAKE_UP;
          }
        }
        ImGui::PushFont(patFont);
        int bytesPerValue=MAX(1,depth/8);
        int columns=MAX(1,regViewColumns/bytesPerValue);
        if (ImGui::BeginTable("Memory",1+columns)) {
          String addrFmt=fmt::sprintf("%X",MAX(size-1,0));
          float widthOne=MAX(ImGui::CalcTextSize(addrFmt.c_str()).x,ImGui::CalcTextSize(bytesPerValue==2?"0000":"00").x);
          ImGui::TableSetupColumn("addr",ImGuiTableColumnFlags_WidthFixed,widthOne);

          ImGui::TableNextRow();
          ImGui::TableNextColumn();
          for (int j=0; j<columns; j++) {
            ImGui::TableNextColumn();
            ImGui::TextColored(uiColors[GUI_COLOR_PATTERN_ROW_INDEX]," %X",j);
          }
          for (int row=0; row<size; row+=columns) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextColored(uiColors[GUI_COLOR_PATTERN_ROW_INDEX],"%.2X",row);
            for (int col=0; col<columns; col++) {
              ImGui::TableNextColumn();
              if (row+col>=size) continue;
              const int index=row+col;
              const float highlight=regViewHighlight?(float)MAX(0.0,1.0-(now-state.changedAt[index])/fadeDuration):0.0f;
              const ImVec4 normalColor=ImGui::GetStyleColorVec4(ImGuiCol_Text);
              const ImVec4 color=ImLerp(normalColor,ImVec4(regViewHighlightColor.x,regViewHighlightColor.y,regViewHighlightColor.z,normalColor.w),highlight);
              switch (depth) {
                case 8: ImGui::TextColored(color,"%.2x",state.values[index]); break;
                case 16: ImGui::TextColored(color,"%.4x",state.values[index]); break;
                default: ImGui::Text("??"); break;
              }
            }
          }

          ImGui::EndTable();
        }
        ImGui::PopFont();
      }
    }
  }
  if (ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows)) curWindow=GUI_WINDOW_REGISTER_VIEW;
  ImGui::End();
}
