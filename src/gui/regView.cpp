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

#include "gui.h"
#include <imgui.h>

void FurnaceGUI::drawRegView() {
  if (nextWindow==GUI_WINDOW_REGISTER_VIEW) {
    channelsOpen=true;
    ImGui::SetNextWindowFocus();
    nextWindow=GUI_WINDOW_NOTHING;
  }
  if (!regViewOpen) return;
  if (ImGui::Begin("Register View",&regViewOpen,globalWinFlags,_("Register View"))) {
    for (int i=0; i<e->song.systemLen; i++) {
      ImGui::Text("%d. %s",i+1,getSystemName(e->song.system[i]));
      int size=0;
      int depth=8;
      unsigned char* regPool=e->getRegisterPool(i,size,depth);
      unsigned short* regPoolW=(unsigned short*)regPool;
      if (regPool==NULL) {
        ImGui::Text(_("- no register pool available"));
      } else {
        ImGui::PushFont(patFont);
        if (depth==16) {
          if (ImGui::BeginTable("Memory",9)) {
            float widthOne=ImGui::CalcTextSize("0").x;
            if (size>0xfff) {
              ImGui::TableSetupColumn("addr",ImGuiTableColumnFlags_WidthFixed,widthOne*4.0f);
            } else if (size>0xff) {
              ImGui::TableSetupColumn("addr",ImGuiTableColumnFlags_WidthFixed,widthOne*3.0f);
            } else {
              ImGui::TableSetupColumn("addr",ImGuiTableColumnFlags_WidthFixed,widthOne*2.0f);
            }

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            for (int j=0; j<8; j++) {
              ImGui::TableNextColumn();
              ImGui::TextColored(uiColors[GUI_COLOR_PATTERN_ROW_INDEX]," %X",j);
            }

            int rows=(size+7)>>3;
            for (int row=0; row<rows; row++) {
              ImGui::TableNextRow();
              ImGui::TableNextColumn();
              ImGui::TextColored(uiColors[GUI_COLOR_PATTERN_ROW_INDEX],"%.2X",row*8);
              for (int col=0; col<8; col++) {
                int idx=row*8+col;
                ImGui::TableNextColumn();
                if (idx>=size) continue;
                ImGui::Text("%.4x",regPoolW[idx]);
              }
            }

            ImGui::EndTable();
          }
        } else {
          if (ImGui::BeginTable("Memory",17)) {
            ImGui::TableSetupColumn("addr",ImGuiTableColumnFlags_WidthFixed);

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            for (int j=0; j<16; j++) {
              ImGui::TableNextColumn();
              ImGui::TextColored(uiColors[GUI_COLOR_PATTERN_ROW_INDEX]," %X",j);
            }

            int rows=(size+15)>>4;
            for (int row=0; row<rows; row++) {
              ImGui::TableNextRow();
              ImGui::TableNextColumn();
              ImGui::TextColored(uiColors[GUI_COLOR_PATTERN_ROW_INDEX],"%.2X",row*16);
              for (int col=0; col<16; col++) {
                int idx=row*16+col;
                ImGui::TableNextColumn();
                if (idx>=size) continue;
                ImGui::Text("%.2x",regPool[idx]);
              }
            }

            ImGui::EndTable();
          }
        }
        ImGui::PopFont();
      }
    }
  }
  if (ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows)) curWindow=GUI_WINDOW_REGISTER_VIEW;
  ImGui::End();
}
