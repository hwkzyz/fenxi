function Assessment=Run_04_SummarizePaperValue()
root=fileparts(mfilename('fullpath'));
A=readtable(fullfile(root,'output','01_energy_matched_phase_aperture','gate.csv'),TextType='string');
B=readtable(fullfile(root,'output','02_conditional_crlb_monte_carlo','summary.csv'),TextType='string');
C=readtable(fullfile(root,'output','03_x01_model_distance_bias_projection','summary.csv'),TextType='string');
gateA=all(A.gate_pass);gateB=all(B.all_same_order(B.snr_db>=15));gateC=C.full_margin_flip(1);
Assessment=table(gateA,gateB,gateC,sum([gateA gateB gateC]),...
    'VariableNames',{'phase_aperture_gate','crlb_gate_high_mid_snr','x01_margin_flip_gate','passed_gate_count'});
out=fullfile(root,'output');writetable(Assessment,fullfile(out,'paper_value_gate_summary.csv'));
fid=fopen(fullfile(out,'PAPER_VALUE_ASSESSMENT.md'),'w');fprintf(fid,'# Route 32 对论文的价值判定\n\n');
fprintf(fid,'- Energy-matched phase-aperture gate: %d\n',gateA);
fprintf(fid,'- Conditional CRLB high/mid-SNR gate: %d\n',gateB);
fprintf(fid,'- X01 optimized-margin-flip gate: %d\n\n',gateC);
if gateA,fprintf(fid,'相位孔径在固定采样与固定响应灵敏度下仍恢复最弱信息方向，可进入主文机制结果。\n\n');
else,fprintf(fid,'相位孔径结果未通过预定门槛，只能保留为 frozen 反事实或补充分析。\n\n');end
if gateB,fprintf(fid,'正确 EO 下 CRLB 与经验方差同量级，可把条件信息用于模型内实用可恢复性。\n\n');
else,fprintf(fid,'CRLB 未稳定预测经验方差，FIM 应降级为局部解释工具。\n\n');end
if gateC,fprintf(fid,'X01 的完整重优化残差裕量翻转得到保留；偏差投影仅用于解释翻转方向。\n');
else,fprintf(fid,'X01 未复现完整裕量翻转，不应作为核心失效机制。\n');end
fclose(fid);disp(Assessment);
end
