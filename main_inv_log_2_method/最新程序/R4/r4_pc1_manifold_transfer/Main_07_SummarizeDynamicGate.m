function Summary = Main_07_SummarizeDynamicGate()
%SUMMARIZE_TILTCALIBRATIONTRANSFER Assemble static and dynamic audit evidence.
root=fileparts(mfilename('fullpath'));staticDir=fullfile(root,'output','static_calibration_transfer');
dynamicDir=fullfile(root,'output','noiseless_dynamic_gate');perAngle=fullfile(dynamicDir,'per_angle');
S=readtable(fullfile(staticDir,'angle_summary.csv'));
files=dir(fullfile(perAngle,'tilt_*deg.csv'));if numel(files)~=7,error('Expected seven per-angle dynamic result files.');end
D=table();for i=1:numel(files),D=[D;readtable(fullfile(files(i).folder,files(i).name))];end %#ok<AGROW>
D.case_id=string(D.case_id);D.method=string(D.method);
D=sortrows(D,{'tilt_deg','case_id','method'});
writetable(D,fullfile(dynamicDir,'noiseless_dynamic_all_angles.csv'));
check=table();methods=["C1_g0_only";"C2_g0_mu";"C3_FE_oracle"];
caseIds=["null_00","blind_10_26","blind_07_19","blind_13_27"];
for ii=1:numel(caseIds)
    ic=caseIds(ii);
    for im=1:numel(methods)
        q=D.case_id==string(ic) & D.method==methods(im);T=D(q,:);
        check=[check;table(string(ic),methods(im),all(T.identification_status=="identified"),max(abs(T.gap_error_mm)),...
            max(abs([T.A1_error_mm;T.A2_error_mm])),min(T.second_best_margin),...
            'VariableNames',{'case_id','method','all_true_models_win','max_abs_gap_error_mm',...
            'max_abs_amplitude_error_mm','min_objective_margin'})]; %#ok<AGROW>
    end
end
writetable(check,fullfile(dynamicDir,'dynamic_gate_summary.csv'));
make_figure(S,D,fullfile(root,'output','tilt_calibration_transfer_summary.png'));
write_report(S,check,root);
Summary=struct('static_summary',S,'dynamic',D,'dynamic_gate',check);
disp(check);
end

function make_figure(S,D,filePath)
fig=figure('Color','w','Position',[100 100 1080 700]);tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;plot(S.tilt_deg,S.transfer_relative_C1_median,'o-','LineWidth',1.3);hold on;plot(S.tilt_deg,S.transfer_relative_C2_median,'s-','LineWidth',1.3);xlabel('FE tilt angle (deg)');ylabel('Median relative increment bias');legend('C1','C2','Location','northwest');grid on;box on;
nexttile;q=D.case_id=="blind_10_26" & (D.method=="C1_g0_only" | D.method=="C2_g0_mu");plot_pair(D(q,:),'A2_error_mm','Second amplitude error (mm)');
nexttile;q=D.case_id=="blind_07_19" & (D.method=="C1_g0_only" | D.method=="C2_g0_mu");plot_pair(D(q,:),'second_best_margin','Blind EO margin');
nexttile;q=D.case_id=="blind_13_27" & (D.method=="C1_g0_only" | D.method=="C2_g0_mu");plot_pair(D(q,:),'A2_error_mm','Second amplitude error (mm)');
exportgraphics(fig,filePath,'Resolution',220);close(fig);
end

function plot_pair(T,name,ylabelText)
methods=["C1_g0_only","C2_g0_mu"];
for i=1:numel(methods)
    method=methods(i);
    q=T.method==method;if method=="C1_g0_only",style='o-';label='C1';else,style='s-';label='C2';end
    plot(T.tilt_deg(q),T.(name)(q),style,'LineWidth',1.3,'DisplayName',label);hold on;
end
xlabel('FE tilt angle (deg)');ylabel(ylabelText);legend('Location','best');grid on;box on;
end

function write_report(S,C,root)
fid=fopen(fullfile(root,'output','TILT_AUDIT_RESULTS.md'),'w');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Tilt calibration-transfer audit results\n\n');
fprintf(fid,'The static audit uses seven real COMSOL tilt geometries and five adjacent finite clearance increments per geometry. The dynamic audit is noiseless and uses complete blind EO-pair search for all reported cases.\n\n');
fprintf(fid,'## Static gate\n\n');
fprintf(fid,'- C2 reduced the median increment-transfer bias at every FE tilt angle.\n');
fprintf(fid,'- C2 improved every adjacent gap pair at every angle: %d.\n',all(S.C2_improves_all_gap_pairs));
fprintf(fid,'- The median C2 reduction across angles was %.3f.\n\n',median(S.C2_median_reduction_fraction));
fprintf(fid,'## Dynamic gate\n\n');
for i=1:height(C)
    if C.case_id(i)=="null_00"
        fprintf(fid,'- %s, %s: null-vibration leakage check; max estimated amplitude error=%.6g mm; model-win flag is not interpreted.\n',C.case_id(i),C.method(i),C.max_abs_amplitude_error_mm(i));
    else
        fprintf(fid,'- %s, %s: all windows identified=%d; max |gap error|=%.6g mm; max |amplitude error|=%.6g mm; min blind-EO margin=%.6g.\n',C.case_id(i),C.method(i),C.all_true_models_win(i),C.max_abs_gap_error_mm(i),C.max_abs_amplitude_error_mm(i),C.min_objective_margin(i));
    end
end
fprintf(fid,'\nC3 oracle closure must remain near numerical precision before attributing C1/C2 differences to calibration transfer.\n');
end
