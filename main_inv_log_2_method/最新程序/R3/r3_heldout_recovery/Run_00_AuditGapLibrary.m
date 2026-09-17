function Result = Run_00_AuditGapLibrary()
%RUN_00_AUDITGAPLIBRARY Audit provenance, continuity and held-out interpolation.

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(fullfile(root,'local_func'),'-begin');P=R3_Protocol();R3=r3_build_context(mainDir,P);
outDir=fullfile(root,'output','00_gap_library_audit');if ~exist(outDir,'dir'),mkdir(outDir);end
x=R3.libCal.xGrid(:);n=numel(P.gAllMm);features=repmat(feature_row(),n,1);
curves=cell(n,1);
for i=1:n
    F=r3_raw_response(R3,P.gAllMm(i));y=F(x);curves{i}=y(:);
    features(i)=struct('gap_mm',P.gAllMm(i),'is_calibration',ismember(P.gAllMm(i),P.gCalMm),...
        'is_held_out_operating',ismember(P.gAllMm(i),P.gOperatingMm),...
        'minimum',min(y),'maximum',max(y),'range',range(y),'mean',mean(y),...
        'area',trapz(x,y),'peak_x_mm',x(find(y==max(y),1)));
end
Features=struct2table(features);

continuity=repmat(continuity_row(),n-1,1);
for i=1:n-1
    d=curves{i+1}-curves{i};scale=max(range(curves{i}),eps);
    continuity(i)=struct('gap_left_mm',P.gAllMm(i),'gap_right_mm',P.gAllMm(i+1),...
        'delta_gap_mm',P.gAllMm(i+1)-P.gAllMm(i),'rmse',rms(d),...
        'normalized_rmse',rms(d)/scale,'max_abs',max(abs(d)),...
        'correlation',corr(curves{i},curves{i+1}));
end
Continuity=struct2table(continuity);

held=repmat(heldout_row(),numel(P.gOperatingMm),1);dx=mean(diff(x));
for i=1:numel(P.gOperatingMm)
    g=P.gOperatingMm(i);raw=r3_raw_response(R3,g);yRaw=raw(x);
    yPred=eval_gap_template(R3.libCal,g,x);d=yPred-yRaw;scale=max(range(yRaw),eps);
    dd=gradient(yPred,dx)-gradient(yRaw,dx);
    held(i)=struct('gap_mm',g,'left_cal_mm',max(P.gCalMm(P.gCalMm<g)),...
        'right_cal_mm',min(P.gCalMm(P.gCalMm>g)),'rmse',rms(d),...
        'normalized_rmse',rms(d)/scale,'max_abs',max(abs(d)),...
        'derivative_rmse_per_mm',rms(dd),'correlation',corr(yRaw,yPred));
end
HeldOut=struct2table(held);
writetable(Features,fullfile(outDir,'gap_response_features.csv'));
writetable(Continuity,fullfile(outDir,'adjacent_response_continuity.csv'));
writetable(HeldOut,fullfile(outDir,'held_out_interpolation_audit.csv'));
Result=struct('P',P,'Features',Features,'Continuity',Continuity,'HeldOut',HeldOut,...
    'maximum_held_out_nrmse',max(HeldOut.normalized_rmse),'outputDir',outDir);
save(fullfile(outDir,'gap_library_audit.mat'),'Result','-v7.3');disp(HeldOut);
end

function r=feature_row(),r=struct('gap_mm',NaN,'is_calibration',false,...
    'is_held_out_operating',false,'minimum',NaN,'maximum',NaN,'range',NaN,...
    'mean',NaN,'area',NaN,'peak_x_mm',NaN);end
function r=continuity_row(),r=struct('gap_left_mm',NaN,'gap_right_mm',NaN,...
    'delta_gap_mm',NaN,'rmse',NaN,'normalized_rmse',NaN,'max_abs',NaN,...
    'correlation',NaN);end
function r=heldout_row(),r=struct('gap_mm',NaN,'left_cal_mm',NaN,...
    'right_cal_mm',NaN,'rmse',NaN,'normalized_rmse',NaN,'max_abs',NaN,...
    'derivative_rmse_per_mm',NaN,'correlation',NaN);end
