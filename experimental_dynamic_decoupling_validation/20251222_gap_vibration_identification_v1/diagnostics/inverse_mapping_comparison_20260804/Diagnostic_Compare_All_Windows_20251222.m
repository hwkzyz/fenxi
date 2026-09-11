%% Diagnostic_Compare_All_Windows_20251222
% Isolated 18-window VP versus inverse-mapping comparison. Formal Main10 is
% not called or modified; saved bundles come from the existing diagnostic.
clear; clc;
thisDir=fileparts(mfilename('fullpath')); packageRoot=fileparts(fileparts(thisDir));
addpath(fullfile(packageRoot,'functions','gap_aware'));
resultFile=fullfile(packageRoot,'diagnostics','projection_trust_region_20260729','results', ...
    'Main_GapAware_VPTopK_FullWave_20251222_B1_S123_diag_projection_allW_fullpoints.mat');
S=load(resultFile,'Result'); W=S.Result.WindowResult; nW=numel(W);
cfg=struct('inverseVoltageToleranceMv',1e-7,'inverseXToleranceMm',1e-9, ...
    'inverseMaxBisectionIter',100,'inverseMinDerivativeMvPerMm',1e-5, ...
    'weightFloor',0.05,'freqSearchHz',[300 1000],'amplitudeLimitMm',0.50, ...
    'dxLimitMm',0.35,'deltaGapLimitMm',0.25,'deltaGapPriorScaleMm',0.25, ...
    'staticRegWeightMv',0.75,'continuousFrequencyRefine',true, ...
    'frequencyRefineHalfWidthHz',2.0);
rows=repmat(struct('windowId',NaN,'inverseValidFraction',NaN,'inverseVoltageRmseMv',NaN, ...
    'querySafeFraction',NaN,'branchSwitchRate',NaN,'candidateRho',NaN, ...
    'vpTopEO',"",'inverseTopEO',"",'vpBestPlainRmseMv',NaN, ...
    'inverseBestPlainRmseMv',NaN,'vpBestEO',NaN,'inverseBestEO',NaN),nW,1);
for iw=1:nW
    bundle=local_decimate_bundle(W(iw).bundle,1200); model=local_make_model(bundle); nSensor=numel(bundle.sensorIds);
    [inv,sm]=step07jcore.inverse_map_experimental_displacement(bundle,zeros(nSensor,1),0,model,cfg);
    [invTable,info]=step07jcore.solve_inverse_mapping_seed(inv,bundle,cfg); vpTable=W(iw).VPSeedTable;
    vpTop=vpTable(1:min(3,height(vpTable)),:); invTop=invTable(1:min(3,height(invTable)),:);
    vpFits=cell(height(vpTop),1); invFits=cell(height(invTop),1);
    for j=1:height(vpTop), c=struct('EO',vpTop.EO(j),'A',vpTop.A(j),'phi',vpTop.phi(j),'dx',vpTop.dx(j)); vpFits{j}=step07jcore.optimize_experimental_full_waveform(bundle,c,model,cfg); end
    for j=1:height(invTop), c=struct('EO',invTop.EO(j),'A',invTop.A(j),'phi',invTop.phi(j),'dx',invTop.dx(j)); invFits{j}=step07jcore.optimize_experimental_full_waveform(bundle,c,model,cfg); end
    [vpRmse,jv]=min(cellfun(@(x)x.plainRmseMv,vpFits)); [inRmse,ji]=min(cellfun(@(x)x.plainRmseMv,invFits));
    rows(iw).windowId=W(iw).windowId; rows(iw).inverseValidFraction=sm.inverse_valid_fraction; rows(iw).inverseVoltageRmseMv=sm.inverse_voltage_rmse;
    rows(iw).querySafeFraction=sm.query_safe_fraction; rows(iw).branchSwitchRate=sm.branch_switch_rate; rows(iw).candidateRho=info.candidate_rho;
    rows(iw).vpTopEO=string(mat2str(vpTop.EO.')); rows(iw).inverseTopEO=string(mat2str(invTop.EO.')); rows(iw).vpBestPlainRmseMv=vpRmse; rows(iw).inverseBestPlainRmseMv=inRmse;
    rows(iw).vpBestEO=vpFits{jv}.EO; rows(iw).inverseBestEO=invFits{ji}.EO;
    fprintf('W%02d valid %.3f rho %.4f VP %d/%.3f INV %d/%.3f\n',W(iw).windowId,sm.inverse_valid_fraction,info.candidate_rho,rows(iw).vpBestEO,vpRmse,rows(iw).inverseBestEO,inRmse);
end
summary=struct2table(rows); outFile=fullfile(thisDir,'InverseMapping_vs_VP_All18_20251222.mat'); save(outFile,'summary','rows','resultFile'); writetable(summary,fullfile(thisDir,'InverseMapping_vs_VP_All18_20251222.csv')); fprintf('Saved %s\n',outFile);

function model=local_make_model(bundle)
nSensor=numel(bundle.sensorIds); model=repmat(struct('xGrid',[],'evaluate',[]),nSensor,1);
for is=1:nSensor, sid=bundle.sensorIds(is); it=find([bundle.Template.Sensor.sensor_id]==sid,1); ic=find([bundle.CorrectedGapLibrary.sensor.sensorId]==sid,1); Tpl=bundle.Template.Sensor(it); corr=bundle.CorrectedGapLibrary.sensor(ic); rs=bundle.responseSurface; model(is).xGrid=Tpl.x_grid(:); model(is).evaluate=@(dg,x)local_eval(Tpl,rs,corr,dg,x); end
end
function [v,aux]=local_eval(Tpl,rs,corr,dg,x)
[v,aux0]=step07jcore.evaluate_experimental_sensor_voltage(Tpl,rs,corr,dg,x,0,0); h=1e-4; [vp,~]=step07jcore.evaluate_experimental_sensor_voltage(Tpl,rs,corr,dg,x+h,0,0); [vm,~]=step07jcore.evaluate_experimental_sensor_voltage(Tpl,rs,corr,dg,x-h,0,0); aux=aux0; aux.dVoltageDx=(vp-vm)/(2*h);
end
function b=local_decimate_bundle(b,maxPoints)
if numel(b.X)<=maxPoints,return;end; keep=false(numel(b.X),1);
for is=1:numel(b.sensorIds), ids=find(b.sensorIndex==is); nk=min(numel(ids),max(1,round(maxPoints*numel(ids)/numel(b.X)))); keep(ids(unique(round(linspace(1,numel(ids),nk)))))=true; end
ff={'X','T','TRel','V','W','Theta','S','sensorIndex','F0','Fx'};
for i=1:numel(ff), f=ff{i}; if isfield(b,f)&&numel(b.(f))==numel(keep), b.(f)=b.(f)(keep); end, end; b.pointCount=nnz(keep);
end
