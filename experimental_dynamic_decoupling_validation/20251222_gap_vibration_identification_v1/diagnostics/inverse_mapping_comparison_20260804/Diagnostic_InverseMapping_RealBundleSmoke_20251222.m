%% Diagnostic_InverseMapping_RealBundleSmoke_20251222
% Uses a saved 20251222 diagnostic bundle and the exact extracted F_s
% evaluator. It does not run Main10 and does not overwrite any result.
clear; clc;
thisDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(fileparts(thisDir));
addpath(fullfile(packageRoot,'functions','gap_aware'));
resultFile = fullfile(packageRoot,'diagnostics','projection_trust_region_20260729', ...
    'results','Main_GapAware_VPTopK_FullWave_20251222_B1_S123_diag_projection_smoke.mat');
S = load(resultFile,'Result');
bundle = S.Result.BestWindow.bundle;
vpTable = S.Result.BestWindow.VPSeedTable;
vpSelectedEO = S.Result.BestWindow.VPSelectedEO;
bundle = local_decimate_bundle(bundle, 1200);
nSensor = numel(bundle.sensorIds);
model = repmat(struct('xGrid',[],'evaluate',[]),nSensor,1);
for is = 1:nSensor
    sid = bundle.sensorIds(is);
    it = find([bundle.Template.Sensor.sensor_id] == sid,1);
    ic = find([bundle.CorrectedGapLibrary.sensor.sensorId] == sid,1);
    assert(~isempty(it) && ~isempty(ic),'Missing sensor %d calibration.',sid);
    Tpl = bundle.Template.Sensor(it);
    corr = bundle.CorrectedGapLibrary.sensor(ic);
    rs = bundle.responseSurface;
    model(is).xGrid = Tpl.x_grid(:);
    model(is).evaluate = @(dg,x) local_eval_with_derivative(Tpl,rs,corr,dg,x);
end
cfg = struct('inverseVoltageToleranceMv',1e-7, ...
    'inverseXToleranceMm',1e-9,'inverseMaxBisectionIter',100, ...
    'inverseMinDerivativeMvPerMm',1e-5,'weightFloor',0.05, ...
    'freqSearchHz',[300 1000],'amplitudeLimitMm',0.50, ...
    'dxLimitMm',0.35,'deltaGapLimitMm',0.25, ...
    'deltaGapPriorScaleMm',0.25,'staticRegWeightMv',0.75, ...
    'continuousFrequencyRefine',true,'frequencyRefineHalfWidthHz',2.0);
[inv, summary] = step07jcore.inverse_map_experimental_displacement( ...
    bundle, zeros(nSensor,1), 0, model, cfg);
[seedTable, info] = step07jcore.solve_inverse_mapping_seed(inv,bundle,cfg);
fprintf('Real-bundle smoke: window=%s, valid=%.3f, query-safe=%.3f, ', ...
    string(bundle.windowId),summary.inverse_valid_fraction,summary.query_safe_fraction);
fprintf('inverse voltage RMSE=%.6g mV, candidates=%d, rho=%.6g.\n', ...
    summary.inverse_voltage_rmse,info.candidate_count,info.candidate_rho);
fprintf('Saved VP Top-3 EO: %s\n',mat2str(vpSelectedEO));
fprintf('Inverse Top-3 EO: %s\n',mat2str(seedTable.EO(1:min(3,height(seedTable))).'));
disp(seedTable(1:min(5,height(seedTable)),:));
vpCandidates = vpTable(1:min(3,height(vpTable)),:);
invCandidates = seedTable(1:min(3,height(seedTable)),:);
vpFits = cell(height(vpCandidates),1);
invFits = cell(height(invCandidates),1);
for i=1:height(vpCandidates)
    c=struct('EO',vpCandidates.EO(i),'A',vpCandidates.A(i), ...
        'phi',vpCandidates.phi(i),'dx',vpCandidates.dx(i));
    vpFits{i}=step07jcore.optimize_experimental_full_waveform(bundle,c,model,cfg);
end
for i=1:height(invCandidates)
    c=struct('EO',invCandidates.EO(i),'A',invCandidates.A(i), ...
        'phi',invCandidates.phi(i),'dx',invCandidates.dx(i));
    invFits{i}=step07jcore.optimize_experimental_full_waveform(bundle,c,model,cfg);
end
fprintf('VP full-wave plain RMSE [mV]: %s\n',mat2str(cellfun(@(x)x.plainRmseMv,vpFits),6));
fprintf('Inverse full-wave plain RMSE [mV]: %s\n',mat2str(cellfun(@(x)x.plainRmseMv,invFits),6));
save(fullfile(thisDir,'real_bundle_smoke_20251222.mat'),'inv','summary','seedTable','info','vpTable','vpSelectedEO','vpFits','invFits');

function b = local_decimate_bundle(b, maxPoints)
if numel(b.X) <= maxPoints, return; end
keep = false(numel(b.X),1);
for is = 1:numel(b.sensorIds)
    ids = find(b.sensorIndex == is);
    nk = min(numel(ids), max(1,round(maxPoints*numel(ids)/numel(b.X))));
    pick = unique(round(linspace(1,numel(ids),nk)));
    keep(ids(pick)) = true;
end
fields = {'X','T','TRel','V','W','Theta','S','sensorIndex','F0','Fx'};
for i=1:numel(fields)
    f=fields{i};
    if isfield(b,f) && numel(b.(f)) == numel(keep), b.(f)=b.(f)(keep); end
end
b.pointCount = nnz(keep);
end

function [v, aux] = local_eval_with_derivative(Tpl,rs,corr,dg,x)
[v, aux0] = step07jcore.evaluate_experimental_sensor_voltage(Tpl,rs,corr,dg,x,0,0);
h = 1e-4;
[vp,~] = step07jcore.evaluate_experimental_sensor_voltage(Tpl,rs,corr,dg,x+h,0,0);
[vm,~] = step07jcore.evaluate_experimental_sensor_voltage(Tpl,rs,corr,dg,x-h,0,0);
aux = aux0;
aux.dVoltageDx = (vp-vm)/(2*h);
end
