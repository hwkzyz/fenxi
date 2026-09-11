%% Production low-speed transfer: frozen registration, compare M0 and M1.
% M0: q=0, nonlinear g0 only. M1: nonlinear [g0,q]. The no-intercept
% cross-platform gain is solved analytically for every candidate.
clear; clc;

thisDir = fileparts(mfilename('fullpath'));
caseDir = fileparts(fileparts(thisDir));
outDir = fullfile(thisDir, 'outputs');
surfaceFile = fullfile(outDir, 'ProductionResponseSurface_20250527.mat');
sourceLibraryFile = fullfile(caseDir, 'results', ...
    'Step06I_OffsetTiltShared_GapLibrary_20250527_B1_S136_nob_harddomain_fixed12pctroi_formal_20260830.mat');
registrationFile = fullfile(caseDir, 'analysis', 'transfer_identifiability_20260830', ...
    'outputs', 'NestedPlatformModels_20250527.csv');

S = load(surfaceFile, 'ProductionResponseSurface', 'GateSummary');
L = load(sourceLibraryFile, 'CorrectedGapLibrary');
Registration = readtable(registrationFile, 'TextType', 'string');
F = S.ProductionResponseSurface;
CorrectedGapLibrary = L.CorrectedGapLibrary;

gDomain = F.directGapDomainMm;
finiteX = all(isfinite(F.coeff), 2);
xDomain = [min(F.xGrid(finiteX)), max(F.xGrid(finiteX))];
qBounds = [-0.08, 0.08];
gGrid = linspace(gDomain(1), gDomain(2), 281);
qGrid = linspace(qBounds(1), qBounds(2), 161);

rows = {};
profileRows = {};
for is = 1:numel(CorrectedGapLibrary.sensor)
    C = CorrectedGapLibrary.sensor(is);
    sid = C.sensorId;
    rowReg = Registration(Registration.sensorId == sid & Registration.model == "P2_registration", :);
    assert(height(rowReg) == 1, 'Missing frozen P2 registration for CH%d.', sid);
    tau = rowReg.tauMm;
    k = rowReg.xScale;
    x = C.x(:);
    v = C.vLowMv(:);
    xLib = k .* (x - tau);
    use = isfinite(x) & isfinite(v) & xLib >= xDomain(1) & xLib <= xDomain(2);
    xUse = x(use); vUse = v(use); xLibUse = xLib(use);
    xCenter = 0.5 * (min(xLibUse) + max(xLibUse));

    [m0, p0] = fit_model_local(0, gGrid, qGrid, false, xLibUse, vUse, xCenter, F, gDomain);
    [m1, p1] = fit_model_local([m0.g0Mm, 0], gGrid, qGrid, true, xLibUse, vUse, xCenter, F, gDomain);
    improvementPct = 100 * (m0.rmseMv - m1.rmseMv) / max(m0.rmseMv, eps);
    m1Boundary = m1.gMarginMm < 1e-4 || m1.qBoundary;
    useM1 = improvementPct >= 5 && ~m1Boundary;
    if useM1, selected = m1; selectedName = "M1_g0_q"; else, selected = m0; selectedName = "M0_g0"; end

    CorrectedGapLibrary.sensor(is).g0Mm = selected.g0Mm;
    CorrectedGapLibrary.sensor(is).tauMm = tau;
    CorrectedGapLibrary.sensor(is).xScale = k;
    CorrectedGapLibrary.sensor(is).muGapPerXMm = selected.q;
    CorrectedGapLibrary.sensor(is).tiltAngleDeg = atan(selected.q) * 180 / pi;
    CorrectedGapLibrary.sensor(is).voltageGain = selected.gain;
    CorrectedGapLibrary.sensor(is).voltageOffsetMv = 0;
    CorrectedGapLibrary.sensor(is).lowFitRmseMv = selected.rmseMv;
    CorrectedGapLibrary.sensor(is).domainFeasible = true;
    CorrectedGapLibrary.sensor(is).gDomainMarginMm = selected.gMarginMm;
    CorrectedGapLibrary.sensor(is).xDomainMarginMm = min(xLibUse-xDomain(1), [], 'omitnan');
    CorrectedGapLibrary.sensor(is).activeBoundary = selected.gMarginMm < 1e-4;
    CorrectedGapLibrary.sensor(is).x = xUse;
    CorrectedGapLibrary.sensor(is).vLowMv = vUse;
    CorrectedGapLibrary.sensor(is).vFitMv = selected.pred;
    CorrectedGapLibrary.sensor(is).xLib = xLibUse;
    CorrectedGapLibrary.sensor(is).gEff = selected.gPath;
    CorrectedGapLibrary.sensor(is).residualMv = vUse-selected.pred;
    CorrectedGapLibrary.sensor(is).transferModel = char(selectedName);

    rows{end+1,1} = table(sid, tau, k, m0.g0Mm, m0.gain, m0.rmseMv, ...
        m1.g0Mm, m1.q, m1.gain, m1.rmseMv, improvementPct, m1.gMarginMm, ...
        m1.qBoundary, selectedName, ...
        'VariableNames', {'sensorId','frozenTauMm','frozenXScale','M0_g0Mm', ...
        'M0_gain','M0_rmseMv','M1_g0Mm','M1_q','M1_gain','M1_rmseMv', ...
        'M1_improvementPct','M1_gMarginMm','M1_qBoundary','selectedModel'}); %#ok<AGROW>
    profileRows{end+1,1} = [add_profile_id_local(p0, sid, "M0_g0"); ...
        add_profile_id_local(p1, sid, "M1_g0_q")]; %#ok<AGROW>
    fprintf('CH%d: M0 %.2f mV; M1 %.2f mV (%+.1f%%), g0=%.4f q=%+.5f, selected=%s\n', ...
        sid, m0.rmseMv, m1.rmseMv, improvementPct, m1.g0Mm, m1.q, selectedName);
end

LowSpeedTransfer = vertcat(rows{:});
LowSpeedProfile = vertcat(profileRows{:});
CorrectedGapLibrary.responseSurface = F;
CorrectedGapLibrary.responseFile = surfaceFile;
CorrectedGapLibrary.method = 'production_one_way_joint_surface_no_intercept';
CorrectedGapLibrary.description = ['Production low-speed transfer using the accepted one-way joint static surface. ' ...
    'Registration is frozen from P2; M0/M1 estimate g0/q with analytic no-intercept gain.'];
CorrectedGapLibrary.formula = ['T_high_nv = T_low + a_s * (F(g_ref + delta_g, x_lib) - ' ...
    'F(g_ref, x_lib)); b_s=0; vibration enters x.'];
CorrectedGapLibrary.productionCalibrationGate = S.GateSummary;
CorrectedGapLibrary.transferModelSelection = LowSpeedTransfer;

outFile = fullfile(outDir, 'ProductionGapLibrary_20250527_B1_S136.mat');
save(outFile, 'CorrectedGapLibrary', 'LowSpeedTransfer', 'LowSpeedProfile', '-v7.3');
writetable(LowSpeedTransfer, fullfile(outDir, 'ProductionLowSpeedTransfer_20250527.csv'));
writetable(LowSpeedProfile, fullfile(outDir, 'ProductionLowSpeedProfiles_20250527.csv'));
fprintf('Saved production low-speed library: %s\n', outFile);

function [best, profile] = fit_model_local(start, gGrid, qGrid, useQ, xLib, v, xCenter, F, domain)
if useQ
    J = nan(numel(gGrid), numel(qGrid));
    for ig = 1:numel(gGrid)
        for iq = 1:numel(qGrid)
            J(ig,iq) = objective_local([gGrid(ig),qGrid(iq)], xLib, v, xCenter, F, domain);
        end
    end
    [~,ix] = min(J(:)); [ig,iq] = ind2sub(size(J),ix); p0=[gGrid(ig),qGrid(iq)];
    if numel(start)==2 && isfinite(objective_local(start,xLib,v,xCenter,F,domain)) && ...
            objective_local(start,xLib,v,xCenter,F,domain)<objective_local(p0,xLib,v,xCenter,F,domain), p0=start; end
    lb=[domain(1),min(qGrid)]; ub=[domain(2),max(qGrid)];
else
    J = nan(numel(gGrid),1);
    for ig=1:numel(gGrid), J(ig)=objective_local([gGrid(ig),0],xLib,v,xCenter,F,domain); end
    [~,ig]=min(J); p0=[gGrid(ig),0]; lb=[domain(1),0]; ub=[domain(2),0];
end
if exist('fmincon','file')==2
    opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',500, ...
        'MaxFunctionEvaluations',3000,'OptimalityTolerance',1e-9,'StepTolerance',1e-10);
    p=fmincon(@(z)objective_local(z,xLib,v,xCenter,F,domain),p0,[],[],[],[],lb,ub,[],opts);
else
    p=p0;
end
[mse,detail]=objective_local(p,xLib,v,xCenter,F,domain);
best=detail; best.rmseMv=sqrt(mse);
best.qBoundary=useQ && (abs(p(2)-lb(2))<1e-4 || abs(p(2)-ub(2))<1e-4);
if useQ
    [GG,QQ]=ndgrid(gGrid,qGrid); profile=table(GG(:),QQ(:),sqrt(J(:)), ...
        'VariableNames',{'g0Mm','q','rmseMv'});
else
    profile=table(gGrid(:),zeros(numel(gGrid),1),sqrt(J(:)), ...
        'VariableNames',{'g0Mm','q','rmseMv'});
end
end

function [mse,detail]=objective_local(p,xLib,v,xCenter,F,domain)
gPath=p(1)+p(2).*(xLib-xCenter);
if any(gPath<domain(1) | gPath>domain(2)), mse=1e12; detail=empty_detail_local(); return; end
raw=eval_surface_local(F,gPath,xLib);
ok=isfinite(raw)&isfinite(v);
if nnz(ok)<10, mse=1e12; detail=empty_detail_local(); return; end
gain=(raw(ok)'*v(ok))/max(raw(ok)'*raw(ok),eps);
gain=max(gain,0);
pred=gain.*raw;
mse=mean((v(ok)-pred(ok)).^2);
detail=struct('g0Mm',p(1),'q',p(2),'gain',gain,'pred',pred,'gPath',gPath, ...
    'gMarginMm',min([min(gPath)-domain(1),domain(2)-max(gPath)]));
end

function Fv=eval_surface_local(F,g,x)
B=nan(numel(x),3);
for k=1:3, B(:,k)=interp1(F.xGrid,F.coeff(:,k),x,'linear',NaN); end
Fv=B(:,1)+B(:,2)./g+B(:,3).*log(g./F.g0Mm);
end

function detail=empty_detail_local()
detail=struct('g0Mm',NaN,'q',NaN,'gain',NaN,'pred',NaN,'gPath',NaN,'gMarginMm',-Inf);
end

function T=add_profile_id_local(T,sid,name)
T.sensorId=repmat(sid,height(T),1); T.model=repmat(name,height(T),1); T=movevars(T,{'sensorId','model'},'Before',1);
end
