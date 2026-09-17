function Result = Run_EO29_ClearanceConfoundingAttribution()
%RUN_EO29_CLEARANCECONFOUNDINGATTRIBUTION Frozen standalone EO29 audit.
% This script does not change the production solver or theory files.

thisDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(thisDir);
addpath(mainDir, '-begin');
addpath(fullfile(mainDir, 'local_func'), '-begin');
outDir = fullfile(thisDir, 'output');
if ~exist(outDir, 'dir'), mkdir(outDir); end

ctx = load_inv_log_2_project_context(mainDir);
trust = load_fixed_trust_domain(mainDir);
cfg = ctx.cfgAna;
cfg.NumRevs_high = 8;
cfg.route30UseAllTrustedSamples = true;

% The audit intentionally uses the absolute static library.  The reference
% response is fixed at g_low: no fit in this script can re-free clearance.
lib = build_gap_template_library(ctx.gapList, ctx.xCell, ctx.yCell, NaN, cfg.xGridN);
lib = restrict_gap_template_domain(lib, trust.domain);
gRef = 0.5;
gHighList = gRef + [0.005, 0.050, 0.300];
eoList = [23, 29, 30];
map = make_common_support_map(cfg, lib.domain);

sensorRows = repmat(sensor_row0(), 0, 1);
legacyRows = repmat(score_row0(), 0, 1);
operatorRows = repmat(operator_row0(), 0, 1);
frozenAuditRows = repmat(frozen_audit_row0(), 0, 1);
fitRows = repmat(fit_row0(), 0, 1);
for ig = 1:numel(gHighList)
    gHigh = gHighList(ig);
    [e, q, fRef, fHigh] = clearance_residual(lib, gRef, gHigh, map.x_v);
    sensorRows = [sensorRows; per_sensor_rows(map, e, q, gHigh-gRef)]; %#ok<AGROW>

    % Retained only to document why the old apparent-shift screen is invalid
    % for causal attribution: the tangent projection explains negligible energy.
    uApp = apparent_shift_by_sensor(map, e, q);
    for ie = 1:numel(eoList)
        eo = eoList(ie);
        legacyRows(end+1,1) = make_score_row(map, e, q, uApp, eo, gHigh-gRef, "geometry_only_deprecated"); %#ok<AGROW>
        legacyRows(end+1,1) = make_score_row(map, e, q, uApp, eo, gHigh-gRef, "frozen_passage_legacy"); %#ok<AGROW>
        legacyRows(end+1,1) = make_score_row(map, e, q, uApp, eo, gHigh-gRef, "full_waveform_legacy"); %#ok<AGROW>
        for ablation = ["full","no_intra_passage","no_angular_diversity","flattened_sensitivity"]
            operatorRows(end+1,1) = operator_gamma(map, e, q, eo, gHigh-gRef, ablation); %#ok<AGROW>
        end
        frozenAuditRows(end+1,1) = frozen_operator_audit(map,e,q,eo,gHigh-gRef); %#ok<AGROW>
        fitRows(end+1,1) = fixed_reference_fit(map, lib, fHigh, gRef, eo, gHigh-gRef); %#ok<AGROW>
    end
end

SensorTable = struct2table(sensorRows);
ScoreTable = struct2table(legacyRows);
OperatorTable = struct2table(operatorRows);
FrozenAuditTable = struct2table(frozenAuditRows);
FitTable = struct2table(fitRows);
write_tables_and_report(outDir, SensorTable, ScoreTable, OperatorTable, FrozenAuditTable, FitTable, cfg, trust, gRef, eoList);
save(fullfile(outDir, 'eo29_clearance_confounding_attribution.mat'), ...
    'SensorTable', 'ScoreTable', 'OperatorTable', 'FrozenAuditTable', 'FitTable', 'cfg', 'trust', 'gRef', 'eoList', 'map', '-v7.3');
Result = struct('sensorTable', SensorTable, 'scoreTable', ScoreTable, ...
    'operatorTable', OperatorTable, 'frozenAuditTable', FrozenAuditTable, 'fitTable', FitTable, 'outputDir', outDir);
end

function map = make_common_support_map(cfg, domain)
% Delegate support and OPR-angle construction to the maintained mapper.
omega = cfg.RPM_high * 2*pi/60;
T = 2*pi/omega;
dt = 1/cfg.fs;
t = (0:dt:(cfg.NumRevs_high+1)*T).';
Tobr = 0.08*T + (0:cfg.NumRevs_high).'*T;
Data = struct('t',t,'T_opr_truth',Tobr,'V_cap',zeros(size(t)));
map = map_highspeed_to_space(Data,cfg.alpha_k,cfg.R_tip,domain,0);
map.theta_center_v = 2*pi*(map.rev_v(:)-1) + cfg.alpha_k(map.S_v(:)).';
map.R_tip = cfg.R_tip;
map.n_passages = cfg.NumRevs_high*numel(cfg.alpha_k);
end

function [e,q,fRef,fHigh] = clearance_residual(lib, gRef, gHigh, x)
fRef = eval_gap_template(lib, gRef, x);
fHigh = eval_gap_template(lib, gHigh, x);
e = fHigh - fRef;
q = -eval_gap_derivative(lib, gRef, x); % dF(x-u)/du at u=0
end

function rows = per_sensor_rows(map,e,q,dg)
ids=unique(map.S_v(:),'stable'); rows=repmat(sensor_row0(),numel(ids),1);
for i=1:numel(ids)
    ix=map.S_v==ids(i); ee=e(ix); qq=q(ix);
    u=(qq'*ee)/(qq'*qq);
    proj=qq*u; residual=ee-proj;
    rows(i)=struct('delta_g_um',1000*dg,'sensor_id',ids(i),'n_samples',nnz(ix), ...
        'clearance_residual_rms_V',rms(ee),'apparent_shift_mm',u, ...
        'translation_explained_fraction',sum(proj.^2)/sum(ee.^2), ...
        'translation_residual_rms_V',rms(residual));
end
end

function u = apparent_shift_by_sensor(map,e,q)
ids=unique(map.S_v(:),'stable'); u=zeros(numel(ids),1);
for i=1:numel(ids)
    ix=map.S_v==ids(i);u(i)=(q(ix)'*e(ix))/(q(ix)'*q(ix));
end
end

function r = make_score_row(map,e,q,uApp,eo,dg,kind)
switch kind
    case "geometry_only_deprecated"
        alpha = unique(map.theta_center_v(:),'stable');
        % unique centre angle includes turns; retain one centre for each sensor.
        alpha = mod(alpha(1:numel(uApp)), 2*pi);
        X = [sin(eo*alpha), cos(eo*alpha)];
        fitted = X*(X\uApp); target=uApp;
        score = sum(fitted.^2)/sum(target.^2);
        rss0=sum(target.^2); rss=sum((target-fitted).^2);
    case "frozen_passage_legacy"
        X=[q, q.*sin(eo*map.theta_center_v), q.*cos(eo*map.theta_center_v)];
        [score,rss0,rss]=conditioned_score(e,X);
    case "full_waveform_legacy"
        X=[q, q.*sin(eo*map.theta_v), q.*cos(eo*map.theta_v)];
        [score,rss0,rss]=conditioned_score(e,X);
    otherwise
        error('Unknown score kind.');
end
r=struct('delta_g_um',1000*dg,'eo',eo,'layer',string(kind), ...
    'explained_fraction',score,'rss_baseline_V2',rss0,'rss_after_eo_V2',rss, ...
    'support_samples',numel(e));
end

function [score,rss0,rss] = conditioned_score(y,X)
X0=X(:,1); Xeo=X(:,2:3); Q0=orth(X0); y0=y-Q0*(Q0'*y); Xe=Xeo-Q0*(Q0'*Xeo);
Qe=orth(Xe); residual=y0-Qe*(Qe'*y0); rss0=sum(y0.^2); rss=sum(residual.^2);
score=max(0,1-rss/max(rss0,eps));
end

function r = operator_gamma(map,e,q,eo,dg,ablation)
% Energy-consistent attribution: all four arms target the same e and
% residualize the same global registration direction q before EO projection.
theta = map.theta_v(:); thetaC = map.theta_center_v(:); x = map.x_v(:);
switch ablation
    case "full"
        thetaUse = theta;
        qUse = q;
        definition = "B=[q sin(k theta), q cos(k theta)]";
    case "no_intra_passage"
        thetaUse = thetaC;
        qUse = q;
        definition = "theta=theta_center; q retained";
    case "no_angular_diversity"
        alphaRef = mod(thetaC(find(map.S_v==1,1,'first')), 2*pi);
        thetaUse = 2*pi*(map.rev_v(:)-1) + alphaRef + x/map.R_tip;
        qUse = q;
        definition = "all sensors use alpha_sensor1; x/R retained";
    case "flattened_sensitivity"
        thetaUse = theta;
        % A constant envelope with the same global RMS as q preserves the
        % total EO-column scale while removing spatial sensitivity shaping.
        qUse = rms(q) * ones(size(q));
        definition = "q_flat=rms(q) over formal support";
    otherwise
        error('Unknown operator ablation.');
end
% q is always the registration nuisance direction; only the EO envelope is
% flattened in arm D.
X0 = q; B = [qUse.*sin(eo*thetaUse), qUse.*cos(eo*thetaUse)];
Q0 = orth(X0); e0 = e-Q0*(Q0'*e); B0 = B-Q0*(Q0'*B);
[U,S,~] = svd(B0,'econ'); sv = diag(S);
tol = max(size(B0))*eps(max([sv; 1])); rankB = sum(sv>tol);
if rankB==0
    gamma=0; residual=e0; condB=Inf;
else
    QB=U(:,1:rankB); projection=QB*(QB'*e0); residual=e0-projection;
    gamma=sum(projection.^2)/max(sum(e0.^2),eps);
    condB=sv(1)/sv(rankB);
end
r=struct('delta_g_um',1000*dg,'eo',eo,'ablation',string(ablation), ...
    'operator_definition',string(definition),'gamma',gamma, ...
    'residualized_energy_V2',sum(e0.^2),'unexplained_energy_V2',sum(residual.^2), ...
    'B_rank',rankB,'B_condition_number',condB,'B_singular_value_1',sv(1), ...
    'B_singular_value_2',get_sv(sv,2),'support_samples',numel(e));
end

function v=get_sv(x,k),if numel(x)>=k,v=x(k);else,v=0;end,end

function r = frozen_operator_audit(map,e,q,eo,dg)
% With theta frozen, B is in the span of the three per-sensor q blocks.
% This direct audit distinguishes that fact from repetition across revolutions.
sid=map.S_v(:); ids=unique(sid,'stable'); thetaC=map.theta_center_v(:);
Qsensor=zeros(numel(q),numel(ids)); inner=zeros(numel(ids),1);
for j=1:numel(ids)
    ix=sid==ids(j); Qsensor(ix,j)=q(ix);
end
Qglobal=orth(q); e0=e-Qglobal*(Qglobal'*e);
for j=1:numel(ids),ix=sid==ids(j);inner(j)=q(ix)'*e0(ix);end
B=[q.*sin(eo*thetaC),q.*cos(eo*thetaC)];
QB=orth(Qsensor); blockResidual=B-QB*(QB'*B);
Qp=orth(Qsensor); ep=e-Qp*(Qp'*e);
r=struct('delta_g_um',1000*dg,'eo',eo,'global_q_residual_energy_V2',sum(e0.^2), ...
    'sensor1_q_dot_global_e0',inner(1),'sensor2_q_dot_global_e0',inner(2), ...
    'sensor3_q_dot_global_e0',inner(3),'max_abs_sensor_q_dot_global_e0',max(abs(inner)), ...
    'frozen_B_outside_per_sensor_q_block_fraction',sum(blockResidual.^2,'all')/max(sum(B.^2,'all'),eps), ...
    'per_sensor_q_block_rank',rank(Qsensor),'per_sensor_residual_energy_V2',sum(ep.^2), ...
    'support_samples',numel(e));
end

function r = fixed_reference_fit(map,lib,y,gRef,eo,dg)
% Same support, dx/amplitude bounds, starts, and iteration budget for all EO.
theta=map.theta_v(:);x=map.x_v(:);starts=[];
for dx=[-.25 0 .25]
    for a=[.05 .25 .50]
        for ph=[0 pi/2 pi 3*pi/2]
            starts(end+1,:)=[dx,a*cos(ph),a*sin(ph)]; %#ok<AGROW>
        end
    end
end
lb=[-.5 -1.5 -1.5];ub=[.5 1.5 1.5];best=inf;zBest=[];exitBest=NaN;iterBest=NaN;
for k=1:size(starts,1)
    z0=starts(k,:); obj=@(z) bounded_sse(z,lb,ub,theta,x,y,lib,gRef,eo);
    opt=optimset('Display','off','MaxIter',800,'MaxFunEvals',5000,'TolX',1e-10,'TolFun',1e-12);
    [z,val,flag,out]=fminsearch(obj,z0,opt);
    if val<best,best=val;zBest=z;exitBest=flag;iterBest=out.iterations;end
end
u=zBest(2)*sin(eo*theta)+zBest(3)*cos(eo*theta);
fit=eval_gap_template(lib,gRef,x-zBest(1)-u);res=y-fit;
r=struct('delta_g_um',1000*dg,'eo',eo,'g_reference_mm',gRef, ...
    'g_free',false,'dx_bound_mm',.5,'amplitude_bound_mm',1.5, ...
    'start_count',size(starts,1),'max_iterations',800,'dx_fit_mm',zBest(1), ...
    'amplitude_fit_mm',hypot(zBest(2),zBest(3)),'rmse_V',rms(res), ...
    'sse_V2',sum(res.^2),'exitflag',exitBest,'iterations',iterBest, ...
    'support_samples',numel(y));
end

function v = bounded_sse(z,lb,ub,theta,x,y,lib,gRef,eo)
if any(z<lb) || any(z>ub), v=1e12+1e10*sum(max(lb-z,0).^2+max(z-ub,0).^2); return; end
u=z(2)*sin(eo*theta)+z(3)*cos(eo*theta); f=eval_gap_template(lib,gRef,x-z(1)-u);
if any(~isfinite(f)),v=1e12;else,v=sum((y-f).^2);end
end

function write_tables_and_report(outDir,S,G,O,A,F,cfg,trust,gRef,eoList)
writetable(S,fullfile(outDir,'per_sensor_clearance_residual.csv'));
% Remove the superseded filename so it cannot be mistaken for attribution.
oldScore=fullfile(outDir,'eo_attribution_scores.csv');
if exist(oldScore,'file'), delete(oldScore); end
writetable(G,fullfile(outDir,'deprecated_apparent_shift_scores.csv'));
writetable(O,fullfile(outDir,'operator_ablation_gamma.csv'));
writetable(A,fullfile(outDir,'frozen_operator_span_audit.csv'));
writetable(F,fullfile(outDir,'fixed_reference_nonlinear_fits.csv'));
fid=fopen(fullfile(outDir,'ANALYSIS_REPORT.md'),'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# EO29 clearance-confounding attribution\n\n');
fprintf(fid,'## Frozen scope\n\n');
fprintf(fid,['Noiseless zero-vibration static data; absolute clearance library; g_ref = %.3f mm; ', ...
    'g_high = %.3f, %.3f, %.3f mm; %d sensors; %d revolutions; fixed spatial common support ', ...
    '(window ratio 0.48); fixed saved trust domain [%.6f, %.6f] mm.\n\n'], ...
    gRef,gRef+.005,gRef+.05,gRef+.3,numel(cfg.alpha_k),cfg.NumRevs_high,trust.domain(1),trust.domain(2));
fprintf(fid,'## Correction of the earlier attribution\n\n');
fprintf(fid,['`per_sensor_clearance_residual.csv` records e=F(g_high,x)-F(g_ref,x) and its tangent projection ', ...
    'onto q=-F_x(g_ref,x). The reported apparent shift is u_app=(q''e)/(q''q), in mm. It is a local ', ...
    'translation-equivalent coordinate, not a physically inferred vibration. The per-sensor tangent ', ...
    'projection explains only the fractions listed in that table (about 1e-8 to 1e-7 here). Therefore ', ...
    '`deprecated_apparent_shift_scores.csv`, including its previous geometry-only ranking, is retained ', ...
    'only as a failed/degenerate diagnostic and is not used for causal attribution.\n\n']);
fprintf(fid,'## Energy-consistent operator ablation\n\n');
fprintf(fid,['`operator_ablation_gamma.csv` keeps the full clearance residual e as the target in every arm. ', ...
    'It first removes the global static-registration direction q=-F_x(x), then reports Gamma=||P_B ', ...
    'M_q e||^2/||M_q e||^2 for each EO. Every row includes rank and condition number of the residualized ', ...
    'two-column EO operator.\n\n']);
fprintf(fid,['A `full`: B=[q sin(k theta),q cos(k theta)]. B `no_intra_passage`: theta=theta_center=2pi(r-1)+alpha_s. ', ...
    'C `no_angular_diversity`: all samples retain x/R but replace alpha_s by sensor-1 alpha. D ', ...
    '`flattened_sensitivity`: q is replaced by a constant rms(q), preserving global EO-column scale while ', ...
    'removing the spatial sensitivity envelope. D is a mathematical counterfactual for separating the ', ...
    'role of F_x shaping; it is not a physical sensor model.\n\n']);
fprintf(fid,['`frozen_operator_span_audit.csv` directly audits the frozen arm. It reports <q_s,e0_s> after ', ...
    'global q residualization and the fraction of B_frozen lying outside the three per-sensor q blocks. ', ...
    'For this shared absolute-template construction, the frozen columns are exactly within that block ', ...
    'span. A per-sensor-static nuisance projection consequently leaves no EO-specific frozen operator ', ...
    'by construction; this is a control, not an additional mechanism claim.\n\n']);
fprintf(fid,['`fixed_reference_nonlinear_fits.csv` is the decisive solver-level check. It compares EO23, EO29 and ', ...
    'EO30 at dg=5, 50 and 300 um with identical samples, dx/amplitude bounds, 36 starts and 800 iterations. ', ...
    'The reference clearance is fixed at %.3f mm: gap is never optimized.\n\n'],gRef);
fprintf(fid,'## Numerical summary\n\n');
for dg=unique(O.delta_g_um).'
    q=F(F.delta_g_um==dg,:);[~,k]=min(q.rmse_V); winner=q.eo(k);
    fprintf(fid,'- dg = %.0f um: nonlinear fixed-reference winner EO%d, RMSE %.6g V.\n',dg,winner,q.rmse_V(k));
    for arm=["full","no_intra_passage","no_angular_diversity","flattened_sensitivity"]
        h=O(O.delta_g_um==dg & O.ablation==arm,:);[~,j]=max(h.gamma);
        full=O(O.delta_g_um==dg & O.ablation=="full" & O.eo==h.eo(j),:);
        fprintf(fid,'  - %s: EO%d Gamma %.6f, delta versus same-EO Full %.6f, rank %d, cond %.4g.\n', ...
            arm,h.eo(j),h.gamma(j),h.gamma(j)-full.gamma,h.B_rank(j),h.B_condition_number(j));
    end
end
fprintf(fid,['\nInterpretation of no-intra-passage: B_frozen=q*[sin(k alpha_s),cos(k alpha_s)] is only a constant ', ...
    'multiple of q within each sensor. In the current shared-template setup, e(x) and q(x) are essentially ', ...
    'the same across sensors, so global q residualization makes <q_s,e0_s> numerically negligible. The ', ...
    'direct block-span audit verifies B_frozen is within the per-sensor q-block span. Eight turns repeat ', ...
    'these samples, but are not the sole or primary explanation for the near-zero frozen projection. ', ...
    'Thus the correct statement is narrower: under this global-q-residualized shared-template linear test, ', ...
    'EO differentiation appears only when intra-passage phase modulation is retained. This does not prove ', ...
    'that passage phase is independently necessary for EO29 in the nonlinear solver or in other layouts. ', ...
    'The no-angular-diversity contrast remains the bounded test of whether installed alpha_s diversity ', ...
    'participates in EO29 relative to EO23/30 advantage.\n']);
fprintf(fid,['\n## Interpretation boundary\n\nEO29 is a risk direction only if it wins a stated layer in this fixed ', ...
    'three-sensor, eight-revolution, support-and-reference configuration. The experiment cannot establish ', ...
    'a universal clearance EO, because the result depends on sensor phases, passage width, reference gap, ', ...
    'static response shape, and the retained support.\n']);
fprintf(fid,'\nEO candidates tested: %s.\n',sprintf('%d ',eoList));
end

function r=sensor_row0(),r=struct('delta_g_um',NaN,'sensor_id',NaN,'n_samples',NaN,'clearance_residual_rms_V',NaN,'apparent_shift_mm',NaN,'translation_explained_fraction',NaN,'translation_residual_rms_V',NaN);end
function r=score_row0(),r=struct('delta_g_um',NaN,'eo',NaN,'layer',"",'explained_fraction',NaN,'rss_baseline_V2',NaN,'rss_after_eo_V2',NaN,'support_samples',NaN);end
function r=operator_row0(),r=struct('delta_g_um',NaN,'eo',NaN,'ablation',"",'operator_definition',"",'gamma',NaN,'residualized_energy_V2',NaN,'unexplained_energy_V2',NaN,'B_rank',NaN,'B_condition_number',NaN,'B_singular_value_1',NaN,'B_singular_value_2',NaN,'support_samples',NaN);end
function r=frozen_audit_row0(),r=struct('delta_g_um',NaN,'eo',NaN,'global_q_residual_energy_V2',NaN,'sensor1_q_dot_global_e0',NaN,'sensor2_q_dot_global_e0',NaN,'sensor3_q_dot_global_e0',NaN,'max_abs_sensor_q_dot_global_e0',NaN,'frozen_B_outside_per_sensor_q_block_fraction',NaN,'per_sensor_q_block_rank',NaN,'per_sensor_residual_energy_V2',NaN,'support_samples',NaN);end
function r=fit_row0(),r=struct('delta_g_um',NaN,'eo',NaN,'g_reference_mm',NaN,'g_free',false,'dx_bound_mm',NaN,'amplitude_bound_mm',NaN,'start_count',NaN,'max_iterations',NaN,'dx_fit_mm',NaN,'amplitude_fit_mm',NaN,'rmse_V',NaN,'sse_V2',NaN,'exitflag',NaN,'iterations',NaN,'support_samples',NaN);end
