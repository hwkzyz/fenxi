function Result = Run_02_ExactFrozen_ConditionalGeometry()
%RUN_02_EXACTFROZEN_CONDITIONALGEOMETRY Analyze the frozen formal geometry.
% This study does not call or modify the identification solver.  It uses the
% maintained static library, clean adaptive-SG low-speed calibration, the
% support-aware contract, and the formal three-probe/eight-revolution map.

analysisDir = fileparts(mfilename('fullpath'));
programDir = fileparts(analysisDir);
addpath(programDir,'-begin');
addpath(fullfile(programDir,'local_func'),'-begin');
addpath(fullfile(analysisDir,'local_func'),'-begin');

[lib,cfg,model,map,contract] = make_frozen_geometry(programDir);
fRot = cfg.RPM_high/60;
cases = [
    struct('id',"regular_10_26",'eo',[10 26],'A',[.25 .15],'phi',[pi/4 -pi/3])
    struct('id',"X01_10_12",'eo',[10 12],'A',[.25 .10],'phi',[0 pi/2])
    struct('id',"X02_15_18",'eo',[15 18],'A',[.15 .10],'phi',[0 pi])];

rows = repmat(empty_case_row(),0,1);
r2rows = repmat(empty_r2_row(),0,1);
scanRows = repmat(empty_scan_row(),0,1);
for i = 1:numel(cases)
    c = cases(i);
    [caseRows,caseR2,caseScan] = analyze_case(c,map,model,cfg,fRot);
    rows = [rows;caseRows]; %#ok<AGROW>
    r2rows = [r2rows;caseR2]; %#ok<AGROW>
    scanRows = [scanRows;caseScan]; %#ok<AGROW>
end

CaseMetrics = struct2table(rows);
PhaseMetrics = struct2table(r2rows);
CompetitorScan = struct2table(scanRows);
Gate = evaluate_gate(CaseMetrics);

outDir = fullfile(analysisDir,'output','02_exact_frozen_conditional_geometry');
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(CaseMetrics,fullfile(outDir,'case_metrics.csv'));
writetable(PhaseMetrics,fullfile(outDir,'phase_metrics.csv'));
writetable(CompetitorScan,fullfile(outDir,'competitor_scan.csv'));
writetable(struct2table(Gate),fullfile(outDir,'gate_assessment.csv'));
save(fullfile(outDir,'exact_frozen_conditional_geometry.mat'), ...
    'CaseMetrics','PhaseMetrics','CompetitorScan','Gate','contract','cfg', ...
    'map','model','cases','-v7.3');
write_report(outDir,CaseMetrics,Gate);
disp(CaseMetrics);
disp(struct2table(Gate));
Result = struct('CaseMetrics',CaseMetrics,'PhaseMetrics',PhaseMetrics, ...
    'CompetitorScan',CompetitorScan,'Gate',Gate,'outputDir',outDir);
end

function [lib,cfg,model,map,contract] = make_frozen_geometry(programDir)
ctx = load_inv_log_2_project_context(programDir);
tr = load_fixed_trust_domain(programDir);
lib = make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell, ...
    NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base = make_inv_log_2_demo_case(ctx,.2);
cfg = base.cfgCase;
cfg.RPM_low = min(cfg.RPM_high,300);
cfg.NumRevs_low = 20;
cfg.NumRevs_high = 8;
cfg.route30ForwardModel = "low_increment";
cfg.route30LowTemplateMethod = "adaptive_sg";
cfg.route30SupportAware = true;
cfg.route30LowTemplateOptions = struct('gridSpacingMm',.02,'minBinCount',5, ...
    'minCoverage',.90,'numFolds',5,'candidateWindowMm', ...
    [.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
cfg.route30PhysicalDomain = struct('dx_bounds_mm',[-.03 .03], ...
    'amplitude_max_mm',[.25 .25],'gap_bounds_mm',[.40 .70], ...
    'source',"Frozen paper simulation envelope");

gLow = .5; gHigh = .6;
low = simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z), ...
    cfg,lib.domain,Inf,911001,'fixed_std');
cal = calibrate_inv_log_2_low_speed(low,lib,cfg);
model = cal.templateModel;

cfg.A_true = [.25 .15];
cfg.f_true = [10 26]*(cfg.RPM_high/60);
cfg.phi_true = [pi/4 -pi/3];
daq = simulate_highspeed_from_low_increment(model,gHigh-gLow,cfg,Inf,'fixed_std',911101);
map0 = map_highspeed_to_space(daq,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
lowProxy = struct('x_v',reshape(cal.low_support_observation.per_sensor_observed_mm.',[],1), ...
    'S_v',repelem(cal.low_support_observation.sensor_ids(:),2));
contract = derive_support_aware_domain(map0,lowProxy,model,cfg.route30PhysicalDomain,cfg);
[map,contract] = restrict_map_to_support_contract(map0,contract);
end

function [rows,r2rows,scanRows] = analyze_case(c,map,model,cfg,~)
theta = map.theta_v(:); x = map.x_v(:); sid = map.S_v(:);
thetaCenter = theta-x/cfg.R_tip;
rows = repmat(empty_case_row(),2,1);
r2rows = repmat(empty_r2_row(),0,1);
scanRows = repmat(empty_scan_row(),0,1);

for modeIndex = 1:2
    if modeIndex == 1, mode = "exact"; phaseTheta = theta;
    else, mode = "frozen"; phaseTheta = thetaCenter; end
    u = harmonic_u(c,phaseTheta);
    z = x-u;
    fx = eval_gap_derivative(model,.6,z,sid);
    h = .002;
    fg = (eval_gap_template(model,.6+h,z,sid)-eval_gap_template(model,.6-h,z,sid))/(2*h);
    q = -fx;
    Jv = [q.*sin(c.eo(1)*phaseTheta),q.*cos(c.eo(1)*phaseTheta), ...
          q.*sin(c.eo(2)*phaseTheta),q.*cos(c.eo(2)*phaseTheta)];
    Jn = [-fg,q];
    valid = all(isfinite([Jv Jn]),2);
    info = conditional_information_metrics(Jv(valid,:),Jn(valid,:));

    [minAngle,worstEO,scan] = scan_replacement_competitors( ...
        c,phaseTheta,q,Jn,valid,mode,info.rank_conditioned==size(Jv,2));
    for j = 1:numel(scan)
        scanRows(end+1,1) = scan(j); %#ok<AGROW>
        scanRows(end).case_id = c.id;
    end

    rows(modeIndex) = pack_case(c,mode,info,minAngle,worstEO,nnz(valid));
    for component = 1:2
        for sensor = [0 unique(sid(:).')]
            if sensor == 0, use = valid; label = "all";
            else, use = valid & sid==sensor; label = "sensor_"+sensor; end
            w = q(use).^2;
            pm = compute_phase_diversity_metrics(c.eo(component)*phaseTheta(use),w);
            rr = empty_r2_row(); rr.case_id=c.id; rr.mode=mode;
            rr.component=component; rr.eo=c.eo(component); rr.sensor=label;
            rr.sample_count=pm.sample_count; rr.weight_sum=pm.weight_sum;
            rr.r2=pm.r2; rr.lambda_min_gram=pm.lambda_min_gram;
            rr.condition_gram=pm.condition_gram;
            r2rows(end+1,1)=rr; %#ok<AGROW>
        end
    end
end
end

function u = harmonic_u(c,phaseTheta)
u = c.A(1)*sin(c.eo(1)*phaseTheta+c.phi(1)) + ...
    c.A(2)*sin(c.eo(2)*phaseTheta+c.phi(2));
end

function M = conditional_information_metrics(Jp,Jn)
[Qn,rankN] = orth_basis(Jn);
Jc = Jp-Qn*(Qn.'*Jp);
Fraw = Jc.'*Jc;
scale = vecnorm(Jp,2,1); scale(scale<=eps)=1;
Jps = Jp./scale;
Jcs = Jps-Qn*(Qn.'*Jps);
Fnorm = Jcs.'*Jcs;
eRaw = sort(real(eig((Fraw+Fraw.')/2)),'ascend');
eNorm = sort(real(eig((Fnorm+Fnorm.')/2)),'ascend');
M = struct('lambda_min_raw',max(0,eRaw(1)), ...
    'lambda_min_scaled',max(0,eNorm(1)), ...
    'condition_scaled',eNorm(end)/max(eNorm(1),eps), ...
    'trace_scaled',trace(Fnorm), ...
    'nuisance_absorption',1-norm(Jc,'fro')^2/max(norm(Jp,'fro')^2,eps), ...
    'rank_nuisance',rankN,'rank_conditioned',rank(Jc,1e-10*norm(Jc,2)));
end

function [minAngle,worstEO,rows] = scan_replacement_competitors(c,phaseTheta,q,Jn,valid,mode,fullModelRank)
sharedEO = c.eo(1); trueEO = c.eo(2); eoGrid = 10:26;
shared = [q.*sin(sharedEO*phaseTheta),q.*cos(sharedEO*phaseTheta)];
base = [Jn shared]; base=base(valid,:); [Qbase,~]=orth_basis(base);
Xt = [q.*sin(trueEO*phaseTheta),q.*cos(trueEO*phaseTheta)]; Xt=Xt(valid,:);
Xt = Xt-Qbase*(Qbase.'*Xt); [Qt,rankTrue]=orth_basis(Xt);
rows = repmat(empty_scan_row(),0,1); minAngle=Inf; worstEO=NaN;
for k = eoGrid
    if k==sharedEO || k==trueEO, continue; end
    Xc=[q.*sin(k*phaseTheta),q.*cos(k*phaseTheta)]; Xc=Xc(valid,:);
    Xc=Xc-Qbase*(Qbase.'*Xc); [Qc,rankComp]=orth_basis(Xc);
    % A principal-angle comparison is not interpretable when the complete
    % four-parameter vibration tangent space is already rank deficient.
    if ~fullModelRank || rankTrue<2 || rankComp<2
        rho=NaN; angle=NaN;
    else
        rho=min(1,norm(Qt.'*Qc,2)); angle=acosd(rho);
    end
    r=empty_scan_row(); r.mode=mode; r.true_pair=string(mat2str(c.eo));
    r.competitor_pair=string(mat2str([sharedEO k])); r.replacement_eo=k;
    r.canonical_correlation=rho; r.min_principal_angle_deg=angle;
    r.true_conditional_rank=rankTrue; r.competitor_conditional_rank=rankComp;
    rows(end+1,1)=r; %#ok<AGROW>
    if isfinite(angle) && angle<minAngle, minAngle=angle; worstEO=k; end
end
if ~isfinite(minAngle), minAngle=NaN; end
end

function [Q,r] = orth_basis(A)
if isempty(A),Q=zeros(size(A,1),0);r=0;return;end
[U,S,~]=svd(A,'econ');s=diag(S);
if isempty(s),r=0;else,r=sum(s>max(size(A))*eps(max(s)));end
Q=U(:,1:r);
end

function r=pack_case(c,mode,m,a,k,n)
r=empty_case_row(); r.case_id=c.id; r.mode=mode;
r.true_pair=string(mat2str(c.eo)); r.valid_samples=n;
r.lambda_min_raw=m.lambda_min_raw; r.lambda_min_scaled=m.lambda_min_scaled;
r.condition_scaled=m.condition_scaled; r.trace_scaled=m.trace_scaled;
r.nuisance_absorption=m.nuisance_absorption;
r.rank_nuisance=m.rank_nuisance; r.rank_conditioned=m.rank_conditioned;
r.worst_replacement_eo=k; r.min_competitor_angle_deg=a;
end

function G=evaluate_gate(T)
ids=unique(T.case_id,'stable'); gainInfo=nan(numel(ids),1); exactAngle=gainInfo;
rankGain=false(numel(ids),1);
for i=1:numel(ids)
    q=T(T.case_id==ids(i),:); e=q(q.mode=="exact",:); f=q(q.mode=="frozen",:);
    gainInfo(i)=e.lambda_min_scaled/max(f.lambda_min_scaled,eps);
    exactAngle(i)=e.min_competitor_angle_deg;
    rankGain(i)=e.rank_conditioned>f.rank_conditioned;
end
G=struct('gate_name',"Gate_1_initial_geometry", ...
    'case_count',numel(ids),'median_conditioned_information_ratio',median(gainInfo), ...
    'max_conditioned_information_ratio',max(gainInfo), ...
    'rank_recovery_case_count',sum(rankGain), ...
    'median_exact_competitor_angle_deg',median(exactAngle,'omitnan'), ...
    'minimum_exact_competitor_angle_deg',min(exactAngle,[],'omitnan'), ...
    'supports_trajectory_main_claim', ...
        all(rankGain) && all(gainInfo>1.10));
end

function write_report(outDir,T,G)
fid=fopen(fullfile(outDir,'README.md'),'w');
fprintf(fid,'# Exact/Frozen Conditional Geometry\n\n');
fprintf(fid,'This is a local information-geometry audit, not an identification success-rate test.\n\n');
fprintf(fid,'- Median exact/frozen conditioned-information ratio: %.6g\n',G.median_conditioned_information_ratio);
fprintf(fid,'- Maximum exact/frozen conditioned-information ratio: %.6g\n',G.max_conditioned_information_ratio);
fprintf(fid,'- Cases with conditioned-rank recovery: %d/%d\n',G.rank_recovery_case_count,G.case_count);
fprintf(fid,'- Median exact-model minimum competitor angle: %.6g deg\n',G.median_exact_competitor_angle_deg);
fprintf(fid,'- Smallest exact-model competitor angle: %.6g deg\n',G.minimum_exact_competitor_angle_deg);
fprintf(fid,'- Frozen-model angles are intentionally undefined because its full vibration tangent space is rank deficient.\n');
fprintf(fid,'- Initial Gate 1 support: %d\n\n',G.supports_trajectory_main_claim);
fprintf(fid,'The gate is deliberately provisional. A positive result justifies matched window and model-distance studies; it does not prove global identifiability or a success-rate advantage.\n\n');
fprintf(fid,'Cases: %s.\n',strjoin(cellstr(unique(T.case_id,'stable')),', '));
fclose(fid);
end

function r=empty_case_row()
r=struct('case_id',"",'mode',"",'true_pair',"",'valid_samples',NaN, ...
    'lambda_min_raw',NaN,'lambda_min_scaled',NaN,'condition_scaled',NaN, ...
    'trace_scaled',NaN,'nuisance_absorption',NaN, ...
    'rank_nuisance',NaN,'rank_conditioned',NaN, ...
    'worst_replacement_eo',NaN,'min_competitor_angle_deg',NaN);
end
function r=empty_r2_row()
r=struct('case_id',"",'mode',"",'component',NaN,'eo',NaN,'sensor',"", ...
    'sample_count',NaN,'weight_sum',NaN,'r2',NaN, ...
    'lambda_min_gram',NaN,'condition_gram',NaN);
end
function r=empty_scan_row()
r=struct('case_id',"",'mode',"",'true_pair',"",'competitor_pair',"", ...
    'replacement_eo',NaN,'canonical_correlation',NaN, ...
    'min_principal_angle_deg',NaN,'true_conditional_rank',NaN, ...
    'competitor_conditional_rank',NaN);
end
