function [Detail,Summary]=Run_DualSyncFunnelV1Generalization(snrDb,seeds,conditionIds,forceTrueEO,shortGnTopK,funnelMode,secondJointK,runTag,supportMode,lowTemplateMethod)
%RUN_DUALSYNCFUNNELV1GENERALIZATION Official dual-sync Funnel V2 screen.
% Each raw low/high DAQ record receives iid Gaussian noise before OPR-based
% mapping.  Its standard deviation is defined from the complete fixed-length
% record AC RMS, so 25/15/5 dB is a sensor-acquisition-level SNR.
if nargin<1||isempty(snrDb),snrDb=15;end
if nargin<2||isempty(seeds),seeds=1:10;end
if nargin<4||isempty(forceTrueEO),forceTrueEO=false;end
if nargin<5||isempty(shortGnTopK),shortGnTopK=40;end
validateattributes(shortGnTopK,{'numeric'},{'scalar','integer','>=',3});
if nargin<6||isempty(funnelMode),funnelMode="v2_all_joint1";end
funnelMode=string(funnelMode);if ~ismember(funnelMode,["v1","v2_all_joint1"]),error('Unknown funnelMode.');end
if nargin<7||isempty(secondJointK),secondJointK=24;end
validateattributes(secondJointK,{'numeric'},{'scalar','integer','>=',1});
if nargin<8||isempty(runTag),runTag="";end
runTag=string(runTag);
if nargin<9||isempty(supportMode),supportMode="legacy";end
supportMode=string(supportMode);
if ~ismember(supportMode,["legacy","support_aware"]),error('Unknown supportMode.');end
if nargin<10||isempty(lowTemplateMethod),lowTemplateMethod="adaptive_sg";end
lowTemplateMethod=lower(string(lowTemplateMethod));
if ~ismember(lowTemplateMethod,["adaptive_sg","bin_mean","local_quadratic","spline"])
    error('Unknown lowTemplateMethod: %s',lowTemplateMethod);
end

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','dual_sync_funnel_v2_generalization_daq_snr');
if ~exist(outDir,'dir'),mkdir(outDir);end

conditions=make_conditions();
runSuffix="";
if nargin>=3&&~isempty(conditionIds)
    conditionIds=string(conditionIds);conditionIds=conditionIds(:);
    conditions=conditions(ismember(conditions.condition_id,conditionIds),:);
    if isempty(conditions),error('No requested condition_id was found.');end
    shortIds=erase(conditions.condition_id,["C0","X0","_"]);
    if numel(shortIds)<=4
        runSuffix="_"+strjoin(shortIds,'-');
    else
        runSuffix="_subset"+numel(shortIds)+"_"+shortIds(1)+"_to_"+shortIds(end);
    end
end
if forceTrueEO,runSuffix=string(runSuffix)+"_fixedtrue";end
if shortGnTopK~=40,runSuffix=string(runSuffix)+"_top"+shortGnTopK;end
if funnelMode=="v2_all_joint1",runSuffix=string(runSuffix)+"_v2alljoint1";end
if funnelMode=="v2_all_joint1"&&secondJointK~=24,runSuffix=string(runSuffix)+"_k2"+secondJointK;end
if strlength(runTag)>0,runSuffix=string(runSuffix)+"_"+runTag;end
if supportMode=="support_aware",runSuffix=string(runSuffix)+"_supportaware";end
if lowTemplateMethod~="adaptive_sg",runSuffix=string(runSuffix)+"_"+char(lowTemplateMethod);end
runSuffix=char(runSuffix);

ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;
cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=20;cfg0.NumRevs_high=8;
cfg0.route30ForwardModel="low_increment";cfg0.route30LowTemplateMethod="adaptive_sg";
cfg0.route30LowTemplateOptions=struct('gridSpacingMm',.02,'minBinCount',5,...
    'minCoverage',.90,'numFolds',5,'candidateWindowMm',...
    [.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
cfg0.route30FrequencyStructureMode="dual_sync_sync";cfg0.route30StructuredFastBudgets=false;
cfg0.route30DualFrequencyRangeHz=[500 1300];cfg0.route30GapHalfWidthMm=.20;
cfg0.dualSyncGapCount=11;cfg0.dualSyncKeepPerGap=80;
cfg0.dualSyncIteratedReplayCount=2100;cfg0.dualSyncRefineCount=50;
cfg0.dualSyncTopPairCount=40;cfg0.dualSyncStartsPerTopPair=2;
cfg0.funnelTopResidualPairCount=shortGnTopK-2;cfg0.funnelShortGnPairCount=shortGnTopK;
cfg0.funnelJointGnSteps=2;cfg0.funnelFullRefineCount=3;
cfg0.funnelAllPairsFirstJointStep=(funnelMode=="v2_all_joint1");
cfg0.funnelSecondJointPairCount=secondJointK;
cfg0.structuredComponentAmplitudeFloorMm=.05;
cfg0.returnDualSyncDiagnostics=true;cfg0.dualSyncSearchStrategy="funnel_v1";
if supportMode=="support_aware"
    cfg0.route30SupportAware=true;
    cfg0.route30PhysicalDomain=struct('dx_bounds_mm',[-.03 .03],...
        'amplitude_max_mm',[.25 .25],'gap_bounds_mm',[.40 .70],...
        'source',"Predeclared paper simulation envelope: amplitude <=0.25 mm; high-gap matrix 0.40-0.70 mm; independent dx injection study -0.03 to 0.03 mm");
end

% Use the baseline only to create the common, fixed spatial sampling geometry.
fRot=cfg0.RPM_high/60;baseline=conditions(conditions.condition_id=="C01_baseline",:);
if isempty(baseline)
    baseline=make_conditions();baseline=baseline(baseline.condition_id=="C01_baseline",:);
end
cfgBase=apply_condition(cfg0,baseline,fRot);
truthBase=make_truth_model(lib,cfg0,baseline.g_low_mm);
baselineData=simulate_highspeed_from_low_increment(truthBase,...
    baseline.g_high_mm-baseline.g_low_mm,cfgBase,Inf,'fixed_std',1);

% Calibrate only once per distinct low-speed operating gap.  The low-speed
% SNR is defined before multi-revolution aggregation, on the full raw record.
gLowSet=unique(conditions.g_low_mm,'stable');calBank=cell(numel(gLowSet),1);
calTime=zeros(numel(gLowSet),1);lowSeed=zeros(numel(gLowSet),1);
lowStats=cell(numel(gLowSet),1);
for ig=1:numel(gLowSet)
    lowSeed(ig)=991001+round(1000*gLowSet(ig));
    low=simulate_low_speed_template(@(z)eval_gap_template(lib,gLowSet(ig),z),...
        cfg0,lib.domain,Inf,lowSeed(ig),'fixed_std');
    sLow=rms(low.V_clean-mean(low.V_clean));
    sigmaLow=sLow/10^(snrDb/20);rng(lowSeed(ig),'twister');
    low.V_cap=low.V_clean+sigmaLow*randn(size(low.V_clean));
    lowNoise=low.V_cap-low.V_clean;
    lowStats{ig}=struct('signal_rms_V',sLow,'noise_std_V',sigmaLow,...
        'realized_snr_db',20*log10(sLow/max(std(lowNoise,1),eps)));
    tc=tic;
    if lowTemplateMethod=="adaptive_sg"
        calBank{ig}=calibrate_inv_log_2_low_speed(low,lib,cfg0);
    else
        opt=cfg0.route30LowTemplateOptions;
        opt.mapWindowRatio=low.spatial_window_ratio;
        if lowTemplateMethod=="bin_mean"
            opt.smoothSpan=1;
            lowTemplate=build_low_speed_templates_binned(low,cfg0.alpha_k,cfg0.R_tip,...
                lib.domain,lib.xGrid,'mean',opt);
        elseif lowTemplateMethod=="spline"
            lowTemplate=build_low_speed_templates_adaptive(low,cfg0.alpha_k,cfg0.R_tip,...
                lib.domain,lib.xGrid,'cv_spline',opt);
        else
            lowTemplate=build_low_speed_templates_local_quadratic(low,cfg0.alpha_k,cfg0.R_tip,...
                lib.domain,lib.xGrid,opt);
        end
        calBank{ig}=calibrate_inv_log_2_low_speed(lowTemplate,lib,cfg0);
    end
    calTime(ig)=toc(tc);
end

checkpoint=fullfile(outDir,sprintf('checkpoint_%gdB%s.csv',snrDb,runSuffix));
if exist(checkpoint,'file'),Detail=align_detail_schema(readtable(checkpoint,'TextType','string'));else,Detail=struct2table(repmat(empty_row(),0,1));end

for ic=1:height(conditions)
    c=conditions(ic,:);cfg=apply_condition(cfg0,c,fRot);
    if forceTrueEO,cfg.dualSyncCandidateEOPairs=[c.eo1_true c.eo2_true];end
    daq=make_condition_daq(baselineData,lib,c,cfg,lib.domain);
    sHigh=rms(daq.V_clean-mean(daq.V_clean));
    sigmaHigh=sHigh/10^(snrDb/20);
    ig=find(abs(gLowSet-c.g_low_mm)<1e-12,1);cal=calBank{ig};
    for seed=seeds(:).'
        if is_done(Detail,c.condition_id,seed),continue;end
        % Keep each condition's noise realization invariant to a requested subset.
        highSeed=992000+1000*c.condition_seed_id+seed;rng(highSeed,'twister');
        highNoise=sigmaHigh*randn(size(daq.V_clean));highNoisy=daq;
        highNoisy.V_cap=daq.V_clean+highNoise;
        highMap=map_highspeed_to_space(highNoisy,cfg0.alpha_k,cfg0.R_tip,lib.domain,.02);
        cleanMap=map_highspeed_to_space(daq,cfg0.alpha_k,cfg0.R_tip,lib.domain,.02);
        staticDaq=daq;staticDaq.V_cap=daq.V_static;
        staticMap=map_highspeed_to_space(staticDaq,cfg0.alpha_k,cfg0.R_tip,lib.domain,.02);
        highSNRRealized=20*log10(sHigh/max(std(highNoise,1),eps));
        actualSnr=20*log10(rms(cleanMap.V_a-staticMap.V_a)/sigmaHigh);
        highMap.snr_db_equiv=highSNRRealized;tc=tic;
        fit=run_inv_log_2_high_with_calibration(cal,highMap,cfg);elapsed=toc(tc);
        row=pack_row(c,seed,snrDb,actualSnr,lowSeed(ig),highSeed,fit,cal,cfg,...
            calTime(ig),elapsed,lowStats{ig},sHigh,sigmaHigh,highSNRRealized);
        Detail=[Detail;struct2table(row)]; %#ok<AGROW>
        writetable(Detail,checkpoint);
        fprintf('%s seed=%d EO=(%d,%d) rank=%.0f dA=%.4g ddg=%.4g %.2fs\n',...
            c.condition_id,seed,fit.eo_id(1),fit.eo_id(2),row.true_eo_rank_gn1,...
            row.max_amplitude_error_mm,row.delta_gap_error_mm,elapsed);
    end
end

Summary=summarize_results(Detail,conditions);
writetable(Detail,fullfile(outDir,sprintf('detail_%gdB%s.csv',snrDb,runSuffix)));
writetable(Summary,fullfile(outDir,sprintf('summary_%gdB%s.csv',snrDb,runSuffix)));
save(fullfile(outDir,sprintf('generalization_%gdB%s.mat',snrDb,runSuffix)),...
    'Detail','Summary','conditions','cfg0','snrDb','seeds','lowStats','-v7.3');
plot_summary(Summary,outDir,snrDb,runSuffix);disp(Summary);
end

function C=make_conditions()
id=["C01_baseline";"C02_eo_10_12";"C03_eo_15_18";"C04_eo_10_18";...
    "C05_amp_25_10";"C06_amp_15_10";"C07_amp_20_20";...
    "C08_gap_50_40";"C09_gap_50_70";"C10_gap_70_50";...
    "C11_phase_00";"C12_phase_0_90";"C13_phase_0_180";...
    "C14_amp_18_12";"C15_amp_20_15";"C16_amp_25_20";...
    "X01_corner_near";"X02_corner_lowamp";"X03_corner_equal";"X04_corner_wide"];
kind=["baseline";repmat("eo",3,1);repmat("amplitude",3,1);repmat("gap",3,1);...
    repmat("phase",3,1);repmat("amplitude",3,1);repmat("corner",4,1)];
eo1=[10;10;15;10;10;10;10;10;10;10;10;10;10;10;10;10;10;15;10;10];
eo2=[26;12;18;18;26;26;26;26;26;26;26;26;26;26;26;26;12;18;18;26];
A1=[.25;.25;.25;.25;.25;.15;.20;.25;.25;.25;.25;.25;.25;.18;.20;.25;.25;.15;.20;.25];
A2=[.15;.15;.15;.15;.10;.10;.20;.15;.15;.15;.15;.15;.15;.12;.15;.20;.10;.10;.20;.10];
gL=[.5;.5;.5;.5;.5;.5;.5;.5;.5;.7;.5;.5;.5;.5;.5;.5;.5;.7;.5;.5];
gH=[.6;.6;.6;.6;.6;.6;.6;.4;.7;.5;.6;.6;.6;.6;.6;.6;.7;.5;.4;.7];
p1=[pi/4;pi/4;pi/4;pi/4;pi/4;pi/4;pi/4;pi/4;pi/4;pi/4;0;0;0;pi/4;pi/4;pi/4;0;0;0;0];
p2=[-pi/3;-pi/3;-pi/3;-pi/3;-pi/3;-pi/3;-pi/3;-pi/3;-pi/3;-pi/3;0;pi/2;pi;-pi/3;-pi/3;-pi/3;pi/2;pi;3*pi/2;0];
C=table((1:numel(id)).',id,kind,eo1,eo2,A1,A2,gL,gH,p1,p2,'VariableNames',...
    {'condition_seed_id','condition_id','factor_group','eo1_true','eo2_true','A1_true_mm','A2_true_mm',...
    'g_low_mm','g_high_mm','phi1_true_rad','phi2_true_rad'});
end

function cfg=apply_condition(cfg,c,fRot)
cfg.f_true=[c.eo1_true c.eo2_true]*fRot;
cfg.A_true=[c.A1_true_mm c.A2_true_mm];cfg.phi_true=[c.phi1_true_rad c.phi2_true_rad];
end

function model=make_truth_model(lib,cfg,gLow)
T=repmat(eval_gap_template(lib,gLow,lib.xGrid),1,numel(cfg.alpha_k));
cal=struct('g0',gLow,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,...
    'sensorIds',1:numel(cfg.alpha_k));model=build_path_template_model(lib,T,cal);
end

function D=make_condition_daq(geometry,lib,c,cfg,domain)
% Generate the complete 8-revolution DAQ record before adding measurement noise.
t=geometry.t(:);omega=cfg.RPM_high*2*pi/60;xHalf=.52*diff(domain);
u=zeros(size(t));
for j=1:numel(cfg.A_true)
    u=u+cfg.A_true(j)*sin(2*pi*cfg.f_true(j)*t+cfg.phi_true(j));
end
base=min(eval_gap_template(lib,c.g_low_mm,lib.xGrid));V=base*ones(size(t));Vstatic=base*ones(size(t));
for rev=1:cfg.NumRevs_high
    for sensor=1:numel(cfg.alpha_k)
        te=geometry.T_opr_truth(rev)+cfg.alpha_k(sensor)/omega;
        idx=find(abs(t-te)<=xHalf/(omega*cfg.R_tip));
        if isempty(idx),continue;end
        xNom=omega*cfg.R_tip*(t(idx)-te);q=xNom-u(idx);
        inside=q>=domain(1)&q<=domain(2);
        if any(inside),V(idx(inside))=eval_gap_template(lib,c.g_high_mm,q(inside),sensor);end
        inside0=xNom>=domain(1)&xNom<=domain(2);
        if any(inside0),Vstatic(idx(inside0))=eval_gap_template(lib,c.g_high_mm,xNom(inside0),sensor);end
    end
end
D=geometry;D.V_clean=V;D.V_static=Vstatic;D.V_cap=V;
end

function tf=is_done(T,id,seed)
tf=~isempty(T)&&any(T.condition_id==string(id)&T.seed==seed);
end

function r=empty_row()
r=struct('condition_id',"",'factor_group',"",'seed',NaN,'nominal_snr_db',NaN,...
    'noise_protocol',"complete_daq_waveform_relative_snr_iid",...
    'actual_dynamic_snr_db',NaN,'low_signal_rms_V',NaN,'low_noise_std_V',NaN,...
    'low_realized_snr_db',NaN,'high_signal_rms_V',NaN,'high_noise_std_V',NaN,...
    'high_realized_snr_db',NaN,'low_seed',NaN,'high_seed',NaN,...
    'eo1_true',NaN,'eo2_true',NaN,'eo1_est',NaN,'eo2_est',NaN,'eo_correct',false,...
    'A1_true_mm',NaN,'A2_true_mm',NaN,'A1_est_mm',NaN,'A2_est_mm',NaN,...
    'A1_signed_error_mm',NaN,'A2_signed_error_mm',NaN,...
    'A1_abs_error_mm',NaN,'A2_abs_error_mm',NaN,...
    'max_amplitude_error_mm',NaN,'g_low_true_mm',NaN,'g_high_true_mm',NaN,...
    'g0_est_mm',NaN,'g0_error_mm',NaN,'delta_gap_true_mm',NaN,'delta_gap_est_mm',NaN,...
    'delta_gap_signed_error_mm',NaN,'delta_gap_error_mm',NaN,'rmse_V',NaN,'calibration_time_s',NaN,'fit_time_s',NaN,...
    'coarse_time_s',NaN,'level1_time_s',NaN,'level2_time_s',NaN,...
    'joint1_time_s',NaN,'joint2_time_s',NaN,'full_refine_time_s',NaN,...
    'joint1_pair_count',NaN,'second_joint_pair_count',NaN,...
    'true_eo_rank_gn1',NaN,'true_eo_rank_joint1',NaN,...
    'frequency_margin_V',NaN,'noise_normalized_margin',NaN,...
    'identification_confident',false,'identification_status',"",...
    'competitor_eo1',NaN,'competitor_eo2',NaN,...
    'support_aware',false,'low_support_min_mm',NaN,'low_support_max_mm',NaN,...
    'high_certified_min_mm',NaN,'high_certified_max_mm',NaN,...
    'support_left_margin_mm',NaN,'support_right_margin_mm',NaN,...
    'support_removed_samples',NaN,'support_runtime_valid',false,...
    'recall_at_40',false,'recall_at_3',false,'success',false);
end

function T=align_detail_schema(T)
schema=struct2table(empty_row());names=schema.Properties.VariableNames;
for i=1:numel(names)
    if ~ismember(names{i},T.Properties.VariableNames)
        value=schema.(names{i});T.(names{i})=repmat(value,height(T),1);
    end
end
extra=setdiff(T.Properties.VariableNames,names,'stable');T=T(:,[names extra]);
end

function r=pack_row(c,seed,snr,actualSnr,lowSeed,highSeed,fit,cal,cfg,calTime,elapsed,lowNoise,highRms,highSigma,highSnr)
r=empty_row();r.condition_id=c.condition_id;r.factor_group=c.factor_group;r.seed=seed;
r.nominal_snr_db=snr;r.actual_dynamic_snr_db=actualSnr;
r.low_signal_rms_V=lowNoise.signal_rms_V;r.low_noise_std_V=lowNoise.noise_std_V;
r.low_realized_snr_db=lowNoise.realized_snr_db;r.high_signal_rms_V=highRms;
r.high_noise_std_V=highSigma;r.high_realized_snr_db=highSnr;
r.low_seed=lowSeed;r.high_seed=highSeed;r.eo1_true=c.eo1_true;r.eo2_true=c.eo2_true;
r.eo1_est=fit.eo_id(1);r.eo2_est=fit.eo_id(2);
eoTrue=[c.eo1_true c.eo2_true];r.eo_correct=all(fit.eo_id==eoTrue);
[~,ord]=sort(fit.eo_id);A=fit.A_id(ord);ATrue=[c.A1_true_mm c.A2_true_mm];
r.A1_true_mm=ATrue(1);r.A2_true_mm=ATrue(2);r.A1_est_mm=A(1);r.A2_est_mm=A(2);
r.A1_signed_error_mm=A(1)-ATrue(1);r.A2_signed_error_mm=A(2)-ATrue(2);
r.A1_abs_error_mm=abs(r.A1_signed_error_mm);r.A2_abs_error_mm=abs(r.A2_signed_error_mm);
r.max_amplitude_error_mm=max(abs(A-ATrue));r.g_low_true_mm=c.g_low_mm;
r.g_high_true_mm=c.g_high_mm;r.g0_est_mm=cal.pathCal.g0;r.g0_error_mm=abs(cal.pathCal.g0-c.g_low_mm);
r.delta_gap_true_mm=c.g_high_mm-c.g_low_mm;r.delta_gap_est_mm=fit.g_used-cal.pathCal.g0;
r.delta_gap_signed_error_mm=r.delta_gap_est_mm-r.delta_gap_true_mm;
r.delta_gap_error_mm=abs(r.delta_gap_signed_error_mm);r.rmse_V=fit.rmse;
r.calibration_time_s=calTime;r.fit_time_s=elapsed;r.coarse_time_s=get_field(fit,'coarse_time_s',NaN);
r.level1_time_s=get_field(fit,'exact_gn1_time_s',NaN);r.level2_time_s=get_field(fit,'short_gn_time_s',NaN);
r.joint1_time_s=get_field(fit,'joint1_time_s',NaN);r.joint2_time_s=get_field(fit,'joint2_time_s',NaN);
r.full_refine_time_s=get_field(fit,'refine_time_s',NaN);
r.joint1_pair_count=get_field(fit,'joint1_pair_count',NaN);r.second_joint_pair_count=get_field(fit,'short_gn_pair_count',NaN);
r.frequency_margin_V=get_field(fit,'frequency_margin',NaN);
r.noise_normalized_margin=get_field(fit,'noise_normalized_margin',NaN);
r.identification_confident=get_field(fit,'identification_confident',false);
r.identification_status=string(get_field(fit,'identification_status',"unvalidated"));
if isfield(fit,'support_contract')
    C=fit.support_contract;r.support_aware=true;
    r.low_support_min_mm=C.low_observed_mm(1);r.low_support_max_mm=C.low_observed_mm(2);
    r.high_certified_min_mm=C.high_certified_mm(1);r.high_certified_max_mm=C.high_certified_mm(2);
    r.support_left_margin_mm=C.left_margin_mm;r.support_right_margin_mm=C.right_margin_mm;
    r.support_removed_samples=C.removed_sample_count;r.support_runtime_valid=C.runtime_domain_valid;
end
competitor=get_field(fit,'competing_frequency_or_order',[NaN NaN]);
competitor=sort(competitor(:).');
if numel(competitor)>=2
    r.competitor_eo1=competitor(1)/(cfg.RPM_high/60);
    r.competitor_eo2=competitor(2)/(cfg.RPM_high/60);
end
eo=fit.diagnostic_pair_eo;score=fit.diagnostic_pair_gn1_score;[~,q]=sort(score);
loc=find(all(eo(q,:)==eoTrue,2),1);if ~isempty(loc),r.true_eo_rank_gn1=loc;end
jointEo=get_field(fit,'diagnostic_joint1_eo',zeros(0,2));jointScore=get_field(fit,'diagnostic_joint1_score',zeros(0,1));
if ~isempty(jointEo)
    [~,q]=sort(jointScore);loc=find(all(jointEo(q,:)==eoTrue,2),1);
    if ~isempty(loc),r.true_eo_rank_joint1=loc;end
end
r.recall_at_40=any(all(fit.diagnostic_short_eo==eoTrue,2));
r.recall_at_3=any(all(fit.diagnostic_full_start_eo==eoTrue,2));
r.success=r.eo_correct&&r.max_amplitude_error_mm<=.01&&r.delta_gap_error_mm<=.01;
end

function S=summarize_results(D,C)
n=height(C);row=struct('condition_id',"",'factor_group',"",'n',0,'actual_snr_db',NaN,...
    'eo_correct_rate',NaN,'recall_at_40_rate',NaN,'recall_at_3_rate',NaN,...
    'success_rate',NaN,'confidence_rate',NaN,'ambiguous_rate',NaN,...
    'eo_correct_given_confident_rate',NaN,'median_noise_normalized_margin',NaN,...
    'max_true_eo_rank_gn1',NaN,'max_true_eo_rank_joint1',NaN,'median_amp_error_mm',NaN,...
    'max_amp_error_mm',NaN,'median_delta_gap_error_mm',NaN,'max_delta_gap_error_mm',NaN,...
    'median_fit_time_s',NaN,'p95_fit_time_s',NaN,'nominal_snr_db',NaN,...
    'support_aware_rate',NaN,'support_runtime_valid_rate',NaN,...
    'median_A1_abs_error_mm',NaN,'p95_A1_abs_error_mm',NaN,'A1_est_std_mm',NaN,...
    'median_A2_abs_error_mm',NaN,'p95_A2_abs_error_mm',NaN,'A2_est_std_mm',NaN,...
    'median_delta_gap_signed_error_mm',NaN,'p95_delta_gap_error_mm',NaN);
rows=repmat(row,n,1);
for i=1:n
    q=D(D.condition_id==C.condition_id(i),:);rows(i).condition_id=C.condition_id(i);
    rows(i).factor_group=C.factor_group(i);rows(i).n=height(q);if isempty(q),continue;end
    rows(i).actual_snr_db=median(q.actual_dynamic_snr_db);
    rows(i).nominal_snr_db=median(q.nominal_snr_db);
    rows(i).eo_correct_rate=mean(q.eo_correct);rows(i).recall_at_40_rate=mean(q.recall_at_40);
    rows(i).recall_at_3_rate=mean(q.recall_at_3);rows(i).success_rate=mean(q.success);
    confident=logical(q.identification_confident);
    rows(i).confidence_rate=mean(confident);
    if ismember('support_aware',q.Properties.VariableNames)
        rows(i).support_aware_rate=mean(q.support_aware);
        rows(i).support_runtime_valid_rate=mean(q.support_runtime_valid);
    end
    rows(i).ambiguous_rate=mean(q.identification_status=="ambiguous");
    qc=q(confident,:);
    if ~isempty(qc),rows(i).eo_correct_given_confident_rate=mean(qc.eo_correct);end
    rows(i).median_noise_normalized_margin=median(q.noise_normalized_margin,'omitnan');
    rows(i).max_true_eo_rank_gn1=max(q.true_eo_rank_gn1);
    rows(i).max_true_eo_rank_joint1=max(q.true_eo_rank_joint1);
    rows(i).median_amp_error_mm=median(q.max_amplitude_error_mm);rows(i).max_amp_error_mm=max(q.max_amplitude_error_mm);
    rows(i).median_delta_gap_error_mm=median(q.delta_gap_error_mm);rows(i).max_delta_gap_error_mm=max(q.delta_gap_error_mm);
    rows(i).median_A1_abs_error_mm=median(q.A1_abs_error_mm);rows(i).p95_A1_abs_error_mm=prctile(q.A1_abs_error_mm,95);rows(i).A1_est_std_mm=std(q.A1_est_mm);
    rows(i).median_A2_abs_error_mm=median(q.A2_abs_error_mm);rows(i).p95_A2_abs_error_mm=prctile(q.A2_abs_error_mm,95);rows(i).A2_est_std_mm=std(q.A2_est_mm);
    rows(i).median_delta_gap_signed_error_mm=median(q.delta_gap_signed_error_mm);rows(i).p95_delta_gap_error_mm=prctile(q.delta_gap_error_mm,95);
    rows(i).median_fit_time_s=median(q.fit_time_s);rows(i).p95_fit_time_s=prctile(q.fit_time_s,95);
end
S=struct2table(rows);
end

function plot_summary(S,outDir,snrDb,runSuffix)
f=figure('Color','w','Units','centimeters','Position',[1 1 17 12.5]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');x=1:height(S);labels=S.condition_id;
nexttile;plot(x,S.eo_correct_rate,'o-','LineWidth',1.0);hold on;
plot(x,S.recall_at_40_rate,'s-','LineWidth',1.0);plot(x,S.recall_at_3_rate,'^-','LineWidth',1.0);
ylabel('Rate');ylim([0 1.05]);legend({'Final EO','Recall@40','Recall@3'},'Location','southwest');title('(a) EO retention');
nexttile;semilogy(x,max(S.max_amp_error_mm,1e-6),'o-','LineWidth',1.0);hold on;
semilogy(x,max(S.max_delta_gap_error_mm,1e-6),'s-','LineWidth',1.0);yline(.01,'--','0.01 mm');
ylabel('Maximum error (mm)');legend({'Amplitude','Gap change'},'Location','best');title('(b) Parameter error');
nexttile;bar(x,S.max_true_eo_rank_gn1,'FaceColor',[.25 .50 .75]);hold on;
plot(x,S.max_true_eo_rank_joint1,'ko-','LineWidth',1.0,'MarkerFaceColor','w');yline(24,'--r','K_2=24');
ylabel('Maximum true-EO rank');legend({'Fixed-gap GN1','After joint-GN1','K_2'},'Location','best');title('(c) EO promotion');
nexttile;plot(x,S.median_fit_time_s,'o-','LineWidth',1.0);hold on;plot(x,S.p95_fit_time_s,'s-','LineWidth',1.0);
ylabel('Time (s)');legend({'Median','P95'},'Location','best');title('(d) Runtime');
ax=findall(f,'Type','axes');for k=1:numel(ax),set(ax(k),'FontName','Times New Roman','FontSize',8,...
        'TickDir','in','Box','on','XTick',x,'XTickLabel',labels,'XTickLabelRotation',45);end
xlabel(tl,'Condition','FontName','Times New Roman','FontSize',9);
exportgraphics(f,fullfile(outDir,sprintf('generalization_%gdB%s.png',snrDb,runSuffix)),'Resolution',300);
exportgraphics(f,fullfile(outDir,sprintf('generalization_%gdB%s.pdf',snrDb,runSuffix)),'ContentType','vector');
try,exportgraphics(f,fullfile(outDir,sprintf('generalization_%gdB%s.emf',snrDb,runSuffix)),'ContentType','vector');catch,end
close(f);
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
