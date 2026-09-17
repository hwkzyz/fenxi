function Result=Run_05_Formal20_Predictive_Association()
%RUN_05_FORMAL20_PREDICTIVE_ASSOCIATION Test locked geometry on formal cases.
analysisDir=fileparts(mfilename('fullpath'));programDir=fileparts(analysisDir);
addpath(programDir,'-begin');addpath(fullfile(programDir,'local_func'),'-begin');
addpath(fullfile(analysisDir,'local_func'),'-begin');
[lib,cfg,model,map]=geometry(programDir);conditions=condition_table();
metricRows=repmat(row0(),height(conditions),1);
for i=1:height(conditions)
    metricRows(i)=case_metric(conditions(i,:),map,model,cfg);
end
Metrics=struct2table(metricRows);
formalDir=fullfile(programDir,'low_high_gap_fast_bridge','output', ...
    'dual_sync_funnel_v2_generalization_daq_snr');
levels=[25 15 5]; Combined=table();
for level=levels
    f=fullfile(formalDir,sprintf('summary_%ddB_v2alljoint1_frozen_main_v1_supportaware.csv',level));
    T=readtable(f,'TextType','string');T=T(:,{'condition_id','nominal_snr_db', ...
        'eo_correct_rate','success_rate','recall_at_3_rate','median_amp_error_mm', ...
        'median_delta_gap_error_mm'});
    Combined=[Combined;innerjoin(T,Metrics,'Keys','condition_id')]; %#ok<AGROW>
end
Association=association_table(Combined,levels);
out=fullfile(analysisDir,'output','05_formal20_predictive_association');
if ~exist(out,'dir'),mkdir(out);end
writetable(Metrics,fullfile(out,'locked_geometry_metrics.csv'));
writetable(Combined,fullfile(out,'metrics_with_formal_outcomes.csv'));
writetable(Association,fullfile(out,'association_summary.csv'));
save(fullfile(out,'formal20_association.mat'),'Metrics','Combined','Association');
write_report(out,Association);disp(Association);
Result=struct('Metrics',Metrics,'Combined',Combined,'Association',Association,'outputDir',out);
end

function [lib,cfg,model,map]=geometry(programDir)
ctx=load_inv_log_2_project_context(programDir);tr=load_fixed_trust_domain(programDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN, ...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
cfg=make_inv_log_2_demo_case(ctx,.2).cfgCase;cfg.RPM_low=min(cfg.RPM_high,300);
cfg.NumRevs_low=20;cfg.NumRevs_high=8;cfg.route30ForwardModel="low_increment";
cfg.route30LowTemplateMethod="adaptive_sg";cfg.route30SupportAware=true;
cfg.route30LowTemplateOptions=struct('gridSpacingMm',.02,'minBinCount',5, ...
    'minCoverage',.90,'numFolds',5,'candidateWindowMm',[.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
cfg.route30PhysicalDomain=struct('dx_bounds_mm',[-.03 .03], ...
    'amplitude_max_mm',[.25 .25],'gap_bounds_mm',[.40 .70],'source',"Frozen paper envelope");
low=simulate_low_speed_template(@(z)eval_gap_template(lib,.5,z),cfg,lib.domain,Inf,921001,'fixed_std');
cal=calibrate_inv_log_2_low_speed(low,lib,cfg);model=cal.templateModel;
cfg.A_true=[.25 .15];cfg.f_true=[10 26]*(cfg.RPM_high/60);cfg.phi_true=[pi/4 -pi/3];
d=simulate_highspeed_from_low_increment(model,.1,cfg,Inf,'fixed_std',921101);
m0=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
proxy=struct('x_v',reshape(cal.low_support_observation.per_sensor_observed_mm.',[],1), ...
    'S_v',repelem(cal.low_support_observation.sensor_ids(:),2));
C=derive_support_aware_domain(m0,proxy,model,cfg.route30PhysicalDomain,cfg);
[map,~]=restrict_map_to_support_contract(m0,C);
end

function C=condition_table()
id=["C01_baseline";"C02_eo_10_12";"C03_eo_15_18";"C04_eo_10_18"; ...
"C05_amp_25_10";"C06_amp_15_10";"C07_amp_20_20";"C08_gap_50_40"; ...
"C09_gap_50_70";"C10_gap_70_50";"C11_phase_00";"C12_phase_0_90"; ...
"C13_phase_0_180";"C14_amp_18_12";"C15_amp_20_15";"C16_amp_25_20"; ...
"X01_corner_near";"X02_corner_lowamp";"X03_corner_equal";"X04_corner_wide"];
eo1=[10;10;15;10;10;10;10;10;10;10;10;10;10;10;10;10;10;15;10;10];
eo2=[26;12;18;18;26;26;26;26;26;26;26;26;26;26;26;26;12;18;18;26];
A1=[.25;.25;.25;.25;.25;.15;.20;.25;.25;.25;.25;.25;.25;.18;.20;.25;.25;.15;.20;.25];
A2=[.15;.15;.15;.15;.10;.10;.20;.15;.15;.15;.15;.15;.15;.12;.15;.20;.10;.10;.20;.10];
gH=[.6;.6;.6;.6;.6;.6;.6;.4;.7;.5;.6;.6;.6;.6;.6;.6;.7;.5;.4;.7];
p1=[repmat(pi/4,10,1);0;0;0;pi/4;pi/4;pi/4;0;0;0;0];
p2=[repmat(-pi/3,10,1);0;pi/2;pi;repmat(-pi/3,3,1);pi/2;pi;3*pi/2;0];
C=table(id,eo1,eo2,A1,A2,gH,p1,p2,'VariableNames', ...
    {'condition_id','eo1','eo2','A1','A2','gHigh','phi1','phi2'});
end

function r=case_metric(c,map,model,cfg)
theta=map.theta_v(:);x=map.x_v(:);sid=map.S_v(:);eo=[c.eo1 c.eo2];
u=c.A1*sin(eo(1)*theta+c.phi1)+c.A2*sin(eo(2)*theta+c.phi2);z=x-u;
fx=eval_gap_derivative(model,c.gHigh,z,sid);h=.002;
fg=(eval_gap_template(model,c.gHigh+h,z,sid)-eval_gap_template(model,c.gHigh-h,z,sid))/(2*h);
q=-fx;Jp=[q.*sin(eo(1)*theta),q.*cos(eo(1)*theta),q.*sin(eo(2)*theta),q.*cos(eo(2)*theta)];
Jn=[-fg,q];v=all(isfinite([Jp Jn]),2);Q=orthb(Jn(v,:));J=Jp(v,:)./max(vecnorm(Jp(v,:),2,1),eps);J=J-Q*(Q.'*J);
e=sort(real(eig((J.'*J+(J.'*J).')/2)),'ascend');
[angle,worst]=worst_angle(eo,theta,q,Jn,v);
r=row0();r.condition_id=c.condition_id;r.lambda_min_scaled=max(0,e(1));
r.condition_scaled=e(end)/max(e(1),eps);r.min_competitor_angle_deg=angle;
r.worst_replacement_eo=worst;r.weak_amplitude_mm=min(c.A1,c.A2);
end
function [best,kbest]=worst_angle(eo,t,q,Jn,v)
base=[Jn,q.*sin(eo(1)*t),q.*cos(eo(1)*t)];Qb=orthb(base(v,:));
Xt=[q.*sin(eo(2)*t),q.*cos(eo(2)*t)];Xt=Xt(v,:)-Qb*(Qb.'*Xt(v,:));Qt=orthb(Xt);
best=Inf;kbest=NaN;
for k=10:26
 if any(k==eo),continue;end
 X=[q.*sin(k*t),q.*cos(k*t)];X=X(v,:)-Qb*(Qb.'*X(v,:));Qc=orthb(X);
 if size(Qt,2)<2||size(Qc,2)<2,continue;end
 a=acosd(min(1,norm(Qt.'*Qc,2)));if a<best,best=a;kbest=k;end
end
end
function Q=orthb(A)
[U,S,~]=svd(A,'econ');s=diag(S);r=sum(s>max(size(A))*eps(max(s)));Q=U(:,1:r);
end
function A=association_table(T,levels)
rows=repmat(struct('nominal_snr_db',NaN,'n_conditions',NaN, ...
 'rho_angle_eo_rate',NaN,'rho_info_success_rate',NaN,'rho_weakamp_success_rate',NaN),numel(levels),1);
for i=1:numel(levels)
 q=T(T.nominal_snr_db==levels(i),:);
 rows(i)=struct('nominal_snr_db',levels(i),'n_conditions',height(q), ...
  'rho_angle_eo_rate',corr(q.min_competitor_angle_deg,q.eo_correct_rate,'Type','Spearman'), ...
  'rho_info_success_rate',corr(q.lambda_min_scaled,q.success_rate,'Type','Spearman'), ...
  'rho_weakamp_success_rate',corr(q.weak_amplitude_mm,q.success_rate,'Type','Spearman'));
end
A=struct2table(rows);
end
function write_report(out,A)
fid=fopen(fullfile(out,'README.md'),'w');fprintf(fid,'# Formal-20 predictive association\n\n');
fprintf(fid,'Metrics were locked to the clean exact geometry and then joined to the pre-existing 600-run formal outcomes. This is an out-of-outcome association test, not yet a new blind simulation set.\n\n');
for i=1:height(A)
fprintf(fid,'- %g dB: rho(angle, EO rate)=%.4f; rho(conditioned information, success)=%.4f; rho(weak amplitude, success)=%.4f.\n', ...
A.nominal_snr_db(i),A.rho_angle_eo_rate(i),A.rho_info_success_rate(i),A.rho_weakamp_success_rate(i));end
fprintf(fid,'\nA weak association means the corresponding local metric must not be advertised as a standalone predictor.\n');fclose(fid);
end
function r=row0()
r=struct('condition_id',"",'lambda_min_scaled',NaN,'condition_scaled',NaN, ...
'min_competitor_angle_deg',NaN,'worst_replacement_eo',NaN,'weak_amplitude_mm',NaN);
end
