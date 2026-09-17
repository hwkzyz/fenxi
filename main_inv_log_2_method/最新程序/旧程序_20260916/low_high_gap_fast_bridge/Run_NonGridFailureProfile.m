function Result=Run_NonGridFailureProfile()
% Fixed-frequency voltage profile for the reproducible 733+1217 failure.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);trust=load_fixed_trust_domain(mainDir);lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),trust);base=make_inv_log_2_demo_case(ctx,.2);c=base.cfgCase;fTrue=[733 1217];
d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,.2,z),c.RPM_high,c.NumRevs_high,c.fs,c.R_tip,c.alpha_k,0,lib.domain,[.25 .15],fTrue,[pi/4 -pi/3],'noise_ratio');sig=rms(d.V_clean-min(d.V_clean));rng(990101,'twister');d.V_cap=d.V_clean+sig/10*randn(size(d.V_clean));m=map_highspeed_to_space(d,c.alpha_k,c.R_tip,lib.domain,.02);m=subset_map(m,900);
pairs=[733 1217;733 1066.6;725 1225;737.5 1212.5];rows=repmat(row0(),0,1);
for i=1:size(pairs,1),tc=tic;[q,s]=fit_pair(pairs(i,:),m,lib);r=row0();r.f1=pairs(i,1);r.f2=pairs(i,2);r.g=q(1);r.dx=q(2);r.A1=hypot(q(3),q(4));r.A2=hypot(q(5),q(6));r.rmse=sqrt(s);r.elapsed_s=toc(tc);rows(end+1,1)=r;end
Audit=struct2table(rows);out=fullfile(root,'output','non_grid_profile');if ~exist(out,'dir'),mkdir(out);end;writetable(Audit,fullfile(out,'non_grid_failure_profile.csv'));save(fullfile(out,'non_grid_failure_profile.mat'),'Audit','-v7.3');Result=struct('Audit',Audit,'outputDir',out);
end
function [bestQ,bestS]=fit_pair(f,m,lib)
t=m.t_v(:);x=m.x_v(:);V=m.V_a(:);opts=optimset('Display','off','MaxIter',300,'TolX',1e-8,'TolFun',1e-13);bestS=inf;seeds=[.2 0 .18 .18 .12 .12;.2 0 .25 .25 .15 .15;.3 0 .18 .18 .12 .12];
for k=1:size(seeds,1),obj=@(z)mean((V-eval_gap_template(lib,z(1),x-z(2)-z(3)*sin(2*pi*f(1)*t)-z(4)*cos(2*pi*f(1)*t)-z(5)*sin(2*pi*f(2)*t)-z(6)*cos(2*pi*f(2)*t))).^2)+1e3*(max(.05-z(1),0)^2+max(z(1)-1.5,0)^2+max(abs(z(2))-.5,0)^2);q=fminsearch(obj,seeds(k,:),opts);s=obj(q);if s<bestS,bestS=s;bestQ=q;end;end
end
function m=subset_map(m,N),if numel(m.t_v)<=N,return;end;n=numel(m.t_v);ii=round(linspace(1,n,N));fn=fieldnames(m);for k=1:numel(fn),v=m.(fn{k});if isnumeric(v)&&isvector(v)&&numel(v)==n,m.(fn{k})=v(ii);end;end;end
function r=row0(),r=struct('f1',NaN,'f2',NaN,'g',NaN,'dx',NaN,'A1',NaN,'A2',NaN,'rmse',NaN,'elapsed_s',NaN);end
