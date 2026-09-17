function Result=Run_02_ConditionalCRLBMonteCarlo(snrList,seeds)
% Fixed-correct-EO local nonlinear least-squares CRLB validation.
if nargin<1||isempty(snrList),snrList=[25 15 5];end
if nargin<2||isempty(seeds),seeds=1:60;end
root=fileparts(mfilename('fullpath'));programDir=fileparts(root);
addpath(programDir,'-begin');addpath(fullfile(programDir,'local_func'),'-begin');
src=fullfile(programDir,'identifiability_mechanism_analysis','output',...
    '02_exact_frozen_conditional_geometry','exact_frozen_conditional_geometry.mat');
S=load(src,'model','cfg','contract');model=S.model;cfg=S.cfg;contract=S.contract;
c=struct('eo',[10 26],'A',[.25 .15],'phi',[pi/4 -pi/3],'gHigh',.6,'deltaGap',.1);
cfg.f_true=c.eo*(cfg.RPM_high/60);cfg.A_true=c.A;cfg.phi_true=c.phi;
cfg.dualSyncCandidateEOPairs=c.eo;cfg.funnelAllPairsFirstJointStep=true;
cfg.funnelSecondJointPairCount=1;cfg.funnelFullRefineCount=1;
cfg.dualSyncMaxIter=600;cfg.structuredComponentAmplitudeFloorMm=.05;
clean=simulate_highspeed_from_low_increment(model,c.deltaGap,cfg,Inf,'fixed_std',930001);
mapClean=map_highspeed_to_space(clean,cfg.alpha_k,cfg.R_tip,model.domain,.02);
[mapClean,~]=restrict_map_to_support_contract(mapClean,contract);
truth=true_parameters(mapClean,c,cfg);rows=repmat(row0(),0,1);predRows=repmat(pred0(),0,1);
Local=make_local_linearization(mapClean,model,truth);
for is=1:numel(snrList)
    snrDb=snrList(is);probe=simulate_highspeed_from_low_increment(model,c.deltaGap,cfg,snrDb,'snr_db',930100+is);
    [CovTarget,Pred]=conditional_crlb(mapClean,model,truth,probe.noise_std);
    predRows(end+1,1)=pack_pred(snrDb,probe.noise_std,CovTarget,Pred,truth); %#ok<AGROW>
    zTruth=[truth.g truth.dx truth.coef];
    rng(940000+1000*is,'twister');
    for k=1:numel(seeds)
        seed=seeds(k);map=mapClean;map.V_a=mapClean.V_a+probe.noise_std*randn(size(mapClean.V_a));
        fit=fixed_pair_local_fit(map,Local,zTruth);
        rr=row0();rr.snr_db=snrDb;rr.seed=seed;rr.noise_std_V=probe.noise_std;
        rr.g_est_mm=fit.theta(1);rr.dx_est_mm=fit.theta(2);
        rr.a1_est_mm=fit.theta(3);rr.b1_est_mm=fit.theta(4);
        rr.a2_est_mm=fit.theta(5);rr.b2_est_mm=fit.theta(6);
        rr.A1_est_mm=fit.A(1);rr.A2_est_mm=fit.A(2);
        rr.phi1_error_deg=wrap180(rad2deg(fit.phi(1)-truth.phi(1)));
        rr.phi2_error_deg=wrap180(rad2deg(fit.phi(2)-truth.phi(2)));
        rr.rmse_V=fit.rmse;rr.converged=isfinite(fit.rmse)&&all(isfinite(fit.theta));
        rows(end+1,1)=rr; %#ok<AGROW>
    end
end

function fit=fixed_pair_local_fit(map,L,z0)
V=map.V_a(:);valid=L.valid;delta=L.J(valid,:)\(V(valid)-L.V0(valid));z=z0(:).'+delta.';
A=[hypot(z(3),z(4)) hypot(z(5),z(6))];phi=[atan2(z(4),z(3)) atan2(z(6),z(5))];
r=V(valid)-(L.V0(valid)+L.J(valid,:)*delta);
fit=struct('theta',z,'A',A,'phi',phi,'rmse',sqrt(mean(r.^2)),...
    'solve_info',struct('iterations',1,'method',"truth-local linearized least squares"));
end
function L=make_local_linearization(map,model,t)
theta=map.theta_v(:);x=map.x_v(:);sid=map.S_v(:);
u=t.coef(1)*sin(t.eo(1)*theta)+t.coef(2)*cos(t.eo(1)*theta)+...
    t.coef(3)*sin(t.eo(2)*theta)+t.coef(4)*cos(t.eo(2)*theta);xq=x-t.dx-u;
V0=eval_gap_template(model,t.g,xq,sid);fx=eval_gap_derivative(model,t.g,xq,sid);h=.002;
fg=(eval_gap_template(model,t.g+h,xq,sid)-eval_gap_template(model,t.g-h,xq,sid))/(2*h);
J=[fg,-fx,-fx.*sin(t.eo(1)*theta),-fx.*cos(t.eo(1)*theta),...
    -fx.*sin(t.eo(2)*theta),-fx.*cos(t.eo(2)*theta)];
L=struct('V0',V0,'J',J,'valid',all(isfinite([V0 J]),2));
end
Detail=struct2table(rows);Prediction=struct2table(predRows);Summary=summarize_mc(Detail,Prediction,truth);
out=fullfile(root,'output','02_conditional_crlb_monte_carlo');if ~exist(out,'dir'),mkdir(out);end
writetable(Detail,fullfile(out,'detail.csv'));writetable(Prediction,fullfile(out,'crlb_prediction.csv'));
writetable(Summary,fullfile(out,'summary.csv'));save(fullfile(out,'crlb_mc.mat'),...
    'Detail','Prediction','Summary','truth','snrList','seeds','-v7.3');write_report(out,Summary);
disp(Summary);Result=struct('Detail',Detail,'Prediction',Prediction,'Summary',Summary,'outputDir',out);
end

function truth=true_parameters(map,c,cfg)
thetaTime=2*pi*(cfg.RPM_high/60)*map.t_v(:);offset=angle(mean(exp(1i*(map.theta_v(:)-thetaTime))));
phi=c.phi-c.eo*offset;truth=struct('g',c.gHigh,'dx',0,'phi',phi,...
    'coef',[c.A(1)*cos(phi(1)),c.A(1)*sin(phi(1)),c.A(2)*cos(phi(2)),c.A(2)*sin(phi(2))],...
    'A',c.A,'eo',c.eo);
end
function [C,P]=conditional_crlb(map,model,t,sigma)
theta=map.theta_v(:);x=map.x_v(:);sid=map.S_v(:);u=t.coef(1)*sin(t.eo(1)*theta)+...
    t.coef(2)*cos(t.eo(1)*theta)+t.coef(3)*sin(t.eo(2)*theta)+t.coef(4)*cos(t.eo(2)*theta);
z=x-t.dx-u;fx=eval_gap_derivative(model,t.g,z,sid);h=.002;
fg=(eval_gap_template(model,t.g+h,z,sid)-eval_gap_template(model,t.g-h,z,sid))/(2*h);
Jp=[fg,-fx.*sin(t.eo(1)*theta),-fx.*cos(t.eo(1)*theta),...
    -fx.*sin(t.eo(2)*theta),-fx.*cos(t.eo(2)*theta)];Jn=-fx;
valid=all(isfinite([Jp Jn]),2);Jp=Jp(valid,:)/sigma;Jn=Jn(valid,:)/sigma;
Q=orthb(Jn);Jc=Jp-Q*(Q.'*Jp);F=Jc.'*Jc;C=pinv((F+F.')/2,1e-12);
P=struct('rank',rank(Jc,1e-10*norm(Jc,2)),'condition',cond(F),...
    'lambda_min',min(real(eig((F+F.')/2))),'sample_count',nnz(valid));
end
function Q=orthb(A)
[U,S,~]=svd(A,'econ');s=diag(S);r=sum(s>1e-10*max(s));Q=U(:,1:r);
end
function r=pack_pred(snr,noise,C,P,t)
r=pred0();r.snr_db=snr;r.noise_std_V=noise;r.rank=P.rank;r.condition=P.condition;
r.lambda_min=P.lambda_min;r.sample_count=P.sample_count;r.g_std_crlb_mm=sqrt(C(1,1));
r.a1_std_crlb_mm=sqrt(C(2,2));r.b1_std_crlb_mm=sqrt(C(3,3));
r.a2_std_crlb_mm=sqrt(C(4,4));r.b2_std_crlb_mm=sqrt(C(5,5));
[r.A1_std_crlb_mm,r.phi1_std_crlb_deg]=polar_std(C(2:3,2:3),t.coef(1:2));
[r.A2_std_crlb_mm,r.phi2_std_crlb_deg]=polar_std(C(4:5,4:5),t.coef(3:4));
end
function [sA,sPhi]=polar_std(C,c)
A=hypot(c(1),c(2));gA=c(:)/A;gp=[-c(2);c(1)]/A^2;
sA=sqrt(max(0,gA.'*C*gA));sPhi=rad2deg(sqrt(max(0,gp.'*C*gp)));
end
function S=summarize_mc(D,P,t)
rows=repmat(sum0(),height(P),1);
for i=1:height(P)
    q=D(D.snr_db==P.snr_db(i)&D.converged,:);r=sum0();r.snr_db=P.snr_db(i);r.n=height(q);
    r.g_std_emp_mm=std(q.g_est_mm);r.g_std_crlb_mm=P.g_std_crlb_mm(i);r.g_ratio=r.g_std_emp_mm/r.g_std_crlb_mm;
    r.A1_std_emp_mm=std(q.A1_est_mm);r.A1_std_crlb_mm=P.A1_std_crlb_mm(i);r.A1_ratio=r.A1_std_emp_mm/r.A1_std_crlb_mm;
    r.A2_std_emp_mm=std(q.A2_est_mm);r.A2_std_crlb_mm=P.A2_std_crlb_mm(i);r.A2_ratio=r.A2_std_emp_mm/r.A2_std_crlb_mm;
    r.phi1_std_emp_deg=std(q.phi1_error_deg);r.phi1_std_crlb_deg=P.phi1_std_crlb_deg(i);r.phi1_ratio=r.phi1_std_emp_deg/r.phi1_std_crlb_deg;
    r.phi2_std_emp_deg=std(q.phi2_error_deg);r.phi2_std_crlb_deg=P.phi2_std_crlb_deg(i);r.phi2_ratio=r.phi2_std_emp_deg/r.phi2_std_crlb_deg;
    r.g_bias_mm=mean(q.g_est_mm-t.g);r.A1_bias_mm=mean(q.A1_est_mm-t.A(1));r.A2_bias_mm=mean(q.A2_est_mm-t.A(2));
    ratios=[r.g_ratio r.A1_ratio r.A2_ratio r.phi1_ratio r.phi2_ratio];
    r.median_std_ratio=median(ratios);r.all_same_order=all(ratios>.33&ratios<3);rows(i)=r;
end
S=struct2table(rows);
end
function write_report(out,S)
fid=fopen(fullfile(out,'README.md'),'w');fprintf(fid,'# Conditional CRLB--Monte Carlo\n\n');
fprintf(fid,'The correct EO pair is fixed and a bounded local nonlinear least-squares fit is initialized at the truth. Generator and inverse share the same frozen forward model. Noise is injected on the complete DAQ time record before OPR mapping.\n\n');
for i=1:height(S),fprintf(fid,'- %g dB: n=%g, median empirical/CRLB std ratio %.6g, all five quantities within [1/3,3]: %d.\n',...
    S.snr_db(i),S.n(i),S.median_std_ratio(i),S.all_same_order(i));end
fprintf(fid,'\nAgreement supports conditional local recoverability only; it does not address EO selection or model discrepancy.\n');fclose(fid);
end
function y=wrap180(x),y=mod(x+180,360)-180;end
function r=row0()
r=struct('snr_db',NaN,'seed',NaN,'noise_std_V',NaN,'g_est_mm',NaN,'dx_est_mm',NaN,...
    'a1_est_mm',NaN,'b1_est_mm',NaN,'a2_est_mm',NaN,'b2_est_mm',NaN,...
    'A1_est_mm',NaN,'A2_est_mm',NaN,'phi1_error_deg',NaN,'phi2_error_deg',NaN,...
    'rmse_V',NaN,'converged',false);
end
function r=pred0()
r=struct('snr_db',NaN,'noise_std_V',NaN,'rank',NaN,'condition',NaN,'lambda_min',NaN,...
    'sample_count',NaN,'g_std_crlb_mm',NaN,'a1_std_crlb_mm',NaN,'b1_std_crlb_mm',NaN,...
    'a2_std_crlb_mm',NaN,'b2_std_crlb_mm',NaN,'A1_std_crlb_mm',NaN,...
    'A2_std_crlb_mm',NaN,'phi1_std_crlb_deg',NaN,'phi2_std_crlb_deg',NaN);
end
function r=sum0()
r=struct('snr_db',NaN,'n',NaN,'g_std_emp_mm',NaN,'g_std_crlb_mm',NaN,'g_ratio',NaN,...
    'A1_std_emp_mm',NaN,'A1_std_crlb_mm',NaN,'A1_ratio',NaN,'A2_std_emp_mm',NaN,...
    'A2_std_crlb_mm',NaN,'A2_ratio',NaN,'phi1_std_emp_deg',NaN,'phi1_std_crlb_deg',NaN,...
    'phi1_ratio',NaN,'phi2_std_emp_deg',NaN,'phi2_std_crlb_deg',NaN,'phi2_ratio',NaN,...
    'g_bias_mm',NaN,'A1_bias_mm',NaN,'A2_bias_mm',NaN,'median_std_ratio',NaN,'all_same_order',false);
end
