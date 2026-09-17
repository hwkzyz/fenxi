function Result=Run_01_EnergyMatchedPhaseAperture()
% Same samples and fixed waveform sensitivities; only phase aperture changes.
root=fileparts(mfilename('fullpath'));programDir=fileparts(root);
addpath(programDir,'-begin');addpath(fullfile(programDir,'local_func'),'-begin');
src=fullfile(programDir,'identifiability_mechanism_analysis','output',...
    '02_exact_frozen_conditional_geometry','exact_frozen_conditional_geometry.mat');
if ~isfile(src),error('route32:MissingRoute31','Run Route 31 geometry first.');end
S=load(src,'map','model','cfg','cases');map=S.map;model=S.model;cfg=S.cfg;cases=S.cases;
gammaGrid=[0 .125 .25 .5 .75 1];rows=repmat(row0(),0,1);
theta=map.theta_v(:);x=map.x_v(:);sid=map.S_v(:);thetaCenter=theta-x/cfg.R_tip;
for ic=1:numel(cases)
    c=cases(ic);
    % Evaluate response derivatives once at the exact physical trajectory.
    u=c.A(1)*sin(c.eo(1)*theta+c.phi(1))+c.A(2)*sin(c.eo(2)*theta+c.phi(2));
    z=x-u;fx=eval_gap_derivative(model,.6,z,sid);h=.002;
    fg=(eval_gap_template(model,.6+h,z,sid)-eval_gap_template(model,.6-h,z,sid))/(2*h);
    valid=isfinite(fx)&isfinite(fg);fx=fx(valid);fg=fg(valid);
    th=theta(valid);thc=thetaCenter(valid);
    local=repmat(row0(),numel(gammaGrid),1);Jstore=cell(numel(gammaGrid),1);
    for ig=1:numel(gammaGrid)
        g=gammaGrid(ig);phase=thc+g*(th-thc);
        Jv=-fx.*[sin(c.eo(1)*phase),cos(c.eo(1)*phase),...
            sin(c.eo(2)*phase),cos(c.eo(2)*phase)];
        Jgap=fg;Jdx=-fx;
        [Jjoint,rankNJoint]=condition_columns([Jgap Jv],Jdx);
        [Jvib,rankNVib]=condition_columns(Jv,[Jgap Jdx]);
        mj=raw_metrics(Jjoint);mv=raw_metrics(Jvib);
        rr=row0();rr.case_id=c.id;rr.gamma=g;rr.sample_count=numel(fx);
        rr.fixed_sensitivity_energy=sum(fx.^2);rr.joint_raw_trace=mj.trace;
        rr.joint_rank=mj.rank;rr.joint_lambda_min_raw=mj.lambdaMin;
        rr.joint_lambda_fraction=mj.lambdaMin/max(mj.trace,eps);
        rr.joint_condition=mj.condition;rr.vibration_rank=mv.rank;
        rr.vibration_lambda_fraction=mv.lambdaMin/max(mv.trace,eps);
        rr.vibration_condition=mv.condition;rr.rank_nuisance_joint=rankNJoint;
        rr.rank_nuisance_vibration=rankNVib;local(ig)=rr;Jstore{ig}=Jjoint;
    end
    traceRef=local(end).joint_raw_trace;
    for ig=1:numel(local)
        m=raw_metrics(Jstore{ig}*sqrt(traceRef/max(local(ig).joint_raw_trace,eps)));
        local(ig).joint_trace_matched=traceRef;
        local(ig).joint_lambda_min_trace_matched=m.lambdaMin;
        rows(end+1,1)=local(ig); %#ok<AGROW>
    end
end
Scan=struct2table(rows);Gate=make_gate(Scan);
out=fullfile(root,'output','01_energy_matched_phase_aperture');if ~exist(out,'dir'),mkdir(out);end
writetable(Scan,fullfile(out,'phase_aperture_scan.csv'));writetable(Gate,fullfile(out,'gate.csv'));
save(fullfile(out,'phase_aperture.mat'),'Scan','Gate','gammaGrid','-v7.3');write_report(out,Gate);
disp(Gate);Result=struct('Scan',Scan,'Gate',Gate,'outputDir',out);
end

function [Jc,rn]=condition_columns(Jp,Jn)
[Q,rn]=orthb(Jn);Jc=Jp-Q*(Q.'*Jp);
end
function M=raw_metrics(J)
F=(J.'*J);e=sort(max(0,real(eig((F+F.')/2))),'ascend');
tol=1e-10*max(norm(J,2),eps);r=rank(J,tol);
M=struct('trace',trace(F),'lambdaMin',e(1),'rank',r,...
    'condition',e(end)/max(e(1),eps));
end
function [Q,r]=orthb(A)
if isempty(A),Q=zeros(size(A,1),0);r=0;return;end
[U,S,~]=svd(A,'econ');s=diag(S);r=sum(s>1e-10*max(s));Q=U(:,1:r);
end
function G=make_gate(T)
ids=unique(T.case_id,'stable');rows=repmat(struct('case_id',"",...
    'rank_gamma0',NaN,'rank_gamma1',NaN,'lambda_gain_gamma1_vs_quarter',NaN,...
    'spearman_gamma_matched_information',NaN,'gate_pass',false),numel(ids),1);
for i=1:numel(ids)
    q=sortrows(T(T.case_id==ids(i),:),'gamma');a=q(q.gamma==0,:);b=q(q.gamma==1,:);
    k=q(q.gamma==.25,:);rho=corr(q.gamma,q.joint_lambda_min_trace_matched,'Type','Spearman');
    rows(i)=struct('case_id',ids(i),'rank_gamma0',a.joint_rank,'rank_gamma1',b.joint_rank,...
        'lambda_gain_gamma1_vs_quarter',b.joint_lambda_min_trace_matched/max(k.joint_lambda_min_trace_matched,eps),...
        'spearman_gamma_matched_information',rho,...
        'gate_pass',b.joint_rank==5&&b.joint_lambda_min_trace_matched>1.25*k.joint_lambda_min_trace_matched&&rho>.8);
end
G=struct2table(rows);
end
function write_report(out,G)
fid=fopen(fullfile(out,'README.md'),'w');fprintf(fid,'# Energy-matched phase-aperture result\n\n');
fprintf(fid,'All rows use identical samples and frozen exact-trajectory response derivatives. Only intra-passage phase aperture changes; the conditional Jacobian trace is then matched to the exact case.\n\n');
for i=1:height(G),fprintf(fid,'- %s: rank %g -> %g; exact/quarter-aperture weakest-information gain %.6g; Spearman %.6g; gate %d.\n',...
    G.case_id(i),G.rank_gamma0(i),G.rank_gamma1(i),G.lambda_gain_gamma1_vs_quarter(i),...
    G.spearman_gamma_matched_information(i),G.gate_pass(i));end
fprintf(fid,'\nThis is a local phase-only intervention. It supports an information-direction mechanism, not global EO uniqueness.\n');fclose(fid);
end
function r=row0()
r=struct('case_id',"",'gamma',NaN,'sample_count',NaN,'fixed_sensitivity_energy',NaN,...
    'joint_raw_trace',NaN,'joint_trace_matched',NaN,'joint_rank',NaN,...
    'joint_lambda_min_raw',NaN,'joint_lambda_min_trace_matched',NaN,...
    'joint_lambda_fraction',NaN,'joint_condition',NaN,'vibration_rank',NaN,...
    'vibration_lambda_fraction',NaN,'vibration_condition',NaN,...
    'rank_nuisance_joint',NaN,'rank_nuisance_vibration',NaN);
end
