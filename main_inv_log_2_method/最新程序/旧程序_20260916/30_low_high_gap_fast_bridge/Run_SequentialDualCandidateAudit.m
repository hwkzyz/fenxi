function Result = Run_SequentialDualCandidateAudit(runMode)
%RUN_SEQUENTIALDUALCANDIDATEAUDIT Audit sequential dual-frequency candidates.
if nargin < 1, runMode = "smoke"; end
root = fileparts(mfilename('fullpath')); mainDir = fileparts(root);
addpath(mainDir,'-begin'); addpath(fullfile(mainDir,'local_func'),'-begin');
ctx = load_inv_log_2_project_context(mainDir); trust = load_fixed_trust_domain(mainDir);
lib = make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),trust);
base = make_inv_log_2_demo_case(ctx,.2);
cases = {[500 1300],[700 1200],[900 1400]};
if runMode == "smoke", gapPairs=[.2 .2;.8 .5]; else, gapPairs=[.2 .2;.5 .5;.8 .8;.5 .2;.8 .5]; end
rows = repmat(row0(),0,1);
for ig=1:size(gapPairs,1)
    gLow=gapPairs(ig,1); gHigh=gapPairs(ig,2); c=base.cfgCase;
    c.RPM_low=min(c.RPM_high,300); c.NumRevs_low=8;
    lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),c,lib.domain,20,930000+ig);
    low=aggregate_low_speed_template(lowData,c.alpha_k,c.R_tip,lib.domain,lib.xGrid);
    st=estimate_highspeed_static_gap_raw(low.mapped,lib,c);
    for ic=1:numel(cases)
        fTrue=cases{ic}; d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gHigh,z),c.RPM_high,c.NumRevs_high,c.fs,c.R_tip,c.alpha_k,0,lib.domain,[.25 .15],fTrue,[pi/4 -pi/3],'noise_ratio');
        sig=rms(d.V_clean-min(d.V_clean)); rng(940000+10*ig+ic,'twister'); d.V_cap=d.V_clean+sig/10*randn(size(d.V_clean));
        m=map_highspeed_to_space(d,c.alpha_k,c.R_tip,lib.domain,.02); m=subset_map(m,600);
        tc=tic; pairs=sequential_candidates(m,lib,st.gHat); elapsed=toc(tc);
        r=row0(); r.case_name=sprintf('%d+%d',fTrue); r.g_low=gLow; r.g_high=gHigh; r.elapsed_s=elapsed;
        err=max(abs(pairs-repmat(fTrue,size(pairs,1),1)),[],2); r.truth_recalled=any(err<=25); r.best_pair_error=min(err); rows(end+1,1)=r; %#ok<AGROW>
    end
end
Audit=struct2table(rows); out=fullfile(root,'output','sequential_audit'); if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'sequential_candidate_audit.csv')); save(fullfile(out,'sequential_candidate_audit.mat'),'Audit','-v7.3'); Result=struct('Audit',Audit,'outputDir',out);
end

function pairs=sequential_candidates(m,lib,g0)
t=m.t_v(:); x=m.x_v(:); V=m.V_a(:); fg=300:25:1500;
gg=linspace(max(.05,g0-.45),min(1.5,g0+.45),31); first=[]; firstScore=[];
for ig=1:numel(gg)
    g=gg(ig); F=eval_gap_template(lib,g,x); Fx=eval_gap_derivative(lib,g,x);
    for ik=1:numel(fg)
        S=sin(2*pi*fg(ik)*t); C=cos(2*pi*fg(ik)*t); B=[Fx,-Fx.*S,-Fx.*C]; q=B\(V-F);
        u=q(2)*S+q(3)*C; s=mean((V-eval_gap_template(lib,g,x-q(1)-u)).^2);
        first(end+1,:)=[g q(1) q(2) q(3) fg(ik)]; firstScore(end+1,1)=s; %#ok<AGROW>
    end
end
[~,ord]=sort(firstScore); ord=ord(1:min(12,numel(ord))); pairs=[]; pairScore=[];
for ii=ord(:)'
    g=first(ii,1); dx=first(ii,2); f1=first(ii,5); S1=sin(2*pi*f1*t); C1=cos(2*pi*f1*t);
    u1=first(ii,3)*S1+first(ii,4)*C1; x1=x-dx-u1; F1=eval_gap_template(lib,g,x1); Fx=eval_gap_derivative(lib,g,x1); R=V-F1;
    for ik=1:numel(fg)
        f2=fg(ik); if abs(f2-f1)<20,continue;end; S=sin(2*pi*f2*t); C=cos(2*pi*f2*t); a=[-Fx.*S,-Fx.*C]\R; u2=a(1)*S+a(2)*C;
        pairScore(end+1,1)=mean((V-eval_gap_template(lib,g,x1-u2)).^2); pairs(end+1,:)=[f1 f2]; %#ok<AGROW>
    end
end
[~,ord]=sort(pairScore); pairs=pairs(ord(1:min(12,numel(ord))),:); pairs=sort(pairs,2); pairs=unique(pairs,'rows','stable');
end

function m=subset_map(m,N)
if numel(m.t_v)<=N,return;end; n=numel(m.t_v); ii=round(linspace(1,n,N)); fn=fieldnames(m);
for k=1:numel(fn), v=m.(fn{k}); if isnumeric(v)&&isvector(v)&&numel(v)==n,m.(fn{k})=v(ii);end;end
end
function r=row0(), r=struct('case_name',"",'g_low',NaN,'g_high',NaN,'elapsed_s',NaN,'truth_recalled',false,'best_pair_error',NaN); end
