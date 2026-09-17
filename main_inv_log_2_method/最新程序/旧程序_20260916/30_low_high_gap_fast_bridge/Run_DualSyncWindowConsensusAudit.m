function Result=Run_DualSyncWindowConsensusAudit(seedOffset)
%RUN_DUALSYNCWINDOWCONSENSUSAUDIT Test independent-window voting on a weak pair.
if nargin<1,seedOffset=2000;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;cfg.NumRevs_high=32;cfg.snrDb=5;
fRot=cfg.RPM_high/60;eoTrue=[6 7];fTrue=eoTrue*fRot;A=[.1 .1];phi=[.13 2.21];
low=simulate_low_speed_template(@(z)eval_gap_template(lib,.2,z),cfg,lib.domain,5,220700+seedOffset);
lowMap=map_highspeed_to_space(low,cfg.alpha_k,cfg.R_tip,lib.domain,0);
state=estimate_highspeed_static_gap_raw(lowMap,lib,cfg);
d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,.8,z),...
    cfg.RPM_high,32,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,A,fTrue,phi,'noise_ratio');
sig=rms(d.V_clean-min(d.V_clean));rng(220750+seedOffset,'twister');
d.V_cap=d.V_clean+sig/10^(5/20)*randn(size(d.V_clean));
map=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
rev=unique(map.rev_v(:)).';windowSize=8;nWindow=floor(numel(rev)/windowSize);
rows=repmat(row0(),nWindow+1,1);
ticFit=tic;whole=run_dual_sync_voltage_vp(map,lib,cfg,state);
rows(1)=pack("whole",whole,toc(ticFit),eoTrue);
for iw=1:nWindow
    useRev=rev((iw-1)*windowSize+(1:windowSize));sub=subset_map(map,ismember(map.rev_v,useRev));
    ticFit=tic;q=run_dual_sync_voltage_vp(sub,lib,cfg,state);
    rows(iw+1)=pack("window_"+iw,q,toc(ticFit),eoTrue);
end
Audit=struct2table(rows);windowRows=Audit(2:end,:);
pairs=windowRows.eo_est;[u,~,ic]=unique(pairs);count=accumarray(ic,1);[votes,k]=max(count);
consensusPair=u(k);consensusSuccess=consensusPair==string(mat2str(eoTrue));
out=fullfile(root,'output',sprintf('dual_sync_window_consensus_seed_%d',seedOffset));
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'window_consensus.csv'));
save(fullfile(out,'window_consensus.mat'),'Audit','consensusPair','votes','consensusSuccess');disp(Audit);
Result=struct('Audit',Audit,'consensusPair',consensusPair,'votes',votes,...
    'consensusSuccess',consensusSuccess,'outputDir',out);
end

function r=pack(label,q,elapsed,eoTrue)
r=struct('window',label,'eo_est',string(mat2str(q.eo_id)),...
    'amplitude_est',string(mat2str(q.A_id,7)),'amplitude_est_min',min(q.A_id),...
    'rmse',q.rmse,'noise_normalized_margin',q.noise_normalized_margin,...
    'confident',q.identification_confident,'elapsed_s',elapsed,...
    'success',isequal(q.eo_id,eoTrue));
end

function out=subset_map(in,mask)
out=in;n=numel(mask);names=fieldnames(in);
for i=1:numel(names)
    v=in.(names{i});
    if isnumeric(v)&&numel(v)==n
        out.(names{i})=v(mask);
    end
end
end

function r=row0()
r=struct('window',"",'eo_est',"",'amplitude_est',"",'amplitude_est_min',NaN,...
    'rmse',NaN,'noise_normalized_margin',NaN,...
    'confident',false,'elapsed_s',NaN,'success',false);
end
