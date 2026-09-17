function R = run_single_frequency_joint(highMap,templateLib,cfg,staticState,mode)
%RUN_SINGLE_FREQUENCY_JOINT  Inverse-frequency plus voltage-domain S/A solver.
% The inverse map supplies frequency candidates only. Final parameters are
% always selected by the complete forward voltage residual.
if nargin<5||isempty(mode),mode="async";end
gapState=estimate_projected_gap_only(highMap,templateLib,cfg,staticState);
if strcmpi(string(mode),'sync') && isfield(cfg,'singleSyncFrequencyHz')
    fGrid=cfg.singleSyncFrequencyHz;
else
    fGrid=cfg.singleFreqGrid;
end
[inverseSeeds,inverseMap]=single_inverse_seed(highMap,templateLib,gapState.g_used,fGrid,cfg);
inverse0=inverseSeeds(1);
if inverse0.valid_fraction<=0,error('invlog2:InverseSeedFailed','No valid inverse-map samples.');end
best=[]; bestRmse=inf;
for is=1:numel(inverseSeeds)
    seed=inverseSeeds(is);
    if strcmpi(string(mode),'sync'),freqHalf=0;else,freqHalf=get_cfg(cfg,'singleAsyncRefineHalfWidthHz',10);end
    direct=fit_single_frequency_ab_multistart(highMap,templateLib,cfg,gapState,seed.f, ...
        struct('amplitudeInit',max(seed.A,get_cfg(cfg,'nonlinearAmpSeedMm',0.18)), ...
        'dxInit',seed.dx,'freqHalfWidth',freqHalf, ...
        'frequencyStartStepHz',get_cfg(cfg,'singleAsyncFineStepHz',0.1), ...
        'inverseMap',inverseMap));
    if direct.rmse<bestRmse,bestRmse=direct.rmse;best=direct;end
end
ident=compute_dx_vibration_identifiability(highMap,templateLib,best.g,best.dx,best.A,best.phi,best.f);
R=struct('method',"inverse_frequency_direct_ab",'mode',string(mode),'g_used',best.g,'dx_used',best.dx, ...
 'f_id',best.f,'A_id',best.A,'phi_id',best.phi,'fit',best,'rmse',bestRmse, ...
 'staticState',staticState,'projectedGapState',gapState, ...
 'inverseSeedFit',inverse0,'inverseCandidateFits',inverseSeeds, ...
 'inverseMap',inverseMap,'identifiability',ident,'vpUsed',false,'directABUsed',best.seed_source=="direct_ab_screen");
end

function [fits,inv]=single_inverse_seed(highMap,lib,g,fGrid,cfg)
inv=inverse_map_local_displacement(highMap,lib,g,0,cfg);
t=highMap.t_v(:); ok=inv.valid(:)&isfinite(inv.u_inv(:))&inv.weight(:)>0;
emptyFit=struct('dx',0,'A',0,'phi',0,'f',fGrid(1),'p',[0,0,fGrid(1)], ...
    'valid_fraction',mean(ok),'inverse_voltage_rmse',inv.rms_voltage_residual, ...
    'linear_sse',inf,'candidate_rho',NaN,'adaptive_expanded',false);
if nnz(ok)<6,fits=emptyFit;return;end
tt=t(ok); uu=inv.u_inv(ok); ww=inv.weight(ok); sw=sqrt(ww);
target=get_cfg(cfg,'singleInverseScanTargetSamples',5000);
ii=pick_evenly_spaced_indices(numel(tt),target);tt=tt(ii);uu=uu(ii);sw=sw(ii);
score=inf(numel(fGrid),1);coef=zeros(numel(fGrid),3);
for k=1:numel(fGrid),f=fGrid(k);
    X=[ones(size(tt)),sin(2*pi*f*tt),cos(2*pi*f*tt)]; b=(X.*sw)\(uu.*sw);
    r=(uu-X*b).*sw;score(k)=dot(r,r);coef(k,:)=b(:).';
end
[~,ord]=sort(score);baseK=get_cfg(cfg,'singleAsyncTopK',3);maxK=max(baseK,get_cfg(cfg,'singleAsyncMaxTopK',5));
snr=NaN;if isfield(highMap,'snr_db_equiv'),snr=highMap.snr_db_equiv;end
if ~isfinite(snr)&&isfield(cfg,'snrDb')&&isscalar(cfg.snrDb),snr=cfg.snrDb;end
if isfinite(snr)
    if snr<get_cfg(cfg,'singleAsyncLowSNRThreshold',20)
        baseK=get_cfg(cfg,'singleAsyncLowSNRTopK',5);maxK=max(maxK,baseK);
    end
    if snr<get_cfg(cfg,'singleAsyncVeryLowSNRThreshold',10)
        baseK=get_cfg(cfg,'singleAsyncVeryLowSNRTopK',8);maxK=max(maxK,baseK);
    end
end
minSep=get_cfg(cfg,'singleAsyncCandidateMinSeparationHz',20);keep=zeros(1,min(maxK,numel(ord)));nk=0;
for i=1:numel(ord)
    if nk==0||all(abs(fGrid(ord(i))-fGrid(keep(1:nk)))>=minSep)
        nk=nk+1;keep(nk)=ord(i);if nk==numel(keep),break;end
    end
end
keep=keep(1:nk);expanded=false;rho=inf;
if nk>baseK
    rho=(score(keep(baseK+1))-score(keep(1)))/max(abs(score(keep(1))),eps);
    expanded=rho<get_cfg(cfg,'singleAsyncExpandRhoThreshold',0.05);
end
if ~expanded,keep=keep(1:min(baseK,numel(keep)));end
fits=repmat(emptyFit,numel(keep),1);
for i=1:numel(keep)
    k=keep(i);b=coef(k,:);fits(i).dx=b(1);fits(i).A=hypot(b(2),b(3));
    fits(i).phi=atan2(b(3),b(2));fits(i).f=fGrid(k);fits(i).p=[fits(i).A,fits(i).phi,fits(i).f];
    fits(i).linear_sse=score(k);fits(i).candidate_rho=rho;fits(i).adaptive_expanded=expanded;
end
end

function v=get_cfg(s,n,d),if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end,end
