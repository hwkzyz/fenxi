function R=run_single_async_voltage_vp(map,lib,cfg,staticState)
%RUN_SINGLE_ASYNC_VOLTAGE_VP Known-order asynchronous single solver.
range=get_field(cfg,'route30SingleFrequencyRangeHz',[5 1500]);
step=get_field(cfg,'route30FrequencyStepHz',5);
cfg.route30SingleGrid=range(1):step:range(2);
cfg.route30GapHalfWidthMm=get_field(cfg,'route30GapHalfWidthMm',.70);
cfg.route30GapCount=get_field(cfg,'route30GapCount',61);
cfg.mainVpGapHalfWidth=cfg.route30GapHalfWidthMm;
cfg.mainVpGapN=cfg.route30GapCount;
seed=fit_single_vp_main(map,lib,cfg,staticState);
candidates=seed.candidates;
if get_field(cfg,'asyncSingleUseDirectReplay',true)
    direct=direct_async_candidates(map,lib,cfg,staticState);
    candidates=[candidates(:);direct(:)];
    [~,candidateOrder]=sort([candidates.J1],'ascend');candidates=candidates(candidateOrder);
end
shortCount=min(get_field(cfg,'asyncSingleShortRefineCount',50),numel(candidates));
cfgShort=cfg;cfgShort.route30SingleMaxIter=get_field(cfg,'asyncSingleShortMaxIter',50);
best=[];bestRmse=Inf;refinedFrequency=[];refinedRmse=[];
for i=1:shortCount
    q=refine_single_async_voltage(map,lib,cfgShort,candidates(i));
    refinedFrequency(end+1,1)=q.f;refinedRmse(end+1,1)=q.rmse; %#ok<AGROW>
    if q.rmse<bestRmse,best=q;bestRmse=q.rmse;end
end
bestSeed=struct('g_used',best.g,'dx_used',best.dx,'A_id',best.A,...
    'phi_id',best.phi,'f_id',best.f);
fit=refine_single_async_voltage(map,lib,cfg,bestSeed);
R=struct('method',"single_async_voltage_vp",'model_order',1,'mode',"single_async",...
    'g_used',fit.g,'dx_used',fit.dx,'A_id',fit.A,'phi_id',fit.phi,...
    'f_id',fit.f,'rmse',fit.rmse,'VFit',fit.VFit,'fit',fit,'seed',seed,...
    'short_refine_count',shortCount,...
    'used_sample_count',numel(map.t_v),'all_trusted_samples_used',true,...
    'staticState',staticState,'frequency_grid_hz',cfg.route30SingleGrid(:));
C=frequency_solution_confidence(fit.f,fit.rmse,refinedFrequency,refinedRmse,numel(map.t_v),cfg);
R=merge_struct(R,C);
end

function a=merge_struct(a,b)
names=fieldnames(b);for i=1:numel(names),a.(names{i})=b.(names{i});end
end

function candidates=direct_async_candidates(map,lib,cfg,staticState)
% Coarse exact-forward replay protects large-amplitude cases from a purely
% first-order VP gate. Frequency remains continuous in the final refinement.
t=map.t_v(:);x=map.x_v(:);V=map.V_a(:);
range=get_field(cfg,'route30SingleFrequencyRangeHz',[5 1500]);
fStep=get_field(cfg,'asyncSingleDirectFrequencyStepHz',10);
fGrid=range(1):fStep:range(2);
gHalf=get_field(cfg,'route30GapHalfWidthMm',.70);gN=get_field(cfg,'asyncSingleDirectGapCount',21);
if isfield(lib,'fixedPathIncrement')&&lib.fixedPathIncrement
    gGrid=linspace(staticState.gHat-gHalf,staticState.gHat+gHalf,gN);
else
    lo=max(.05,min(lib.gapTrain)-cfg.rawGapSearchMargin);hi=max(lib.gapTrain)+cfg.rawGapSearchMargin;
    gGrid=linspace(max(lo,staticState.gHat-gHalf),min(hi,staticState.gHat+gHalf),gN);
end
Agrid=get_field(cfg,'asyncSingleAmplitudeSeedGridMm',[.1 .25 .4]);
Pgrid=get_field(cfg,'asyncSinglePhaseSeedGrid',0:pi/4:(2*pi-pi/4));
[AA,PP]=ndgrid(Agrid,Pgrid);avec=AA(:).';pvec=PP(:).';dx0=staticState.dx0;
rows=repmat(struct('g_used',NaN,'dx_used',NaN,'f_id',NaN,'A_id',NaN,...
    'phi_id',NaN,'p',[],'J1',Inf,'linear_sse',Inf),numel(gGrid)*numel(fGrid),1);
at=0;
for g=gGrid
    for f=fGrid
        at=at+1;U=sin(2*pi*f*t+pvec).*avec;
        VFit=eval_gap_template(lib,g,x-dx0-U);r=V-VFit;
        valid=isfinite(r);r(~valid)=0;sse=sum(r.^2,1)./max(sum(valid,1),1);
        [best,k]=min(sse);A=avec(k);phi=pvec(k);
        rows(at)=struct('g_used',g,'dx_used',dx0,'f_id',f,'A_id',A,...
            'phi_id',phi,'p',[A,phi,f],'J1',best,'linear_sse',best);
    end
end
[~,ord]=sort([rows.J1],'ascend');keep=min(get_field(cfg,'asyncSingleDirectCandidateCount',30),numel(ord));
candidates=rows(ord(1:keep));
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
