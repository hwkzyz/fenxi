function R = run_single_near_sync_voltage_vp(map,lib,cfg,staticState)
%RUN_SINGLE_NEAR_SYNC_VOLTAGE_VP Blind single near-synchronous solver.
% u(t)=a*sin(EO*Theta(t)+2*pi*df*(t-tc)) +
%      b*cos(EO*Theta(t)+2*pi*df*(t-tc)).
% EO is enumerated, delta-f is profiled on a coarse grid, and the retained
% candidates are refined with the complete low-increment voltage model.

t=map.t_v(:); x=map.x_v(:); V=map.V_a(:); fRot=cfg.RPM_high/60;
if isfield(map,'theta_v') && numel(map.theta_v)==numel(t)
    theta=map.theta_v(:); angleSource="measured_revolution_timing";
else
    theta=2*pi*fRot*t; angleSource="nominal_constant_speed";
end
tc=median(t);
range=get_field(cfg,'route30SingleFrequencyRangeHz',[5 1500]);
eoGrid=ceil(range(1)/fRot):floor(range(2)/fRot);
if isfield(cfg,'nearSyncEOGrid') && ~isempty(cfg.nearSyncEOGrid)
    eoGrid=unique(round(cfg.nearSyncEOGrid(:).'));
end
eoGrid=eoGrid(eoGrid>=1);
dfRange=get_field(cfg,'nearSyncDeltaFreqBoundsHz',[-2 2]);
dfStep=get_field(cfg,'nearSyncDeltaFreqStepHz',.25);
dfGrid=dfRange(1):dfStep:dfRange(2);
if isempty(dfGrid), dfGrid=0; end

gHalf=get_field(cfg,'route30GapHalfWidthMm',.70); gN=get_field(cfg,'route30GapCount',61);
if isfield(lib,'fixedPathIncrement')&&lib.fixedPathIncrement
    gGrid=linspace(staticState.gHat-gHalf,staticState.gHat+gHalf,gN);
else
    lo=max(.05,min(lib.gapTrain)-cfg.rawGapSearchMargin);
    hi=max(lib.gapTrain)+cfg.rawGapSearchMargin;
    gGrid=linspace(max(lo,staticState.gHat-gHalf),min(hi,staticState.gHat+gHalf),gN);
end

bank=repmat(empty_candidate(),0,1); tc0=tic;
for g=gGrid
    F0=eval_gap_template(lib,g,x); Fx=eval_gap_derivative(lib,g,x);
    valid=isfinite(F0)&isfinite(Fx)&isfinite(V)&abs(Fx)>eps;
    if nnz(valid)<8, continue; end
    q=-Fx(valid); y=V(valid)-F0(valid); tv=t(valid); th=theta(valid);
    for eo=eoGrid
        for df=dfGrid
            phase=eo*th+2*pi*df*(tv-tc);
            X=[q,q.*sin(phase),q.*cos(phase)]; beta=X\y; r=y-X*beta;
            bank(end+1)=pack_candidate(g,eo,df,beta,mean(r.^2)); %#ok<AGROW>
        end
    end
end
if isempty(bank), error('invlog2:NoNearSyncCandidate','No near-synchronous candidate was found.'); end
bank=exact_replay_per_frequency(bank,eoGrid,dfGrid,t,theta,tc,x,V,lib,cfg);
[~,ord]=sort([bank.sse],'ascend'); keep=min(get_field(cfg,'nearSyncTopK',12),numel(ord)); ord=ord(1:keep);
refined=repmat(empty_fit(),0,1); refineTic=tic;
for ii=1:numel(ord)
    c=bank(ord(ii)); A=hypot(c.beta(2),c.beta(3)); phi=atan2(c.beta(3),c.beta(2));
    z0=[c.g,c.dx_guess,A,phi,c.df];
    gh=get_field(cfg,'route30SingleGapHalfWidthMm',.12); dh=get_field(cfg,'route30DxHalfWidthMm',.5);
    au=get_field(cfg,'route30AmplitudeUpperMm',1.5); dfh=max(dfStep,get_field(cfg,'nearSyncRefineHalfWidthHz',2));
    lb=[c.g-gh,-dh,0,-pi,max(dfRange(1),c.df-dfh)]; ub=[c.g+gh,dh,au,pi,min(dfRange(2),c.df+dfh)];
    if ~(isfield(lib,'fixedPathIncrement')&&lib.fixedPathIncrement)
        lb(1)=max(lb(1),min(lib.gapTrain)-cfg.rawGapSearchMargin);
        ub(1)=min(ub(1),max(lib.gapTrain)+cfg.rawGapSearchMargin);
    end
    res=@(z)near_sync_residual(z,c.eo,t,theta,tc,V,x,lib);
    % Preserve the exact-replay point as a valid candidate. A bounded
    % nonlinear solver can occasionally move away from a good coarse point
    % along a nearly flat (gap, amplitude, delta-f) valley.
    r0=res(z0); refined(end+1)=pack_fit(z0,c.eo,fRot,t,theta,tc,V,x,lib,struct('solver',"exact_replay"),sqrt(mean(r0.^2))); %#ok<AGROW>
    [z,si]=solve_lsq_bounded(res,z0,lb,ub,get_field(cfg,'nearSyncMaxIter',220),1e-12,1e-12);
    % The nonlinear objective is not guaranteed to be convex. Re-start once
    % from the VP linear estimate and from the zero-detuning center; this is
    % an estimator-only multistart and uses no truth information.
    starts=[z0; z0(1),z0(2),max(0.5*A),phi,0];
    for is=1:size(starts,1)
        [zs,sis]=solve_lsq_bounded(res,starts(is,:),lb,ub,get_field(cfg,'nearSyncMaxIter',220),1e-12,1e-12);
        if norm(res(zs))<norm(res(z)), z=zs; si=sis; end
    end
    r=res(z); refined(end+1)=pack_fit(z,c.eo,fRot,t,theta,tc,V,x,lib,si,sqrt(mean(r.^2))); %#ok<AGROW>
end

function bank=exact_replay_per_frequency(bank,eoGrid,dfGrid,t,theta,tc,x,V,lib,cfg)
% Keep one exact-forward candidate for every (EO,df). This prevents a
% large-amplitude waveform from being discarded by the small-motion VP
% approximation before the complete nonlinear model is evaluated.
Agrid=get_field(cfg,'nearSyncAmplitudeSeedGridMm',[.05 .1 .2 .3 .4]);
Pgrid=get_field(cfg,'nearSyncPhaseSeedGrid',0:pi/4:(2*pi-pi/4));
[AA,PP]=ndgrid(Agrid,Pgrid); avec=AA(:).'; pvec=PP(:).';
for eo=eoGrid
    for df=dfGrid
        q=find([bank.eo]==eo & abs([bank.df]-df)<1e-12);
        if isempty(q),continue;end
        [~,jq]=sort([bank(q).sse],'ascend'); q=q(jq(1:min(2,numel(jq))));
        for iq=q(:).'
            c=bank(iq); phase=eo*theta+2*pi*df*(t-tc);
            U=sin(phase+pvec).*avec;
            Vfit=eval_gap_template(lib,c.g,x-c.dx_guess-U); r=V-Vfit;
            valid=isfinite(r); r(~valid)=0; sse=sum(r.^2,1)./max(sum(valid,1),1);
            [c.sse,k]=min(sse); c.beta(2)=avec(k)*cos(pvec(k)); c.beta(3)=avec(k)*sin(pvec(k)); bank(iq)=c;
        end
    end
end
end
[~,ix]=sort([refined.rmse]); best=refined(ix(1)); sortedRmse=[refined(ix).rmse];
if numel(sortedRmse)>1, margin=(sortedRmse(2)-sortedRmse(1))/max(sortedRmse(1),eps); else, margin=NaN; end
if best.A < get_field(cfg,'nearSyncAmplitudeFloorMm',.02) || ~isfinite(margin) || margin<get_field(cfg,'nearSyncAmbiguityMargin',.01)
    status="ambiguous_near_sync";
else
    status="identified";
end
R=struct('method',"single_near_sync_voltage_vp",'model_order',1,'mode',"single_near_sync",...
    'g_used',best.g,'dx_used',best.dx,'A_id',best.A,'phi_id',best.phi,'f_id',best.f,...
    'eo_id',best.eo,'delta_f_hz',best.df,'rmse',best.rmse,'VFit',best.VFit,'fit',best,...
    'angle_source',angleSource,'eo_grid',eoGrid(:),'delta_f_grid_hz',dfGrid(:),...
    'candidate_count',numel(bank),'refine_count',numel(refined),'second_best_margin',margin,...
    'identification_status',status,'used_sample_count',numel(t),'all_trusted_samples_used',true,...
    'staticState',staticState,'solver_time_s',toc(tc0),'refine_time_s',toc(refineTic));
end

function c=pack_candidate(g,eo,df,beta,sse)
c=struct('g',g,'eo',eo,'df',df,'beta',beta(:).','sse',sse,'dx_guess',beta(1));
end
function c=empty_candidate()
c=struct('g',NaN,'eo',NaN,'df',NaN,'beta',nan(1,3),'sse',Inf,'dx_guess',NaN);
end

function f=pack_fit(z,eo,fRot,t,theta,tc,V,x,lib,si,rmse)
f=empty_fit(); f.g=z(1);f.dx=z(2);f.A=z(3);f.phi=z(4);f.df=z(5);f.eo=eo;
f.f=eo*fRot+z(5);
u=z(3)*sin(eo*theta+2*pi*z(5)*(t-tc)+z(4)); f.VFit=eval_gap_template(lib,z(1),x-z(2)-u);
f.rmse=rmse; f.solve_info=si;
end

function r=near_sync_residual(z,eo,t,theta,tc,V,x,lib)
u=z(3)*sin(eo*theta+2*pi*z(5)*(t-tc)+z(4));
r=V-eval_gap_template(lib,z(1),x-z(2)-u); r(~isfinite(r))=10*max(std(V),1e-3);
end

function f=empty_fit()
f=struct('g',NaN,'dx',NaN,'A',NaN,'phi',NaN,'df',NaN,'eo',NaN,'f',NaN,'VFit',[],'rmse',Inf,'solve_info',struct());
end
function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
