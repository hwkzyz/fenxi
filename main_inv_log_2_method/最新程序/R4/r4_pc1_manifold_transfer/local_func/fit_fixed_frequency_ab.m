function R = fit_fixed_frequency_ab(highMap, templateLib, cfg, staticState, freqPair, opts)
%FIT_FIXED_FREQUENCY_AB  Fit gap/offset and sine-cosine coefficients at fixed frequencies.
% This is a diagnostic structured-frequency route. Frequencies are held
% exactly at the supplied pair, while g and dx remain free within bounds.
if nargin < 6, opts = struct(); end
t = highMap.t_v(:); x = highMap.x_v(:); V = highMap.V_a(:);
f = sort(freqPair(:).');
ampUpper = get_opt(opts,'ampUpperBound',1.5);
gapHalf = get_opt(opts,'gapHalfWidth',0.12);
dxHalf = get_opt(opts,'dxHalfWidth',0.20);
maxIter = get_opt(opts,'maxIter',250);
invalidPenalty = 10*max(std(V),1e-3);
amp0 = get_opt(opts,'amplitudeInit',[0.20,0.12]);
phase0 = get_opt(opts,'phaseInit',[0,0]);
lb = [max(0.05,staticState.gHat-gapHalf), staticState.dx0-dxHalf, -ampUpper*ones(1,4)];
ub = [staticState.gHat+gapHalf, staticState.dx0+dxHalf, ampUpper*ones(1,4)];
lb(1)=max(lb(1),min(templateLib.gapTrain)-cfg.rawGapSearchMargin);
ub(1)=min(ub(1),max(templateLib.gapTrain)+cfg.rawGapSearchMargin);
resFun = @(th) fixed_residual(th,t,x,V,templateLib,f,invalidPenalty);
gSeeds = get_opt(opts,'gapSeeds',staticState.gHat + [-0.08,-0.03,0,0.03,0.08]);
dxSeeds = get_opt(opts,'dxSeeds',staticState.dx0 + [-dxHalf, -0.5*dxHalf, 0, 0.5*dxHalf, dxHalf]);
bestSse = inf; theta = []; info = struct();
for ig = 1:numel(gSeeds)
    for id = 1:numel(dxSeeds)
        theta0 = [gSeeds(ig), dxSeeds(id), amp0(1)*cos(phase0(1)), amp0(1)*sin(phase0(1)), ...
            amp0(2)*cos(phase0(2)), amp0(2)*sin(phase0(2))];
        theta0 = min(max(theta0,lb),ub);
        [thetaTry,infoTry] = solve_lsq_bounded(resFun,theta0,lb,ub,maxIter);
        rTry = resFun(thetaTry); sseTry = dot(rTry,rTry);
        if sseTry < bestSse, bestSse = sseTry; theta = thetaTry; info = infoTry; end
    end
end
p = [hypot(theta(3),theta(4)),atan2(theta(4),theta(3)),f(1), ...
     hypot(theta(5),theta(6)),atan2(theta(6),theta(5)),f(2)];
VFit = eval_gap_template(templateLib,theta(1),x-theta(2)-fit_u(p,t));
valid = isfinite(VFit); rmse = sqrt(mean((V(valid)-VFit(valid)).^2));
R = struct('g_used',theta(1),'dx_used',theta(2),'f_id',f,'A_id',p([1,4]), ...
    'phi_id',p([2,5]),'p_id',p,'theta',theta,'VFit',VFit,'rmse',rmse, ...
    'invalid_count',sum(~valid),'solve_info',info,'fixed_frequency',true);
end

function r = fixed_residual(th,t,x,V,lib,f,penalty)
p = [hypot(th(3),th(4)),atan2(th(4),th(3)),f(1), ...
     hypot(th(5),th(6)),atan2(th(6),th(5)),f(2)];
y = eval_gap_template(lib,th(1),x-th(2)-fit_u(p,t)); r = V-y;
bad = ~isfinite(r); if any(bad), r(bad)=penalty; end
end

function v = get_opt(S,n,d)
if isstruct(S)&&isfield(S,n)&&~isempty(S.(n)), v=S.(n); else, v=d; end
end
