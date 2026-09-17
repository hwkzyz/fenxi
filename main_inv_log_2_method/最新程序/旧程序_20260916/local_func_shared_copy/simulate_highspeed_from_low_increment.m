function Data = simulate_highspeed_from_low_increment(pathModel, deltaGap, cfg, noiseLevel, noiseMode, seed)
%SIMULATE_HIGHSPEED_FROM_LOW_INCREMENT  Generate the theoretical high-speed model.
% The low-speed measured template is retained as the baseline and the static
% response surface contributes only the high-low clearance increment.

if nargin < 6 || isempty(seed), seed = 1; end
if nargin < 5 || isempty(noiseMode), noiseMode = 'snr_db'; end
if nargin < 4 || isempty(noiseLevel), noiseLevel = Inf; end
rng(seed, 'twister');

RPM = cfg.RPM_high; nRev = cfg.NumRevs_high; fs = cfg.fs;
Omega = RPM * 2*pi / 60; T = 2*pi / Omega; dt = 1/fs;
t = (0:dt:(nRev+1)*T)';
oprDelay = 0.08*T; T_opr = oprDelay + (0:nRev)'*T;
domain = pathModel.templateLib.domain;
xHalf = 0.52*diff(domain);
A = cfg.A_true; f = cfg.f_true; phi = cfg.phi_true;
u = zeros(size(t));
for q=1:numel(A), u = u + A(q)*sin(2*pi*f(q)*t + phi(q)); end

baselineTemplate = min(pathModel.templateLow(:));
Vclean = baselineTemplate*ones(size(t));
Vstatic = baselineTemplate*ones(size(t));
analysisMask = false(size(t));
for m=1:nRev
    for k=1:numel(cfg.alpha_k)
        Te=T_opr(m)+cfg.alpha_k(k)/Omega;
        idx=abs(t-Te)<=xHalf/(Omega*cfg.R_tip);
        xNom=Omega*cfg.R_tip*(t(idx)-Te); xPhys=xNom-u(idx);
        inside=xPhys>=domain(1)&xPhys<=domain(2); idxAll=find(idx); idxUse=idxAll(inside);
        if isempty(idxUse),continue;end
        Vclean(idxUse)=eval_path_increment_template(pathModel,deltaGap,xPhys(inside),k);
        insideStatic=xNom>=domain(1)&xNom<=domain(2);idxStatic=idxAll(insideStatic);
        Vstatic(idxStatic)=eval_path_increment_template(pathModel,deltaGap,xNom(insideStatic),k);
        analysisMask(idxStatic)=true;
    end
end

signalRange=max(Vclean)-min(Vclean); baseline=min(Vclean);
if isinf(noiseLevel), noiseStd=0; snrDb=Inf;
elseif strcmpi(noiseMode,'snr_db')
    active=abs(Vclean-baseline)>0.02*signalRange; if ~any(active),active=true(size(Vclean));end
    signalRms=rms(Vclean(active)-baseline); noiseStd=signalRms/(10^(noiseLevel/20)); snrDb=noiseLevel;
elseif strcmpi(noiseMode,'fixed_std')
    if ~isscalar(noiseLevel)||~isfinite(noiseLevel)||noiseLevel<0,error('fixed_std noiseLevel must be nonnegative.');end
    noiseStd=noiseLevel;active=abs(Vclean-baseline)>0.02*signalRange;if ~any(active),active=true(size(Vclean));end
    signalRms=rms(Vclean(active)-baseline);snrDb=20*log10(signalRms/max(noiseStd,eps));
elseif strcmpi(noiseMode,'vibration_snr_db')
    vibrationRms=rms(Vclean(analysisMask)-Vstatic(analysisMask));noiseStd=vibrationRms/(10^(noiseLevel/20));snrDb=noiseLevel;
else
    noiseStd=noiseLevel*signalRange; snrDb=20*log10(signalRange/max(noiseStd,eps));
end
Data=struct('t',t,'V_cap',Vclean+noiseStd*randn(size(Vclean)),'V_clean',Vclean, ...
    'T_opr_truth',T_opr,'u_truth',u,'noise_std',noiseStd,'signal_range',signalRange, ...
    'noise_mode',char(noiseMode),'noise_level',noiseLevel,'snr_db_equiv',snrDb,...
    'V_static',Vstatic,'vibration_rms',rms(Vclean(analysisMask)-Vstatic(analysisMask)),...
    'vibration_analysis_mask',analysisMask);
end
