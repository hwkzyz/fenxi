function Data = simulate_rotating_waveform_from_template(Fx, RPM, NumRevs, fs, R_tip, alpha_k, noiseLevel, domain, A, f, phi, noiseMode)
%simulate_rotating_waveform_from_template  Generate high-speed waveform data.

if nargin < 12 || isempty(noiseMode)
    noiseMode = 'noise_ratio';
end

Omega = RPM * 2*pi / 60;
T = 2*pi / Omega;
dt = 1 / fs;
tEnd = (NumRevs + 1) * T;
t = (0:dt:tEnd)';

oprDelay = 0.08 * T;
T_opr = oprDelay + (0:NumRevs)' * T;
xHalf = 0.52 * diff(domain);

probe = linspace(domain(1), domain(2), 400)';
yProbe = Fx(probe);
baseline = min(yProbe);
V_clean = baseline * ones(size(t));

u_t = zeros(size(t));
for im = 1:numel(A)
    u_t = u_t + A(im) * sin(2*pi*f(im)*t + phi(im));
end

for m = 1:NumRevs
    for k = 1:numel(alpha_k)
        Te = T_opr(m) + alpha_k(k) / Omega;
        idx = abs(t - Te) <= xHalf / (Omega * R_tip);
        xNom = Omega * R_tip * (t(idx) - Te);
        xPhys = xNom - u_t(idx);
        inDomain = xPhys >= domain(1) & xPhys <= domain(2);
        idxAll = find(idx);
        idxUse = idxAll(inDomain);
        if isempty(idxUse)
            continue;
        end
        V_clean(idxUse) = Fx(xPhys(inDomain));
    end
end

signalRange = max(V_clean) - min(V_clean);
switch lower(string(noiseMode))
    case "noise_ratio"
        noiseStd = noiseLevel * signalRange;
        snrDb = 20 * log10(signalRange / max(noiseStd, eps));
    case "snr_db"
        activeMask = abs(V_clean - baseline) > 0.02 * signalRange;
        if ~any(activeMask)
            activeMask = true(size(V_clean));
        end
        signalRms = rms(V_clean(activeMask) - baseline);
        noiseStd = signalRms / (10^(noiseLevel / 20));
        snrDb = noiseLevel;
    case "fixed_std"
        if ~isscalar(noiseLevel) || ~isfinite(noiseLevel) || noiseLevel < 0
            error('fixed_std noiseLevel must be a finite nonnegative voltage standard deviation.');
        end
        noiseStd = noiseLevel;
        activeMask = abs(V_clean - baseline) > 0.02 * signalRange;
        if ~any(activeMask), activeMask = true(size(V_clean)); end
        signalRms = rms(V_clean(activeMask) - baseline);
        snrDb = 20 * log10(signalRms / max(noiseStd, eps));
    otherwise
        error('Unknown noiseMode: %s', noiseMode);
end

Data = struct();
Data.t = t;
Data.V_cap = V_clean + noiseStd * randn(size(V_clean));
Data.V_clean = V_clean;
Data.T_opr_truth = T_opr;
Data.u_truth = u_t;
Data.noise_std = noiseStd;
Data.signal_range = signalRange;
Data.noise_mode = char(noiseMode);
Data.noise_level = noiseLevel;
Data.snr_db_equiv = snrDb;
% Acquisition geometry metadata. Template builders may use this declared
% support instead of duplicating the generator's window constant.
Data.spatial_window_ratio = 0.52;
Data.spatial_window_half_width_mm = xHalf;
end
