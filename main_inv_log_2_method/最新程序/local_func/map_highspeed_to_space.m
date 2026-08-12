function highMap = map_highspeed_to_space(Data_High, alpha_k, R_tip, domain, activeLevel, windowRatio)
%map_highspeed_to_space  Convert time-domain passages to fixed x-space windows.
%
% The voltage threshold is retained only as diagnostic metadata. The fitted
% waveform samples are selected by the spatial window so that noise does not
% decide which parts of the waveform enter the gap-vibration model.

if nargin<6||isempty(windowRatio),windowRatio=.48;end
validateattributes(windowRatio,{'numeric'},{'scalar','finite','positive'});
T_opr = Data_High.T_opr_truth(:);
t = Data_High.t(:);
V = Data_High.V_cap(:);
xHalf = windowRatio * diff(domain);

base = local_percentile(V, 5);
span = local_percentile(V, 99) - base;
activeThreshold = base + activeLevel * span;

tAll = [];
vAll = [];
xAll = [];
revAll = [];
sAll = [];
thetaAll = [];
for m = 1:numel(T_opr)-1
    Omega = 2*pi / (T_opr(m+1) - T_opr(m));
    for k = 1:numel(alpha_k)
        Te = T_opr(m) + alpha_k(k) / Omega;
        idx = find(abs(t - Te) <= xHalf / (Omega * R_tip));
        if isempty(idx)
            continue;
        end
        xLocal = Omega * R_tip * (t(idx) - Te);
        keep = xLocal >= domain(1) & xLocal <= domain(2);
        idx = idx(keep);
        xLocal = xLocal(keep);
        if isempty(idx)
            continue;
        end
        tAll = [tAll; t(idx)]; %#ok<AGROW>
        vAll = [vAll; V(idx)]; %#ok<AGROW>
        xAll = [xAll; xLocal]; %#ok<AGROW>
        revAll = [revAll; m * ones(numel(idx), 1)]; %#ok<AGROW>
        sAll = [sAll; k * ones(numel(idx), 1)]; %#ok<AGROW>
        thetaLocal = 2*pi*(m-1) + Omega*(t(idx)-T_opr(m));
        thetaAll = [thetaAll; thetaLocal]; %#ok<AGROW>
    end
end

highMap = struct('t_v', tAll, 'V_a', vAll, 'x_v', xAll, ...
    'rev_v', revAll, 'S_v', sAll, 'theta_v', thetaAll, ...
    'activeThreshold', activeThreshold, ...
    'selectionMode', "fixed_spatial_window",'window_ratio',windowRatio);
if isfield(Data_High,'snr_db_equiv')
    highMap.snr_db_equiv=Data_High.snr_db_equiv;
elseif isfield(Data_High,'snrDb')
    highMap.snr_db_equiv=Data_High.snrDb;
else
    highMap.snr_db_equiv=NaN;
end
end
