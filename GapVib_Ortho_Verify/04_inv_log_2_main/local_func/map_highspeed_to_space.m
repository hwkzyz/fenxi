function highMap = map_highspeed_to_space(Data_High, alpha_k, R_tip, domain, activeLevel)
%map_highspeed_to_space  Convert time-domain passages to fixed x-space windows.
%
% The voltage threshold is retained only as diagnostic metadata. The fitted
% waveform samples are selected by the spatial window so that noise does not
% decide which parts of the waveform enter the gap-vibration model.

T_opr = Data_High.T_opr_truth(:);
t = Data_High.t(:);
V = Data_High.V_cap(:);
xHalf = 0.48 * diff(domain);

base = local_percentile(V, 5);
span = local_percentile(V, 99) - base;
activeThreshold = base + activeLevel * span;

tAll = [];
vAll = [];
xAll = [];
revAll = [];
sAll = [];
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
    end
end

highMap = struct('t_v', tAll, 'V_a', vAll, 'x_v', xAll, ...
    'rev_v', revAll, 'S_v', sAll, 'activeThreshold', activeThreshold, ...
    'selectionMode', "fixed_spatial_window");
end
