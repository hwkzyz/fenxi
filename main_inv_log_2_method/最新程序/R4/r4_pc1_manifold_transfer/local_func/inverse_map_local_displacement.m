function inv = inverse_map_local_displacement(highMap, templateLib, gQuery, dxNom, cfg)
%INVERSE_MAP_LOCAL_DISPLACEMENT  Local branch inverse of a fixed-path waveform.
%
% The inverse is used only to generate vibration seeds.  The final estimate
% must still be obtained from the forward voltage-domain objective.  Each
% sample is initialized at its nominal spatial coordinate and iterated on
% the local response branch, which avoids the left/right ambiguity of a
% complete pulse inversion.

if nargin < 5 || isempty(cfg), cfg = struct(); end
xNom = highMap.x_v(:) - dxNom;
vObs = highMap.V_a(:);
maxIter = get_cfg(cfg, 'inverseMapMaxIter', 8);
halfWidth = get_cfg(cfg, 'inverseMapHalfWidthMm', 0.30);
derivFloor = get_cfg(cfg, 'inverseMapDerivativeFloor', 1e-5);
voltageTol = get_cfg(cfg, 'inverseMapVoltageTolerance', inf);
if isfinite(voltageTol)
    solveTol = max(voltageTol, 1e-10);
else
    % Inf disables the post-fit voltage rejection; it must not make the
    % bisection accept its first midpoint unconditionally.
    solveTol = 1e-10;
end

% The fitted path parameter tau is not generally the voltage-peak location
% after the gap increment is applied.  Branch selection must use the actual
% response extremum at the queried gap; otherwise valid samples are inverted
% on the opposite side of the pulse and the recovered displacement collapses.
if isfield(templateLib,'xTemplate') && ~isempty(templateLib.xTemplate)
    xProbe=templateLib.xTemplate(:);
elseif isfield(templateLib,'xGrid') && ~isempty(templateLib.xGrid)
    xProbe=templateLib.xGrid(:);
else
    xProbe=linspace(templateLib.domain(1),templateLib.domain(2),801)';
end
yProbe=eval_gap_template(templateLib,gQuery,xProbe);
nEdge=max(3,round(0.05*numel(xProbe)));
edgeBase=median([yProbe(1:nEdge);yProbe(end-nEdge+1:end)],'omitnan');
[~,iCenter]=max(abs(yProbe-edgeBase));
xCenter=xProbe(iCenter);
xInv = xNom;
valid = isfinite(xNom) & isfinite(vObs);
% Solve both monotone branches and retain the root closest to the nominal
% coordinate.  A vibration can move a sample across the pulse extremum, so
% selecting a branch from xNom alone creates a systematic inverse bias.
xDomain=[min(xProbe),max(xProbe)]; epsCenter=max(1e-7,1e-5*diff(xDomain));
for i = find(valid).'
    lo0=max(xDomain(1),xNom(i)-halfWidth); hi0=min(xDomain(2),xNom(i)+halfWidth);
    intervals=[lo0,min(hi0,xCenter-epsCenter);max(lo0,xCenter+epsCenter),hi0]; roots=nan(2,1);
    for ib=1:2
        lo=intervals(ib,1); hi=intervals(ib,2); if hi<=lo,continue;end
        fLo=eval_gap_template(templateLib,gQuery,lo)-vObs(i); fHi=eval_gap_template(templateLib,gQuery,hi)-vObs(i);
        if ~isfinite(fLo)||~isfinite(fHi)||fLo*fHi>0,continue;end
        for k=1:max(12,2*maxIter)
            mid=0.5*(lo+hi); fMid=eval_gap_template(templateLib,gQuery,mid)-vObs(i);
            if ~isfinite(fMid),break;end
            if abs(fMid)<=solveTol
                lo=mid; hi=mid; break;
            elseif fLo*fMid<=0
                hi=mid;
            else
                lo=mid; fLo=fMid;
            end
        end
        roots(ib)=0.5*(lo+hi);
    end
    roots=roots(isfinite(roots)); if isempty(roots),valid(i)=false;else,[~,ir]=min(abs(roots-xNom(i)));xInv(i)=roots(ir);end
end

vFit = eval_gap_template(templateLib, gQuery, xInv);
dFit = eval_gap_derivative(templateLib, gQuery, xInv);
residual = vObs - vFit;
valid = valid & isfinite(vFit) & isfinite(dFit) & abs(dFit) > derivFloor & abs(xInv - xNom) <= halfWidth;
if isfinite(voltageTol)
    valid = valid & abs(residual) <= voltageTol;
end

uInv = xNom - xInv;
weight = abs(dFit).^2;
weight(~valid | ~isfinite(weight)) = 0;
if any(weight > 0)
    weight = weight ./ mean(weight(weight > 0));
end

inv = struct();
inv.x_nom = xNom;
inv.x_inv = xInv;
inv.u_inv = uInv;
inv.v_fit = vFit;
inv.v_residual = residual;
inv.sensitivity = dFit;
inv.weight = weight;
inv.valid = valid;
inv.valid_fraction = mean(valid);
inv.rms_voltage_residual = sqrt(mean(residual(valid).^2, 'omitnan'));
end

function value = get_cfg(cfg, name, defaultValue)
if isfield(cfg, name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
