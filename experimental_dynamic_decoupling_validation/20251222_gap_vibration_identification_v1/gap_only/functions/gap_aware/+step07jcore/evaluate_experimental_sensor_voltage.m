function [model, info] = evaluate_experimental_sensor_voltage(Tpl, responseSurface, corr, dg, xQuery, dmu, dtau)
%EVALUATE_EXPERIMENTAL_SENSOR_VOLTAGE Exact diagnostic copy of Main10's F_s.
% Units and clipping follow Main10: voltage is mV relative to baseline,
% xQuery is mm, and dg/dmu/dtau are the calibrated sensor increments.
if nargin < 6, dmu = 0; end
if nargin < 7, dtau = 0; end
xQuery = xQuery(:);
xWarpRaw = xQuery - dtau;
xLo = min(Tpl.x_grid(:)); xHi = max(Tpl.x_grid(:));
xLow = min(max(xWarpRaw,xLo),xHi);
overLow = max(xLo-xWarpRaw,0) + max(xWarpRaw-xHi,0);
vLow = interp1(Tpl.x_grid(:),(Tpl.v_grid(:)-Tpl.baseline)*1000,xLow,'pchip',NaN);
if abs(dg) <= eps && abs(dmu) <= eps
    model = vLow; over = overLow;
else
    [fd,od] = local_raw(responseSurface,corr.g0Mm+dg,corr.tauMm, ...
        corr.xScale,corr.muGapPerXMm+dmu,xWarpRaw);
    [fb,ob] = local_raw(responseSurface,corr.g0Mm,corr.tauMm, ...
        corr.xScale,corr.muGapPerXMm,xWarpRaw);
    model = vLow + corr.voltageGain .* (fd-fb);
    over = max(max(od,ob),overLow);
end
info = struct('overshoot',over,'querySafe',over <= 0);
end

function [f,over] = local_raw(rs,g0,tau,k,mu,x)
xRaw = k .* (x(:)-tau); gRaw = g0 + mu .* (x(:)-tau);
xMin=min(rs.xGrid(:)); xMax=max(rs.xGrid(:));
gMin=min(rs.gTrainMm(:)); gMax=max(rs.gTrainMm(:));
xq=min(max(xRaw,xMin),xMax); gq=min(max(gRaw,gMin),gMax);
over=hypot(max(xMin-xRaw,0)+max(xRaw-xMax,0), ...
    max(gMin-gRaw,0)+max(gRaw-gMax,0));
f=local_response(rs,gq,xq);
end

function f = local_response(rs,g,xq)
b0=local_interp(rs.xGrid,rs.coeff(:,1),xq);
b1=local_interp(rs.xGrid,rs.coeff(:,2),xq);
b2=local_interp(rs.xGrid,rs.coeff(:,3),xq);
f=b0+b1./g+b2.*log(g./rs.g0Mm);
end

function yq = local_interp(xg,yg,xq)
ok=isfinite(xg(:)) & isfinite(yg(:));
if nnz(ok)<2, yq=nan(size(xq)); return; end
x=xg(ok); y=yg(ok); xq=min(max(xq,min(x)),max(x));
yq=interp1(x,y,xq,'linear',NaN);
end
