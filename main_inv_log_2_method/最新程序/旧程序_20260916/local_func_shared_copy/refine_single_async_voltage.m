function fit=refine_single_async_voltage(map,lib,cfg,seed)
%REFINE_SINGLE_ASYNC_VOLTAGE Complete-voltage continuous-frequency refinement.
t=map.t_v(:);x=map.x_v(:);V=map.V_a(:);
[~,component]=max(seed.A_id);f0=seed.f_id(component);
a0=seed.A_id(component)*cos(seed.phi_id(component));
b0=seed.A_id(component)*sin(seed.phi_id(component));z0=[seed.g_used,seed.dx_used,a0,b0,f0];
fGrid=cfg.route30SingleGrid;fMin=min(fGrid);fMax=max(fGrid);
fh=get_field(cfg,'route30SingleRefineHalfWidthHz',20);
gh=get_field(cfg,'route30SingleGapHalfWidthMm',.12);
dh=get_field(cfg,'route30DxHalfWidthMm',.5);au=get_field(cfg,'route30AmplitudeUpperMm',1.5);
lb=[max(.05,seed.g_used-gh),-dh,-au,-au,max(fMin,f0-fh)];
ub=[seed.g_used+gh,dh,au,au,min(fMax,f0+fh)];
if ~(isfield(lib,'fixedPathIncrement')&&lib.fixedPathIncrement)
    lb(1)=max(lb(1),min(lib.gapTrain)-cfg.rawGapSearchMargin);
    ub(1)=min(ub(1),max(lib.gapTrain)+cfg.rawGapSearchMargin);
end
z0=min(max(z0,lb),ub);penalty=10*max(std(V),1e-3);
[z,si]=solve_lsq_bounded(@(q)residual(q,t,V,x,lib,penalty),z0,lb,ub,...
    get_field(cfg,'route30SingleMaxIter',160),1e-12,1e-12);
A=hypot(z(3),z(4));phi=atan2(z(4),z(3));
VFit=eval_gap_template(lib,z(1),x-z(2)-A*sin(2*pi*z(5)*t+phi));
r=V-VFit;good=isfinite(r);
fit=struct('g',z(1),'dx',z(2),'A',A,'phi',phi,'f',z(5),...
    'p',[A,phi,z(5)],'theta',z,'VFit',VFit,...
    'rmse',sqrt(mean(r(good).^2)),'solve_info',si,'seed_frequency_hz',f0);
end

function r=residual(z,t,V,x,lib,penalty)
A=hypot(z(3),z(4));phi=atan2(z(4),z(3));
r=V-eval_gap_template(lib,z(1),x-z(2)-A*sin(2*pi*z(5)*t+phi));
r(~isfinite(r))=penalty;
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
