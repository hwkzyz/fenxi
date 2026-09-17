function C=derive_support_aware_domain(highMap,lowMap,templateLib,physical,cfg)
%DERIVE_SUPPORT_AWARE_DOMAIN Certify the joint low/high forward-model domain.
% physical requires dx_bounds_mm and amplitude_max_mm.  The latter are
% physical amplitude bounds, not numerical coefficient bounds.
if nargin<5||isempty(cfg),cfg=struct();end
required={'dx_bounds_mm','amplitude_max_mm'};
for i=1:numel(required),if ~isfield(physical,required{i}),error('support:MissingPhysicalBound','physical.%s is required.',required{i});end,end
d=sort(physical.dx_bounds_mm(:).');a=physical.amplitude_max_mm(:).';
if numel(d)~=2||any(~isfinite(d))||any(~isfinite(a))||any(a<0),error('support:InvalidPhysicalBound','Physical support bounds must be finite.');end
u=sum(a);
if isfield(physical,'numeric_margin_mm')&&~isempty(physical.numeric_margin_mm)
    margin=struct('numeric_margin_mm',physical.numeric_margin_mm,'rule','explicit audited operator margin');
else
    margin=derive_support_numeric_margin(templateLib,cfg);
end
m=margin.numeric_margin_mm;
xh=support_interval(highMap);xl=support_interval(lowMap);xr=static_x_support(templateLib);
s=[max(xl(1),xr(1)) min(xl(2),xr(2))];safe=[s(1)+m s(2)-m];
q=[xh(1)-d(2)-u xh(2)-d(1)+u];
xcert=[safe(1)+d(2)+u safe(2)+d(1)-u];
Cgap=gap_support(templateLib,physical);
C=struct('high_observed_mm',xh,'low_observed_mm',xl,'static_x_support_mm',xr,...
    'common_support_mm',s,'safe_support_mm',safe,'reachable_query_mm',q,...
    'high_certified_mm',xcert,'dx_bounds_mm',d,'amplitude_max_mm',a,...
    'u_max_mm',u,'numeric_margin_mm',m,'numeric_margin_audit',margin,...
    'static_g_support_mm',Cgap.static,'physical_g_bounds_mm',Cgap.physical,...
    'gap_domain_supported',Cgap.supported,...
    'full_high_domain_supported',q(1)>=safe(1)&&q(2)<=safe(2),...
    'certified_domain_nonempty',xcert(1)<=xcert(2));
end
function x=support_interval(map)
if isfield(map,'S_v')&&numel(map.S_v)==numel(map.x_v)
    ids=unique(map.S_v(:));limits=nan(numel(ids),2);
    for i=1:numel(ids),q=map.S_v(:)==ids(i);limits(i,:)=[min(map.x_v(q)) max(map.x_v(q))];end
    x=[max(limits(:,1)) min(limits(:,2))];
else
    x=[min(map.x_v) max(map.x_v)];
end
end
function x=static_x_support(lib)
if isfield(lib,'pathCache')&&isfield(lib.pathCache,'xGrid')
    grid=lib.pathCache.xGrid;
elseif isfield(lib,'xGrid')
    grid=lib.xGrid;
else
    error('support:MissingStaticXSupport','Template library has no x support metadata.');
end
x=[min(grid) max(grid)];
end
function G=gap_support(lib,physical)
if isfield(lib,'pathCache')&&isfield(lib.pathCache,'gGrid')
    grid=lib.pathCache.gGrid;
elseif isfield(lib,'gapTrain')
    grid=lib.gapTrain;
else
    grid=[-Inf Inf];
end
gs=[min(grid) max(grid)];gp=get_field(physical,'gap_bounds_mm',[NaN NaN]);
ok=all(isfinite(gp))&&numel(gp)==2;
if ok,gp=sort(gp(:).');supported=gp(1)>=gs(1)&&gp(2)<=gs(2);else,supported=NaN;end
G=struct('static',gs,'physical',gp,'supported',supported);
end
function v=get_field(s,n,d),if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end,end
