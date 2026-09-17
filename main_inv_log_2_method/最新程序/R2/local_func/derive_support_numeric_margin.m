function M=derive_support_numeric_margin(templateLib,cfg)
%DERIVE_SUPPORT_NUMERIC_MARGIN Audit the x support consumed by operators.
% Values are operator half-widths, not empirical safety factors.
if nargin<2||isempty(cfg),cfg=struct();end
fixedPath=isfield(templateLib,'fixedPathIncrement')&&templateLib.fixedPathIncrement;
if fixedPath
    templateDerivative=get_field(cfg,'supportTemplateDerivativeHalfWidthMm',1e-4);
else
    templateDerivative=get_field(cfg,'supportTemplateDerivativeHalfWidthMm',0);
end
staticRx=get_field(cfg,'supportStaticRxHalfWidthMm',0);
staticRgX=get_field(cfg,'supportStaticRgSpatialHalfWidthMm',0);
interpGuard=get_field(cfg,'supportInterpolationGuardWidthMm',0);
values=[templateDerivative staticRx staticRgX interpGuard];
if any(~isfinite(values))||any(values<0)
    error('support:InvalidNumericMargin','Numerical support half-widths must be finite and nonnegative.');
end
M=struct('template_derivative_half_width_mm',templateDerivative,...
    'static_rx_half_width_mm',staticRx,...
    'static_rg_spatial_half_width_mm',staticRgX,...
    'interpolation_guard_width_mm',interpGuard,...
    'numeric_margin_mm',max(values),...
    'rule','maximum x-direction half-width of all forward-model operators');
end
function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
