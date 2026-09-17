function pathModel = make_low_increment_path_model(templateLib, g0, templateLow)
%MAKE_LOW_INCREMENT_PATH_MODEL Freeze a low-speed baseline for delta-g use.
% Use eval_path_increment_template for delta_g queries.  eval_gap_template
% on this path model retains its legacy absolute-gap query contract.
if ~isscalar(g0) || ~isfinite(g0), error('g0 must be a finite scalar.'); end
templateLib = prepare_gap_template_library(templateLib);
if nargin < 3 || isempty(templateLow)
    templateLow = eval_gap_template(templateLib,g0,templateLib.xGrid(:));
end
pathCal = struct('g0',g0,'mu',0,'tau',0,'zeta',1,'kappa',1);
pathModel = build_path_template_model(templateLib,templateLow,pathCal);
end
