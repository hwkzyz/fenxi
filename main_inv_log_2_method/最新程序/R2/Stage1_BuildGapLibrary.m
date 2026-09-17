function lib = Stage1_BuildGapLibrary(forceRebuild)
%STAGE1_BUILDGAPLIBRARY 构造并缓存 R2 间隙响应模板库。
% 该阶段只依赖静态响应面，不生成振动波形，也不执行辨识。
if nargin < 1, forceRebuild = false; end
here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(fullfile(root,'local_func'),'-begin'); addpath(root,'-begin');
outDir = fullfile(here,'r2_profiled_recovery','output','stages');
if ~exist(outDir,'dir'), mkdir(outDir); end
cacheFile = fullfile(outDir,'R2_gap_library.mat');
if exist(cacheFile,'file') && ~forceRebuild
    S = load(cacheFile,'lib'); lib = S.lib; return
end
ctx = load_inv_log_2_project_context(root);
trust = load_fixed_trust_domain(root);
% 当前模型定义固定为三阶 field-basis；直接写在阶段配置中，避免为一行
% 表格构造再增加一个包装函数。
modelDef = table("field_basis_inv_log_2", "field_basis", "inv_log_2", NaN, 3, ...
    'VariableNames', {'model_label','response_model','basis_name','parameter','basis_order'});
lib = make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell, ...
    NaN,ctx.cfgAna.xGridN,modelDef,trust);
lib.stage = 'R2_gap_library'; lib.projectRoot = root;
save(cacheFile,'lib','-v7.3');
fprintf('R2 间隙库已保存：%s\n',cacheFile);
end
