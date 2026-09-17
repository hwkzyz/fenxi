function Result = Run_R3_BlindSingleFrequency(nCases)
%RUN_R3_BLINDSINGLEFREQUENCY  兼容入口（推荐运行 Main_R3）。
% 主程序只负责“路径—检查—运行—回显”；留出间隙的生成和辨识算法
% 保留在 r3_heldout_recovery 中，便于单独检查和复现实验。
if nargin < 1
    nCases = 3;
end
validateattributes(nCases, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'nCases');

%% 1. 项目路径
here = fileparts(mfilename('fullpath'));       % .../最新程序/R3
projectRoot = fileparts(here);                 % .../最新程序
addpath(fullfile(here, 'r3_heldout_recovery'), '-begin');
addpath(fullfile(here, 'local_func'), '-begin');
addpath(fullfile(projectRoot, 'R2', 'r2_profiled_recovery'), '-begin');
addpath(fullfile(projectRoot, 'R4', 'r4_pc1_manifold_transfer', ...
    'blind_unknown_gap'), '-begin');

%% 2. 运行前检查
assert_latest_formal_path(projectRoot);

%% 3. 运行 R3 正式辨识
% 保留此文件是为了兼容旧的命令行调用；新的主流程直接调用后端。
Result = Run_R3_UnifiedBlindBackend(nCases);

%% 4. 简短回显
disp(Result.table);
end
