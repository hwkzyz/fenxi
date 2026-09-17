%% R4 主程序：多响应面 PC1 流形迁移盲验证
% 使用方法：在 MATLAB 当前目录运行本文件，或直接运行 Main_R4。
% 这里保留少量实验设计参数；R4 的 PC1/多响应面算法集中在唯一后端
% 文件中，避免拆成很多只转发数据的子程序。
clc;
clear;
close all;

% 将 R4 唯一数值入口目录加入路径；不再通过额外的 Run 包装层转发。
addpath(fullfile(fileparts(mfilename('fullpath')), ...
    'r4_pc1_manifold_transfer', 'blind_unknown_gap'), '-begin');

%% 分析设置
opts = struct();
opts.gridN = 401;           % 响应面 x 方向插值网格数
opts.maxLowCandidates = 3;  % 低速候选保留数

%% 执行 R4 正式流程
fprintf('R4：开始运行多响应面流形迁移盲验证...\n');
Result = Main_10_UnifiedBlindBackendValidation(opts);

%% 结果摘要
disp('R4：运行完成，案例结果如下：');
disp(Result.cases);
fprintf('结果文件夹：%s\n', Result.output_dir);
