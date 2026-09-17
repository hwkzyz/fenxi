%% R3 主程序：同一响应面族内的留出间隙验证
% 使用方法：在 MATLAB 当前目录运行本文件，或直接运行 Main_R3。
% 只在本文件顶部修改案例数；路径、数据加载和算法由正式入口统一管理。
clc;
clear;
close all;

%% 分析设置
nCases = 3;                 % 要运行的留出间隙案例数（1~3）

%% 执行 R3 正式流程
% 这里直接调用唯一的 R3 数值后端；不再经过 Run_R3_BlindSingleFrequency
% 这一层纯转发包装。
here = fileparts(mfilename('fullpath'));
projectRoot = fileparts(here);
addpath(fullfile(here,'r3_heldout_recovery'),'-begin');
addpath(fullfile(here,'local_func'),'-begin');
addpath(fullfile(projectRoot,'local_func'),'-begin');
fprintf('R3：开始运行 %d 个留出间隙案例...\n', nCases);
Result = Run_R3_UnifiedBlindBackend(nCases);

%% 结果摘要
disp('R3：运行完成，结果如下：');
disp(Result.table);
fprintf('结果文件夹：%s\n', Result.outputDir);
