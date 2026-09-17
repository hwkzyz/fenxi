%% R2 主程序：同一响应面上的盲间隙/振动辨识
% 使用方法：在 MATLAB 当前目录运行本文件，或直接运行 Main_R2。
% 参数集中放在本文件顶部，阶段按理论顺序独立运行；数值核心函数
% 只在对应阶段直接调用，不再通过多层“主程序→包装程序→后端”转发。
clc; clear; close all;
%% 运行设置（library / waveform / identify / visualize / all）
stage = "all";
caseID = 1;
nCases = 6;
switch lower(stage)
    case "library"
        Result = Stage1_BuildGapLibrary(false);
    case "waveform"
        Result = Stage2_GenerateWaveform(caseID,false);
    case "identify"
        for k=1:nCases, Result = Stage3_IdentifyVibration(k,false); end
    case "visualize"
        Result = Stage4_AnalyzeResults(1:nCases);
    case "all"
        Stage1_BuildGapLibrary(false);
        for k=1:nCases, Stage2_GenerateWaveform(k,false); end
        for k=1:nCases, Result = Stage3_IdentifyVibration(k,false); end
        Stage4_AnalyzeResults(1:nCases);
    otherwise
        error('stage 必须为 library、waveform、identify、visualize 或 all');
end
