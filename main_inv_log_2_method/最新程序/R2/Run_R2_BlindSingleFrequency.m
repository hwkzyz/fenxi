function Result = Run_R2_BlindSingleFrequency(nCases)
%RUN_R2_BLINDSINGLEFREQUENCY  R2 兼容入口，转到分阶段流程。
%
% 代码按试验程序的阅读顺序组织：
%   1) 设置案例数；2) 设置项目路径；3) 检查正式路径；
%   4) 调用一次 R2 计算；5) 返回并显示结果。
% 具体数值算法仍放在 r2_profiled_recovery/local_func 中，避免主程序
% 被数值细节淹没，也避免在不同入口复制同一段算法。
if nargin < 1
    nCases = 6;
end
validateattributes(nCases, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'nCases');

%% 分阶段运行（旧的统一后端不再作为正式入口）
here = fileparts(mfilename('fullpath'));
addpath(here,'-begin');
Stage1_BuildGapLibrary(false);
for k=1:nCases
    Stage2_GenerateWaveform(k,false);
    Stage3_IdentifyVibration(k,false);
end
Result.table = Stage4_AnalyzeResults(1:nCases);
Result.outputDir = fullfile(here,'r2_profiled_recovery','output','stages');
end
