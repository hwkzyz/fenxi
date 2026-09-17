function Result = Run_R3_UnifiedBlindBackend(nCases)
%RUN_R3_UNIFIEDBLINDBACKEND  R3 数值后端。
% 主流程按“配置→留出数据→盲辨识→汇总”展开；核心数学仍由
% local_func 中的函数完成，便于单独测试和复用。
if nargin < 1, nCases = 3; end

%% 1. 项目路径和 R3 协议
root = fileparts(mfilename('fullpath'));
mainDir = fileparts(fileparts(root));
addpath(mainDir, '-begin');
addpath(fullfile(fileparts(root), 'local_func'), '-begin');
addpath(fullfile(root, 'local_func'), '-begin');
addpath(fullfile(mainDir, 'local_func'), '-begin');
% R3 固定试验协议。协议只有常量配置，直接放在主流程，避免再跳转
% 到一个只返回结构体的 R3_Protocol 子程序。
P0 = struct();
P0.gAllMm = 0.2:0.1:1.5;
P0.gCalMm = [0.2 0.3 0.4 0.6 0.8 1.1 1.3 1.4 1.5];
P0.gOperatingMm = [0.5 0.7 0.9 1.0 1.2];
P0.gReferenceMm = 0.8;
P0.gMainMm = [P0.gOperatingMm P0.gReferenceMm];
P0.gapBoundsMm = [0.2 1.5];
P0.adaptiveSearchBoundsMm = [0.4 1.3];
P0.adaptiveCoarseGapGridMm = 0.4:0.1:1.3;
P0.adaptiveSearchScope = "tested numerical search range, not an engineering maximum";
P0.amplitudeMainMm = [0.10 0.15 0.20 0.25 0.30 0.37];
P0.amplitudeNullMm = 0;
P0.eoTrue = 10;
P0.deltaFTrueHz = 0;
P0.frequencyRangeHz = [300 1500];
P0.phaseValuesRad = 2*pi*(((1:75)-0.5)/75);
P0.referenceSnrDb = 15;
P0.referenceAmplitudeMm = 0.10;
P0.referencePhaseRad = pi/4;
P0.amplitudeAbsToleranceMm = 0.020;
P0.amplitudeRelTolerance = 0.10;
P0.phaseToleranceRad = deg2rad(20);
P0.gapToleranceMm = 0.020;
P0.mainLowSeed = 2026082501;
P0.mainHighSeeds = 2026082600 + (1:75);
% 结果行字段固定，直接写出，避免再增加 row0 这种纯初始化函数。
emptyRow = struct('case_id',NaN,'eo_true',NaN,'eo_hat',NaN, ...
    'gap_true_mm',NaN,'gap_hat_mm',NaN,'gap_error_mm',NaN, ...
    'A_true_mm',NaN,'A_hat_mm',NaN,'delta_f_true_hz',NaN, ...
    'delta_f_hat_hz',NaN,'rmse_V',NaN,'status',"");
rows = repmat(emptyRow, 0, 1);

%% 2. 逐个留出间隙案例生成观测并执行盲辨识
for ic=1:min(nCases,numel(P0.gOperatingMm))
 P = P0;
 R3 = r3_build_context(mainDir, P);
 [C0, ~] = r3_freeze_reference_calibration(R3, 0, P.mainLowSeed + ic);
 g = P.gOperatingMm(ic);
 A = 0.20;
 phi = P.referencePhaseRad;
 [~, map, truth] = r3_generate_high_case(R3, g, A, phi, 0, P.mainHighSeeds(ic));
 % R3 truth is the raw COMSOL curve at a held-out gap node.  Do not rebuild
 % the voltage from the calibration interpolant here: that would turn R3
 % into an in-sample response-surface test and defeat the held-out design.
 bundle = struct('t',map.t_v,'x',map.x_v,'V',map.V_a, ...
     'theta',map.theta_v,'sensorId',map.S_v,'turnId',map.rev_v, ...
     'Wvp',ones(size(map.V_a)),'Wfull',ones(size(map.V_a)));
 response = struct('evalF',@(gg,x)eval_gap_template(C0.templateModel,gg,x), ...
     'evalFx',@(gg,x)eval_gap_derivative(C0.templateModel,gg,x));
 gb = R3.P.adaptiveSearchBoundsMm;
 bc = struct('gapBounds',gb,'gapGrid',linspace(gb(1),gb(2),31), ...
     'eoGrid',1:30,'topK',5,'rotHz',R3.cfg.RPM_high/60, ...
     'minGradient',1e-6,'minSamples',30,'dxBound',.2, ...
     'ampCoeffBound',1.5,'deltaFBound',2,'maxIter',250, ...
     'ambiguityMargin',.01);
 B = run_unified_blind_dynamic_backend(bundle,response,bc);
 b = B.best;
 rows(end+1) = struct('case_id',ic,'eo_true',truth.eo_true, ...
     'eo_hat',b.eo,'gap_true_mm',g,'gap_hat_mm',b.gap, ...
     'gap_error_mm',b.gap-g,'A_true_mm',A,'A_hat_mm',b.A, ...
     'delta_f_true_hz',truth.delta_f_true_hz,'delta_f_hat_hz',b.deltaF, ...
     'rmse_V',b.fullWaveRmseV,'status',B.status); %#ok<AGROW>
end
%% 3. 保存并返回结果
T = struct2table(rows);
out = fullfile(root,'output','unified_blind_backend');
if ~exist(out,'dir'), mkdir(out); end
writetable(T,fullfile(out,'r3_unified_blind_backend.csv'));
save(fullfile(out,'r3_unified_blind_backend.mat'),'T','-v7.3');
Result = struct('table',T,'outputDir',out);
disp(T);
end
