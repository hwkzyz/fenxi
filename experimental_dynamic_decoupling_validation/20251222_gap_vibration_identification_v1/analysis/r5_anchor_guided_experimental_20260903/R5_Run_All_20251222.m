function result = R5_Run_All_20251222(varargin)
%R5_RUN_ALL_20251222 Sidecar R5 smoke/entry point. Production files are read-only.
p = inputParser; addParameter(p,'DynamicMapFile',''); addParameter(p,'OutputDir',''); addParameter(p,'WindowId',[]); addParameter(p,'TargetBlade',[]); addParameter(p,'SensorId',[]); addParameter(p,'RunStaticValidation',true); parse(p,varargin{:});
sidecar = fileparts(mfilename('fullpath')); addpath(genpath(sidecar));
cfg = R5_Config_20251222();
if ~isempty(p.Results.OutputDir)
    cfg.outputs.root = p.Results.OutputDir;
    cfg.paths.r5Results = p.Results.OutputDir;
    cfg.paths.r5Tables = fullfile(p.Results.OutputDir,'tables');
    cfg.paths.r5Figures = fullfile(p.Results.OutputDir,'figures');
end
if ~exist(cfg.outputs.root,'dir'), mkdir(cfg.outputs.root); end
S = R5_Load_Main05_Surface(cfg.inputs.responseSurface);
T = R5_Load_Main07_LowSpeedTemplate(cfg.inputs.lowSpeedTemplate);
[dmFile, dmCandidates] = R5_Find_Main08_DynamicMap(cfg, p.Results.DynamicMapFile);
if ~isempty(dmFile), D = R5_Load_Main08_DynamicMap(dmFile); else, D = []; end
contract = R5_Validate_InputContract(cfg,S,T,dmFile);
staticMetrics = struct('status','skipped');
if p.Results.RunStaticValidation
    staticMetrics.sameBlade = R5_V1_SameBlade_Migration(S,cfg.paths.r5Tables);
    staticMetrics.leaveGap = R5_V2_LeaveGap_Migration(S,cfg.paths.r5Tables);
    staticMetrics.leaveBlade = R5_V4_LeaveBlade(S,cfg.paths.r5Tables);
    staticMetrics.status = 'completed';
end
% The first executable R5 stage is deliberately a static anchor-to-state audit.
target = cfg.case.targetBlade; if ~isempty(p.Results.TargetBlade), target=p.Results.TargetBlade; end
sensorId = cfg.case.analysisSensors(1); if ~isempty(p.Results.SensorId), sensorId=p.Results.SensorId; end
anchor = R5_Make_Main07_Anchor(T,target,sensorId);
if target == 2
    axisType = 'physical_gap_mm'; stateValues = S.physical_gap_mm;
else
    axisType = 'effective_state'; stateValues = 1:size(S.waveforms_mv,3);
end
L = R5_Build_TransitionLibrary(S,axisType,stateValues);
candidateOpts = struct('topKDonor',cfg.donor.topK, ...
    'distanceFloor',cfg.donor.distanceFloor, ...
    'maxManifoldDistance',cfg.donor.manifoldMaxDistance);
C = R5_Generate_CandidateStates(anchor,L,stateValues, candidateOpts);
plotFile = fullfile(cfg.outputs.root,'figures','R5_anchor_candidate_waveforms.png');
if ~exist(fileparts(plotFile),'dir'), mkdir(fileparts(plotFile)); end
P = R5_Plot_CandidateStates(anchor.x_mm,anchor.y_mv,C,struct('stateValue',NaN),plotFile);
profile = []; profilePlot = '';
profileCurvePlot = '';
gate = struct('status','not_run');
if ~isempty(D) && ~isempty(p.Results.WindowId)
    if isfield(D,'TargetBlade') && isfinite(double(D.TargetBlade)) && double(D.TargetBlade) ~= target
        error('R5:TargetBladeMismatch','DynamicMap target B%d does not match requested target B%d.',double(D.TargetBlade),target);
    end
    window = R5_Make_Main08_Window(D,p.Results.WindowId,sensorId,T);
    % No extrapolation: the high-speed objective uses only the migrated
    % waveform's common support and reports how much of Main08 was retained.
    fitDomain = anchor.x_mm(anchor.supportMask);
    inDomain = window.x_mm >= min(fitDomain) & window.x_mm <= max(fitDomain);
    window.inputPointCount = numel(window.x_mm);
    window.validMask = window.validMask & inDomain;
    window.fitCoverageFraction = nnz(window.validMask) / window.inputPointCount;
    profileOpts = struct('frequencyGridHz',cfg.profile.frequencyGridHz, ...
        'dxBoundsMm',cfg.profile.dxBoundsMm,'amplitudeBoundsMm',cfg.profile.amplitudeBoundsMm, ...
        'phaseBoundsRad',cfg.profile.phaseBoundsRad,'maxIterations',cfg.profile.maxIterations, ...
        'minSupportFraction',cfg.gates.minSupportFraction,'maxProfileSamples',600);
    profile = R5_ProfileFit_StateVibration(window,C,profileOpts);
    gate = R5_IdentifiabilityGate_StateVibration(profile,C,cfg);
    profilePlot = fullfile(cfg.outputs.root,'figures',sprintf('R5_profile_window_%03d.png',p.Results.WindowId));
    if isfield(profile.selected,'predicted_mv')
        R5_Plot_ProfileFit(profile.selected,profilePlot);
    else
        profilePlot = '';
    end
    profileCurvePlot = fullfile(cfg.outputs.root,'figures',sprintf('R5_profile_objective_window_%03d.png',p.Results.WindowId));
    R5_Plot_ProfileCurve(profile,profileCurvePlot);
end
status = 'low_speed_anchor_candidate_ready';
if ~isempty(profile) && isfield(gate,'status') && strcmp(gate.status,'identifiable_window_only'), status='high_speed_profile_identifiable_window_only'; end
if ~isempty(profile) && isfield(gate,'status') && strcmp(gate.status,'diagnostic_only'), status='high_speed_profile_diagnostic_only'; end
result = struct('status',status,'config',cfg,'contract',contract, ...
    'dynamicMapFile',dmFile,'dynamicMapCandidates',{dmCandidates},'dynamicMapLoaded',~isempty(D), ...
    'anchor',anchor,'transitions',L,'candidates',C,'plotFile',P, ...
    'profile',profile,'identifiabilityGate',gate,'profilePlotFile',profilePlot, ...
    'profileCurvePlotFile',profileCurvePlot, ...
    'staticMigrationMetrics',staticMetrics);
save(fullfile(cfg.outputs.root,'R5_static_anchor_result.mat'),'result','-v7.3');
end
