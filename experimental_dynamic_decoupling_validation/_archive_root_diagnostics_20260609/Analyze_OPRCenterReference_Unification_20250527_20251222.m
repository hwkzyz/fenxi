%% Analyze OPR center reference unification for 20250527 and 20251222
% This script audits whether low-speed rotating calibration and high-speed
% legacy extraction use the same OPR pulse-center coordinate reference.
%
% It does not modify any calibration or identification result.

clear; clc;

rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir, 'analysis_outputs', 'opr_center_reference');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

% Shared geometry in these two data sets.
pulsesPerRev = 6;
rTipMm = 62.0;

cases = struct([]);
cases(end+1).dataset = '20250527';
cases(end).sensorTag = 'S136';
cases(end).templateFile = fullfile(rootDir, '20250527_low_speed_rotating_calibration', ...
    'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S136_GradientXRange030_OPRCenterStd_20250527.mat');
cases(end).dynamicFile = fullfile(rootDir, '20250527_low_speed_rotating_calibration', ...
    'output', 'dynamic_maps', ...
    'DynamicMap_B1_S136_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20250527.mat');
cases(end).highSpeedOprFile = fullfile(rootDir, '20250527', 'legacy', 'output', ...
    '20250526_2500-3500_t400', 'jiluOPR.mat');
cases(end).gapHighMapFile = fullfile(rootDir, '20250527_low_speed_gap_prior_decoupling', ...
    'outputs', 'Step06_HighMap_20250527_B1_S136.mat');

cases(end+1).dataset = '20251222';
cases(end).sensorTag = 'S123';
cases(end).templateFile = fullfile(rootDir, '20251222_low_speed_rotating_calibration', ...
    'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S123_GradientXRange030_OPRCenterStd_20251222.mat');
cases(end).dynamicFile = fullfile(rootDir, '20251222_low_speed_rotating_calibration', ...
    'output', 'dynamic_maps', ...
    'DynamicMap_B1_S123_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20251222.mat');
cases(end).highSpeedOprFile = fullfile(rootDir, '20251222', 'legacy', 'output', ...
    '1000_2500_3500', 'jiluOPR.mat');
cases(end).gapHighMapFile = '';

rows = {};
for ic = 1:numel(cases)
    C = cases(ic);

    rows(end+1, :) = make_saved_reference_row_local( ...
        C.dataset, 'low_speed_template', C.templateFile, 'Template', ...
        pulsesPerRev, rTipMm); %#ok<SAGROW>

    rows(end+1, :) = make_saved_reference_row_local( ...
        C.dataset, 'low_speed_dynamic_map', C.dynamicFile, 'DynamicMap', ...
        pulsesPerRev, rTipMm); %#ok<SAGROW>

    rows(end+1, :) = make_jilu_reference_row_local( ...
        C.dataset, 'high_speed_legacy_jiluOPR', C.highSpeedOprFile, ...
        pulsesPerRev, rTipMm, 'computed_from_jiluOPR'); %#ok<SAGROW>

    if ~isempty(C.gapHighMapFile)
        rows(end+1, :) = make_highmap_reference_row_local( ...
            C.dataset, 'gap_prior_highMap_existing_output', C.gapHighMapFile, ...
            pulsesPerRev, rTipMm); %#ok<SAGROW>
    end
end

T = cell2table(rows, 'VariableNames', { ...
    'dataset', 'source', 'status', 'filePath', 'jiluColumns', ...
    'centerMinusStartUs', 'pulseWidthUs', 'phaseShiftDeg', 'phaseShiftMm', ...
    'standardAngleReference', 'comment'});

disp(T(:, {'dataset', 'source', 'status', 'centerMinusStartUs', ...
    'phaseShiftDeg', 'phaseShiftMm', 'standardAngleReference', 'comment'}));

outCsv = fullfile(outDir, 'OPR_Center_Reference_Audit_20250527_20251222.csv');
writetable(T, outCsv);
fprintf('\nSaved audit table: %s\n', outCsv);

%% Local helper functions
function row = make_saved_reference_row_local(dataset, sourceName, matFile, varName, pulsesPerRev, rTipMm)
if exist(matFile, 'file') ~= 2
    row = empty_row_local(dataset, sourceName, 'missing_file', matFile);
    return;
end
S = load(matFile);
if ~isfield(S, varName)
    row = empty_row_local(dataset, sourceName, 'missing_variable', matFile);
    row{end} = sprintf('Variable %s not found.', varName);
    return;
end
Obj = S.(varName);
if isfield(Obj, 'OPRReference')
    ref = Obj.OPRReference;
    phaseDeg = read_field_or_nan_local(ref, 'phase_shift_deg');
    phaseMm = read_field_or_nan_local(ref, 'phase_shift_mm');
    centerMinusStartUs = read_field_or_nan_local(ref, 'center_minus_start_s') * 1e6;
    pulseWidthUs = NaN;
    if isfield(ref, 'standard_angle_reference')
        angleRef = string(ref.standard_angle_reference);
    else
        angleRef = "opr_pulse_center_inferred_from_program";
    end
    row = {dataset, sourceName, 'ok_opr_reference_field', matFile, NaN, ...
        centerMinusStartUs, pulseWidthUs, phaseDeg, phaseMm, char(angleRef), ...
        'Saved result carries OPRReference; current readers use OPR-center standard angles directly and only convert legacy start-edge fields.'};
    return;
end
if isfield(Obj, 'jiluOPR')
    row = make_jilu_reference_row_local(dataset, sourceName, matFile, ...
        pulsesPerRev, rTipMm, 'computed_from_saved_jiluOPR');
    return;
end
row = empty_row_local(dataset, sourceName, 'no_opr_reference_field', matFile);
row{end} = 'Output exists but does not expose OPRReference.';
end

function row = make_jilu_reference_row_local(dataset, sourceName, matFile, pulsesPerRev, rTipMm, commentPrefix)
if exist(matFile, 'file') ~= 2
    row = empty_row_local(dataset, sourceName, 'missing_file', matFile);
    return;
end
S = load(matFile, 'jiluOPR', 'metadata');
if ~isfield(S, 'jiluOPR')
    row = empty_row_local(dataset, sourceName, 'missing_jiluOPR', matFile);
    return;
end
jiluOPR = S.jiluOPR;
ref = build_opr_reference_from_jilu_local(jiluOPR, pulsesPerRev, rTipMm);
comment = sprintf('%s; jiluOPR(:,1)=center, (:,2)=legacy start, (:,3)=legacy end if three columns.', commentPrefix);
if size(jiluOPR, 2) < 3
    comment = [comment, ' Only two columns found; interpretation is less certain.'];
end
row = {dataset, sourceName, 'ok_computed_from_jiluOPR', matFile, size(jiluOPR, 2), ...
    ref.center_minus_start_s * 1e6, ref.pulse_width_s * 1e6, ...
    ref.phase_shift_deg, ref.phase_shift_mm, 'opr_pulse_center_if_standard_angle_shifted', comment};
end

function row = make_highmap_reference_row_local(dataset, sourceName, matFile, ~, ~)
if exist(matFile, 'file') ~= 2
    row = empty_row_local(dataset, sourceName, 'missing_file', matFile);
    return;
end
S = load(matFile, 'highMap');
if ~isfield(S, 'highMap')
    row = empty_row_local(dataset, sourceName, 'missing_highMap', matFile);
    return;
end
H = S.highMap;
if isfield(H, 'experiment_case') && isfield(H.experiment_case, 'Diagnostics') && ...
        isfield(H.experiment_case.Diagnostics, 'OPRReference')
    ref = H.experiment_case.Diagnostics.OPRReference;
    angleRef = read_string_field_local(H.experiment_case.Diagnostics, ...
        'standard_angle_reference', 'unknown');
    row = {dataset, sourceName, 'ok_highMap_has_diagnostics', matFile, NaN, ...
        read_field_or_nan_local(ref, 'center_minus_start_s') * 1e6, ...
        read_field_or_nan_local(ref, 'pulse_width_s') * 1e6, ...
        read_field_or_nan_local(ref, 'phase_shift_deg'), ...
        read_field_or_nan_local(ref, 'phase_shift_mm'), angleRef, ...
        'Existing highMap was generated after helper diagnostics were available.'};
    return;
end
row = empty_row_local(dataset, sourceName, 'existing_output_lacks_opr_diagnostics', matFile);
row{end} = ['This output cannot prove whether x_points used OPR-center standard-angle conversion. ' ...
    'Rebuild Step06 after helper update before comparing final gap/vibration results.'];
end

function ref = build_opr_reference_from_jilu_local(jiluOPR, pulsesPerRev, rTipMm)
ref = struct('center_minus_start_s', NaN, 'pulse_width_s', NaN, ...
    'phase_shift_deg', NaN, 'phase_shift_mm', NaN);
if size(jiluOPR, 2) < 2 || pulsesPerRev < 1
    return;
end
centerTime = jiluOPR(:, 1);
startTime = jiluOPR(:, 2);
if size(jiluOPR, 2) >= 3
    endTime = jiluOPR(:, 3);
else
    endTime = NaN(size(centerTime));
end
validSpan = (1:numel(centerTime) - pulsesPerRev).';
if isempty(validSpan)
    return;
end
revPeriod = centerTime(validSpan + pulsesPerRev) - centerTime(validSpan);
centerMinusStart = centerTime(validSpan) - startTime(validSpan);
pulseWidth = endTime(validSpan) - startTime(validSpan);
valid = isfinite(revPeriod) & revPeriod > 0 & isfinite(centerMinusStart);
if ~any(valid)
    return;
end
phaseShiftDeg = centerMinusStart(valid) ./ revPeriod(valid) * 360;
ref.center_minus_start_s = median(centerMinusStart(valid), 'omitnan');
ref.pulse_width_s = median(pulseWidth(isfinite(pulseWidth)), 'omitnan');
ref.phase_shift_deg = median(phaseShiftDeg, 'omitnan');
ref.phase_shift_mm = ref.phase_shift_deg * (pi / 180) * rTipMm;
end

function val = read_field_or_nan_local(S, fieldName)
val = NaN;
if isstruct(S) && isfield(S, fieldName)
    raw = S.(fieldName);
    if isnumeric(raw) && isscalar(raw)
        val = raw;
    end
end
end

function txt = read_string_field_local(S, fieldName, defaultText)
txt = defaultText;
if isstruct(S) && isfield(S, fieldName)
    raw = S.(fieldName);
    if isstring(raw) || ischar(raw)
        txt = char(raw);
    end
end
end

function row = empty_row_local(dataset, sourceName, statusText, matFile)
row = {dataset, sourceName, statusText, matFile, NaN, NaN, NaN, NaN, NaN, ...
    'unknown', ''};
end
