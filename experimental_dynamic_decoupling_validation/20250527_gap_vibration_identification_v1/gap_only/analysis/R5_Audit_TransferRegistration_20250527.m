function A = R5_Audit_TransferRegistration_20250527(r5File, legacyCsv, outputFile)
%R5_AUDIT_TRANSFERREGISTRATION_20250527 Compare R5 and legacy registrations.
% This is an audit only. It does not substitute legacy parameters into R5.
cfg = Config_20250527();
if nargin < 1 || isempty(r5File)
    r5File = fullfile(cfg.paths.results, 'r5_b2_anchor_registration.mat');
end
if nargin < 2 || isempty(legacyCsv)
    legacyCsv = cfg.files.legacyP2Registration;
end
if nargin < 3 || isempty(outputFile)
    outputFile = fullfile(cfg.paths.results, 'r5_vs_legacy_registration_audit.csv');
end
S = load(r5File, 'G');
G = S.G;
L = readtable(legacyCsv, 'TextType', 'string');
L = L(strcmpi(string(L.model), 'P2_registration'), :);
rows = cell(numel(G.registration), 1);
for k = 1:numel(G.registration)
    sid = G.registration(k).sensor_id;
    j = find(L.sensorId == sid, 1);
    assert(~isempty(j), 'R5:LegacyRegistrationMissing', 'No legacy registration for S%d.', sid);
    rows{k} = table(sid, G.registration(k).tau_mm, L.tauMm(j), ...
        G.registration(k).tau_mm - L.tauMm(j), ...
        G.registration(k).x_scale, L.xScale(j), ...
        G.registration(k).x_scale - L.xScale(j), ...
        G.registration(k).voltage_gain, L.voltageGain(j), ...
        G.registration(k).voltage_gain - L.voltageGain(j), ...
        G.registration(k).rmse_mv, L.rmseMv(j), ...
        'VariableNames', {'sensor_id','r5_tau_mm','legacy_tau_mm','delta_tau_mm', ...
        'r5_x_scale','legacy_x_scale','delta_x_scale','r5_voltage_gain', ...
        'legacy_voltage_gain','delta_voltage_gain','r5_rmse_mv','legacy_rmse_mv'});
end
A = vertcat(rows{:});
if ~exist(fileparts(outputFile), 'dir'), mkdir(fileparts(outputFile)); end
writetable(A, outputFile);
save(strrep(outputFile, '.csv', '.mat'), 'A', 'r5File', 'legacyCsv', '-v7.3');
end
