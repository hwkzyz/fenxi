function assert_latest_formal_path(root)
%ASSERT_LATEST_FORMAL_PATH Prevent archived implementations from shadowing the formal backend.
if nargin < 1, root = fileparts(fileparts(mfilename('fullpath'))); end
names = {'run_unified_blind_dynamic_backend', ...
    'estimate_gap_init_vib_basis_projected', ...
    'Run_R2_UnifiedBlindBackend', 'Run_R3_UnifiedBlindBackend', ...
    'Main_10_UnifiedBlindBackendValidation'};
for i = 1:numel(names)
    p = which(names{i});
    if isempty(p) || ~startsWith(string(p), string(root), 'IgnoreCase', true)
        error('formal:pathShadowed', ...
            'Formal function %s is not resolved from latest program: %s', names{i}, p);
    end
end
end
