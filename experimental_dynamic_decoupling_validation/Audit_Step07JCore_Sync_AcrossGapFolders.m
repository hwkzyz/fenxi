%% Audit folder-local Step07J helper packages
% This script checks that the three gap-aware folders keep identical
% +step07jcore implementation files. It is an audit only: Step07J does not
% call a shared root package at runtime.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
gapFolders = {
    '20241106_low_speed_gap_prior_decoupling'
    '20250527_low_speed_gap_prior_decoupling'
    '20251222_low_speed_gap_prior_decoupling'
    };
coreFiles = {
    'load_joint_static_eta_preview.m'
    'fixed_sensor_eta_by_position.m'
    'sensor_eta_vector.m'
    'solve_gradient_displacement_vp_seed.m'
    };

rows = {};
allMatch = true;
for ifile = 1:numel(coreFiles)
    relFile = fullfile('+step07jcore', coreFiles{ifile});
    refPath = fullfile(thisDir, gapFolders{1}, relFile);
    if ~isfile(refPath)
        error('Missing reference helper file: %s', refPath);
    end
    refText = fileread(refPath);
    refHash = sha256_text_local(refText);
    for ifolder = 1:numel(gapFolders)
        path = fullfile(thisDir, gapFolders{ifolder}, relFile);
        if ~isfile(path)
            rows(end+1, :) = {string(gapFolders{ifolder}), string(coreFiles{ifile}), ...
                false, "missing", string(path)}; %#ok<SAGROW>
            allMatch = false;
            continue;
        end
        txt = fileread(path);
        hash = sha256_text_local(txt);
        isMatch = strcmp(hash, refHash);
        rows(end+1, :) = {string(gapFolders{ifolder}), string(coreFiles{ifile}), ...
            isMatch, string(hash), string(path)}; %#ok<SAGROW>
        allMatch = allMatch && isMatch;
    end
end

AuditTable = cell2table(rows, 'VariableNames', ...
    {'folder', 'helper_file', 'matches_reference', 'sha256', 'path'});
disp(AuditTable);

if ~allMatch
    error(['Step07J helper packages are not synchronized. Keep dataset folders ' ...
        'independent, but copy intentional helper changes to all three +step07jcore folders.']);
end

fprintf('\nStep07J helper packages are synchronized across all gap-aware folders.\n');

function hash = sha256_text_local(txt)
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(uint8(txt(:)));
bytes = typecast(md.digest(), 'uint8');
hash = lower(reshape(dec2hex(bytes).', 1, []));
end
