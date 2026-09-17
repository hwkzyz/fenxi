function Data = Main_01_BuildSimulationSurface()
%MAIN_01_BUILDSIMULATIONSURFACE Assemble raw COMSOL tilt-gap waveforms.
% response_V is indexed as x sample, clearance, and FE tilt angle.

thisDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(thisDir);
addpath(mainDir, '-begin');
addpath(fullfile(mainDir, 'local_func'), '-begin');

gapList = (0.5:0.2:1.5).';
tiltList = [0.5 1 1.5 2 2.5 3 3.5].';
nX = 161;
nGap = numel(gapList);
nTilt = numel(tiltList);
x_mm = nan(nX, 1);
response_V = nan(nX, nGap, nTilt);
sourceFiles = strings(nTilt, 1);

for iTilt = 1:nTilt
    filePath = tilt_file(thisDir, tiltList(iTilt));
    sourceFiles(iTilt) = string(filePath);
    [gRead, xCell, yCell] = load_stacked_gap_curves(filePath, gapList);
    assert(isequal(size(gRead), size(gapList)) && all(abs(gRead-gapList) < 1e-12), ...
        'Unexpected clearance list in %s.', filePath);
    for iGap = 1:nGap
        x = xCell{iGap}(:);
        y = yCell{iGap}(:);
        assert(numel(x) == nX && numel(y) == nX, ...
            'Unexpected sample count for tilt %.1f deg, gap %.1f mm.', ...
            tiltList(iTilt), gapList(iGap));
        if iTilt == 1 && iGap == 1
            x_mm = x;
        else
            assert(max(abs(x-x_mm)) < 1e-12, ...
                'x-grid mismatch at tilt %.1f deg, gap %.1f mm.', ...
                tiltList(iTilt), gapList(iGap));
        end
        response_V(:, iGap, iTilt) = y;
    end
end

outDir = fullfile(thisDir, 'output', 'waveform_families');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
matPath = fullfile(outDir, 'TiltGap_WaveformSurfaceData.mat');
Data = struct();
Data.x_mm = x_mm;
Data.gap_mm = gapList;
Data.tilt_deg = tiltList;
Data.response_V = response_V;
Data.source_files = sourceFiles;
Data.index_definition = 'response_V(x_index, gap_index, tilt_index)';
Data.description = 'Raw COMSOL static sensor-response waveforms; no interpolation or normalization.';
Data.created_by = mfilename;
save(matPath, 'Data', '-v7.3');

fprintf('Saved %d x %d x %d tilt-gap waveform dataset to:\n%s\n', ...
    nX, nGap, nTilt, matPath);
end

function filePath = tilt_file(thisDir, angleDeg)
fileName = sprintf('直叶片2mm_不同间隙0.5_0.2_1.5_倾斜角%g.txt', angleDeg);
filePath = fullfile(thisDir, 'input_fe_tilt', fileName);
assert(exist(filePath, 'file') == 2, 'Missing COMSOL file: %s', filePath);
end
