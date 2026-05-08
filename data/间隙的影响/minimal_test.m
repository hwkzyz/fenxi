% Minimal test
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');

data = readmatrix(dataFile, 'FileType', 'text', 'CommentStyle', '%');
data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
xAll = data(:, 1);
yAll = data(:, 2);
breakIdx = [find(diff(xAll) < 0); numel(xAll)];
startIdx = [1; breakIdx(1:end-1) + 1];
nGap = numel(breakIdx);
fprintf('nGap = %d\n', nGap);

gapList = (0.2:0.2:1.4)';
xCell = cell(nGap, 1);
yCell = cell(nGap, 1);
for i = 1:nGap
    idx = startIdx(i):breakIdx(i);
    xCell{i} = xAll(idx);
    yCell{i} = yAll(idx);
    fprintf('  curve %d: %d points\n', i, numel(xCell{i}));
end

qLevels = (0.10:0.05:0.90)';
V_tip = 3.0e5;

% Build library
xAllCat = cell2mat(xCell);
xcFine = linspace(min(xAllCat)-mean(xAllCat), max(xAllCat)-mean(xAllCat), 2000)';
fprintf('xcFine: %d points, range [%.4f, %.4f]\n', numel(xcFine), min(xcFine), max(xcFine));

Y_lib_raw = zeros(nGap, numel(xcFine));
for ig = 1:nGap
    xg = xCell{ig}(:);
    yg = yCell{ig}(:);
    base = min(yg);
    yBc = yg - base;
    [~, iPeak] = max(yBc);
    xc = xg - xg(iPeak);
    yn = yBc / max(yBc);
    tmp = interp1(xc, yn, xcFine, 'linear', 'extrap');
    tmp = max(0, min(1, tmp));
    Y_lib_raw(ig, :) = tmp;
    fprintf('  lib %d: range [%.4f, %.4f], nFinite=%d\n', ig, min(tmp), max(tmp), sum(isfinite(tmp)));
end

% Test one case
ig = 4;
gTrue = gapList(ig);
xTrue = xCell{ig}(:);
yTrue = yCell{ig}(:);
A = 0.5; f = 600;

x0 = xTrue;
x0Center = mean(xTrue);
for iter = 1:4
    t = (x0 - x0Center) / V_tip;
    u = A * sin(2*pi*f*t);
    x0 = linspace(min(xTrue)+max(u), max(xTrue)+min(u), numel(xTrue))';
end
t = (x0 - x0Center) / V_tip;
u = A * sin(2*pi*f*t);
xi = x0 - u;

% Check xi
fprintf('xi: diff_min=%.6e, range=[%.4f,%.4f] vs xTrue=[%.4f,%.4f]\n', ...
    min(diff(xi)), min(xi), max(xi), min(xTrue), max(xTrue));

yObs = interp1(xTrue, yTrue, xi, 'linear', 'extrap');
fprintf('yObs: range [%.4f, %.4f], nFinite=%d\n', min(yObs), max(yObs), sum(isfinite(yObs)));

baseObs = min(yObs);
yObsBc = yObs - baseObs;
[~, iPeakObs] = max(yObsBc);
xcObs = x0 - x0(iPeakObs);
ynObs = yObsBc / max(yObsBc);
ynObs = max(0, min(1, ynObs));
fprintf('ynObs: range [%.4f, %.4f]\n', min(ynObs), max(ynObs));

yObsOnGrid = interp1(xcObs, ynObs, xcFine, 'linear', 'extrap');
yObsOnGrid = max(0, min(1, yObsOnGrid));
fprintf('yObsOnGrid: range [%.4f, %.4f], nFinite=%d\n', ...
    min(yObsOnGrid), max(yObsOnGrid), sum(isfinite(yObsOnGrid)));

% Try the correlation
for j = 1:nGap
    valid = isfinite(yObsOnGrid) & isfinite(Y_lib_raw(j, :));
    nv = sum(valid);
    fprintf('j=%d, nValid=%d\n', j, nv);
    if nv > 50
        a = yObsOnGrid(valid);
        b = Y_lib_raw(j, valid);
        fprintf('  a: [%.6f,%.6f] std=%.6f\n', min(a), max(a), std(a));
        fprintf('  b: [%.6f,%.6f] std=%.6f\n', min(b), max(b), std(b));
        if std(a) > 0 && std(b) > 0
            am = a - mean(a);
            bm = b - mean(b);
            cnt = (am' * bm) / (sqrt(am'*am)*sqrt(bm'*bm));
            fprintf('  corr = %.6f\n', cnt);
        end
    end
end

disp('ALL DONE');
