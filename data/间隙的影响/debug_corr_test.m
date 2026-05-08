% Debug script v2
scriptDir = fileparts(mfilename('fullpath'));
data = readmatrix(fullfile(scriptDir, '直叶片2mm_不同间隙.txt'), 'FileType', 'text', 'CommentStyle', '%');
data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
x = data(:, 1); y = data(:, 2);
breakIdx = [find(diff(x) < 0); numel(x)];
startIdx = [1; breakIdx(1:end-1) + 1];
nCurve = numel(breakIdx);
gapList = (0.2:0.2:1.4)';
xCell = cell(nCurve, 1); yCell = cell(nCurve, 1);
for i = 1:nCurve
    idx = startIdx(i):breakIdx(i);
    xCell{i} = x(idx); yCell{i} = y(idx);
end
nGap = numel(gapList);
V_tip = 3.0e5;

% Build library
xAll = cell2mat(xCell);
xcFine = linspace(min(xAll) - mean(xAll), max(xAll) - mean(xAll), 2000)';
Y_lib_raw = zeros(nGap, numel(xcFine));
for ig = 1:nGap
    xg = xCell{ig}(:); yg = yCell{ig}(:);
    base = min(yg); yBc = yg - base;
    [~, iPeak] = max(yBc);
    xc = xg - xg(iPeak); yn = yBc / max(yBc);
    Y_lib_raw(ig, :) = interp1(xc, yn, xcFine, 'pchip', 'extrap');
    fprintf('Library %d: nFinite=%d/%d\n', ig, sum(isfinite(Y_lib_raw(ig,:))), numel(xcFine));
end

% Test
gTrue = gapList(4); xTrue = xCell{4}(:); yTrue = yCell{4}(:);
x0Center = mean(xTrue); x0 = xTrue;
A = 0.5; f = 600;
for iter = 1:4
    t = (x0 - x0Center) / V_tip;
    u = A * sin(2*pi*f*t);
    x0 = linspace(min(xTrue) + max(u), max(xTrue) + min(u), numel(xTrue))';
end
t = (x0 - x0Center) / V_tip;
u = A * sin(2*pi*f*t);
xi = x0 - u;
yObs = interp1(xTrue, yTrue, xi, 'pchip', 'extrap');
baseObs = min(yObs); yObsBc = yObs - baseObs;
[~, iPeakObs] = max(yObsBc);
xcObs = x0 - x0(iPeakObs); ynObs = yObsBc / max(yObsBc);
yObsOnGrid = interp1(xcObs, ynObs, xcFine, 'pchip', 'extrap');

for j = 1:nGap
    valid = isfinite(yObsOnGrid) & isfinite(Y_lib_raw(j, :));
    nv = sum(valid);
    fprintf('j=%d, nValid=%d\n', j, nv);
    if nv > 50
        a = yObsOnGrid(valid); b = Y_lib_raw(j, valid);
        fprintf('  a: [%.6f, %.6f], std=%.6f\n', min(a), max(a), std(a));
        fprintf('  b: [%.6f, %.6f], std=%.6f\n', min(b), max(b), std(b));
        if std(a) == 0 || std(b) == 0
            fprintf('  SKIP: zero variance\n');
        else
            try
                ct = corrcoef(a(:), b(:));
                fprintf('  corr=%.6f\n', ct(1,2));
            catch ME
                fprintf('  corrcoef ERROR: %s\n', ME.message);
                % manual fallback
                cManual = sum((a-mean(a)).*(b-mean(b))) / (std(a)*std(b)*(numel(a)-1));
                fprintf('  manual corr=%.6f\n', cManual);
            end
        end
    end
end
disp('done');
