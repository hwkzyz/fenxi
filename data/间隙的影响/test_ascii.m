% ASCII-only test script
d = fileparts(mfilename('fullpath'));
dat = readmatrix(fullfile(d, 'calibration_data.txt'), 'FileType', 'text', 'CommentStyle', '%');
dat = dat(all(isfinite(dat(:,1:2)),2), 1:2);
xAll = dat(:,1); yAll = dat(:,2);
bIdx = [find(diff(xAll)<0); numel(xAll)];
sIdx = [1; bIdx(1:end-1)+1];
nGap = numel(bIdx);
gapList = (0.2:0.2:1.4)';
fprintf('nGap = %d\n', nGap);

xCell = cell(nGap,1);
yCell = cell(nGap,1);
for i = 1:nGap
    idx = sIdx(i):bIdx(i);
    xCell{i} = xAll(idx);
    yCell{i} = yAll(idx);
end

qLevels = (0.10:0.05:0.90)';
nQ = numel(qLevels);
V_tip = 3.0e5;

% Library
xAllCat = cell2mat(xCell);
xcFine = linspace(min(xAllCat)-mean(xAllCat), max(xAllCat)-mean(xAllCat), 2000)';

Ylib = zeros(nGap, numel(xcFine));
for ig = 1:nGap
    xg = xCell{ig}(:); yg = yCell{ig}(:);
    base = min(yg); yBc = yg - base;
    [~, iPeak] = max(yBc);
    xc = xg - xg(iPeak); yn = yBc / max(yBc);
    tmp = interp1(xc, yn, xcFine, 'linear', 'extrap');
    tmp = max(0, min(1, tmp));
    Ylib(ig, :) = tmp;
end

% Single test case
igTest = 4; gTrue = gapList(igTest);
xTrue = xCell{igTest}(:); yTrue = yCell{igTest}(:);
A = 0.5; f = 600;

x0 = xTrue; x0Center = mean(xTrue);
for iter = 1:4
    t = (x0 - x0Center) / V_tip;
    u = A * sin(2*pi*f*t);
    x0 = linspace(min(xTrue)+max(u), max(xTrue)+min(u), numel(xTrue))';
end
t = (x0 - x0Center) / V_tip; u = A * sin(2*pi*f*t);
xi = x0 - u;

if any(diff(xi) <= 0) || min(xi) < min(xTrue)-1e-8 || max(xi) > max(xTrue)+1e-8
    error('Invalid xi');
end

yObs = interp1(xTrue, yTrue, xi, 'linear', 'extrap');
baseObs = min(yObs); yObsBc = yObs - baseObs;
[~, iPeakObs] = max(yObsBc);
xcObs = x0 - x0(iPeakObs); ynObs = yObsBc / max(yObsBc);
ynObs = max(0, min(1, ynObs));
yObsOnGrid = interp1(xcObs, ynObs, xcFine, 'linear', 'extrap');
yObsOnGrid = max(0, min(1, yObsOnGrid));

fprintf('yObsOnGrid size: %d x %d, nFinite: %d\n', ...
    size(yObsOnGrid,1), size(yObsOnGrid,2), sum(isfinite(yObsOnGrid)));

% Try loop
fprintf('About to loop nGap=%d iterations...\n', nGap);
for jj = 1:nGap
    valid = isfinite(yObsOnGrid) & isfinite(Ylib(jj, :));
    nv = sum(valid);
    fprintf('jj=%d, nValid=%d\n', jj, nv);
    if nv > 50
        aa = yObsOnGrid(valid);
        bb = Ylib(jj, valid);
        fprintf('  aa: [%.4f, %.4f] std=%.4f\n', min(aa), max(aa), std(aa));
        fprintf('  bb: [%.4f, %.4f] std=%.4f\n', min(bb), max(bb), std(bb));
        if std(aa) > 0 && std(bb) > 0
            am = aa - mean(aa); bm = bb - mean(bb);
            cc = (am' * bm) / (sqrt(am'*am)*sqrt(bm'*bm));
            fprintf('  corr = %.6f\n', cc);
        end
    end
end

disp('SUCCESS');
