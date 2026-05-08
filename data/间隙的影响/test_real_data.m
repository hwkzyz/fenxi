% Test with real data, very basic
d = fileparts(mfilename('fullpath'));
dat = readmatrix(fullfile(d, 'calibration_data.txt'), 'FileType', 'text', 'CommentStyle', '%');
dat = dat(all(isfinite(dat(:,1:2)),2), 1:2);
xAll = dat(:,1); yAll = dat(:,2);
bIdx = [find(diff(xAll)<0); numel(xAll)];
sIdx = [1; bIdx(1:end-1)+1];
NG = numel(bIdx);
GL = (0.2:0.2:1.4)';
fprintf('NG = %d\n', NG);

xC = cell(NG,1);
yC = cell(NG,1);
for i = 1:NG
    idx = sIdx(i):bIdx(i);
    xC{i} = xAll(idx);
    yC{i} = yAll(idx);
end

qL = (0.10:0.05:0.90)';
Vt = 3.0e5;

% Library
xAllCat = cell2mat(xC);
xF = linspace(min(xAllCat)-mean(xAllCat), max(xAllCat)-mean(xAllCat), 2000)';
fprintf('xF: %d points, [%.4f, %.4f]\n', numel(xF), min(xF), max(xF));

YL = zeros(NG, numel(xF));
for ig = 1:NG
    xg = xC{ig}(:); yg = yC{ig}(:);
    base = min(yg); yBc = yg - base;
    [~, iPk] = max(yBc);
    xc = xg - xg(iPk); yn = yBc / max(yBc);
    tmp = interp1(xc, yn, xF, 'linear', 'extrap');
    tmp = max(0, min(1, tmp));
    YL(ig, :) = tmp;
    fprintf('  lib %d: [%.4f, %.4f], nFin=%d\n', ig, min(tmp), max(tmp), sum(isfinite(tmp)));
end

% Test case
itest = 4;
xT = xC{itest}(:); yT = yC{itest}(:);
Amp = 0.5; Frq = 600;

x0 = xT; x0c = mean(xT);
for iter = 1:4
    t = (x0 - x0c) / Vt;
    u = Amp * sin(2*pi*Frq*t);
    x0 = linspace(min(xT)+max(u), max(xT)+min(u), numel(xT))';
end
t = (x0 - x0c) / Vt; u = Amp * sin(2*pi*Frq*t);
xi = x0 - u;

yObs = interp1(xT, yT, xi, 'linear', 'extrap');
bObs = min(yObs); yBcO = yObs - bObs;
[~, iPkO] = max(yBcO);
xcO = x0 - x0(iPkO); ynO = yBcO / max(yBcO);
ynO = max(0, min(1, ynO));
yOG = interp1(xcO, ynO, xF, 'linear', 'extrap');
yOG = max(0, min(1, yOG));

fprintf('yOG: class=%s, size=[%d,%d], nFin=%d\n', class(yOG), size(yOG,1), size(yOG,2), sum(isfinite(yOG)));
fprintf('YL: class=%s, size=[%d,%d]\n', class(YL), size(YL,1), size(YL,2));

fprintf('About to loop NG=%d...\n', NG);
for ii = 1:NG
    fprintf('  ii=%d, NG=%d...', ii, NG);
    valid = isfinite(yOG) & isfinite(YL(ii, :));
    nv = sum(valid);
    fprintf(' nv=%d\n', nv);
    if ii > 10
        fprintf('BREAKING: ii exceeded 10!\n');
        break;
    end
end

disp('DONE');
