% Ultra-simple loop test with real data
d = fileparts(mfilename('fullpath'));
dat = readmatrix(fullfile(d, 'calibration_data.txt'), 'FileType', 'text', 'CommentStyle', '%');
dat = dat(all(isfinite(dat(:,1:2)),2), 1:2);
xAll = dat(:,1); yAll = dat(:,2);
bIdx = [find(diff(xAll)<0); numel(xAll)];
sIdx = [1; bIdx(1:end-1)+1];
NG = numel(bIdx);
GL = (0.2:0.2:1.4)';
fprintf('NG = %d, GL = %s\n', NG, mat2str(GL'));

xC = cell(NG,1); yC = cell(NG,1);
for i = 1:NG
    idx = sIdx(i):bIdx(i);
    xC{i} = xAll(idx); yC{i} = yAll(idx);
end

xAllCat = cell2mat(xC);
xF = linspace(min(xAllCat)-mean(xAllCat), max(xAllCat)-mean(xAllCat), 2000)';
fprintf('xF: %d points\n', numel(xF));

YL = zeros(NG, numel(xF));
for ig = 1:NG
    xg = xC{ig}(:); yg = yC{ig}(:);
    base = min(yg); yBc = yg - base;
    [~, iPk] = max(yBc);
    xc = xg - xg(iPk); yn = yBc / max(yBc);
    tmp = interp1(xc, yn, xF, 'linear', 'extrap');
    tmp = max(0, min(1, tmp));
    YL(ig, :) = tmp;
end

% Test case
itest = 4; Vt = 3.0e5;
xT = xC{itest}(:); yT = yC{itest}(:);
x0 = xT; x0c = mean(xT);
for iter = 1:4
    t = (x0 - x0c) / Vt;
    u = 0.5 * sin(2*pi*600*t);
    x0 = linspace(min(xT)+max(u), max(xT)+min(u), numel(xT))';
end
t = (x0 - x0c) / Vt; u = 0.5 * sin(2*pi*600*t);
xi = x0 - u;
yObs = interp1(xT, yT, xi, 'linear', 'extrap');
bObs = min(yObs); yBcO = yObs - bObs;
[~, iPkO] = max(yBcO);
xcO = x0 - x0(iPkO); ynO = yBcO / max(yBcO);
ynO = max(0, min(1, ynO));
yOG = interp1(xcO, ynO, xF, 'linear', 'extrap');
yOG = max(0, min(1, yOG));

fprintf('yOG size: [%d %d]\n', size(yOG,1), size(yOG,2));
fprintf('YL size: [%d %d]\n', size(YL,1), size(YL,2));

fprintf('Entering loop, NG=%d...\n', NG);
for ii = 1:NG
    fprintf('ii=%d ', ii);
    % Just print dimensions
    row = YL(ii, :);
    fprintf('row_size=[%d %d] ', size(row,1), size(row,2));
    valid = isfinite(yOG) & isfinite(row);
    fprintf('valid_size=[%d %d] ', size(valid,1), size(valid,2));
    nv = sum(valid);
    fprintf('nv_size=[%d %d]\n', size(nv,1), size(nv,2));
end

disp('LOOP DONE');
