function [gapList, xCell, yCell] = load_stacked_gap_curves(filePath, gapListInput)
%load_stacked_gap_curves  读取静态间隙波形库。
%
% 支持两种 COMSOL 文本导出排列方式：
%   1) 曲线堆叠格式：
%      x1 y(g1), x2 y(g1), ..., xN y(g1),
%      x1 y(g2), x2 y(g2), ..., xN y(g2), ...
%   2) 空间点交错格式：
%      x1 y(g1), x1 y(g2), ..., x1 y(gM),
%      x2 y(g1), x2 y(g2), ..., x2 y(gM), ...

data = readmatrix(filePath, 'FileType', 'text', 'CommentStyle', '%');
data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
x = data(:, 1);
y = data(:, 2);

if nargin < 2
    gapListInput = [];
end
gapListHeader = read_gap_list_from_header(filePath);

if any(diff(x) < 0)
    [xCell, yCell] = split_curve_stacked_layout(x, y);
else
    [xCell, yCell] = split_point_interleaved_layout(x, y);
end

nCurve = numel(xCell);
if ~isempty(gapListInput)
    gapList = gapListInput(:);
    if numel(gapList) ~= nCurve
        error('Specified gapList has %d values, but data contains %d curves.', ...
            numel(gapList), nCurve);
    end
elseif ~isempty(gapListHeader)
    gapList = gapListHeader(:);
    if numel(gapList) ~= nCurve
        error('Header gapList_mm has %d values, but data contains %d curves.', ...
            numel(gapList), nCurve);
    end
elseif nCurve == 7
    gapList = (0.2:0.2:1.4)';
elseif nCurve == 6 && contains(filePath, '0.3_1.5')
    warning(['检测到 0.3_1.5 文件中有 6 条曲线，但文本头没有保存中间间隙值。' ...
        '当前返回曲线序号。请用 load_stacked_gap_curves(filePath, gapList) 显式传入真实间隙序列。']);
    gapList = (1:nCurve)';
else
    gapList = (1:nCurve)';
end
end

function gapList = read_gap_list_from_header(filePath)
gapList = [];
fid = fopen(filePath, 'r');
if fid < 0
    return;
end
cleanupObj = onCleanup(@() fclose(fid));

while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    line = strtrim(line);
    if startsWith(line, '%')
        if contains(line, 'gapList_mm')
            token = regexp(line, 'gapList_mm\s*:\s*(.*)$', 'tokens', 'once');
            if ~isempty(token)
                gapList = sscanf(token{1}, '%f').';
            end
            return;
        end
    elseif ~isempty(line)
        return;
    end
end
end

function [xCell, yCell] = split_curve_stacked_layout(x, y)
breakIdx = [find(diff(x) < 0); numel(x)];
startIdx = [1; breakIdx(1:end-1) + 1];
nCurve = numel(breakIdx);
xCell = cell(nCurve, 1);
yCell = cell(nCurve, 1);
for i = 1:nCurve
    idx = startIdx(i):breakIdx(i);
    [xUnique, ia] = unique(x(idx), 'stable');
    yUnique = y(idx);
    xCell{i} = xUnique(:);
    yCell{i} = yUnique(ia);
end
end

function [xCell, yCell] = split_point_interleaved_layout(x, y)
tol = max(1, max(abs(x))) * 1e-12;
firstX = x(1);
nCurve = find(abs(x - firstX) > tol, 1, 'first') - 1;
if isempty(nCurve)
    error('Could not split point-interleaved data: all x values are identical.');
end
if mod(numel(x), nCurve) ~= 0
    error('Point-interleaved data length is not divisible by the number of curves.');
end

nPoint = numel(x) / nCurve;
xMat = reshape(x, nCurve, nPoint).';
yMat = reshape(y, nCurve, nPoint).';

xRef = xMat(:, 1);
if any(max(abs(xMat - xRef), [], 2) > tol)
    error('Point-interleaved data is not grouped by identical x positions.');
end

xCell = cell(nCurve, 1);
yCell = cell(nCurve, 1);
for i = 1:nCurve
    xCell{i} = xRef(:);
    yCell{i} = yMat(:, i);
end
end
