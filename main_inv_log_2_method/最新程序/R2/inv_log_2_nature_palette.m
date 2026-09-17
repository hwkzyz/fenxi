function colors = inv_log_2_nature_palette(keys)
%inv_log_2_nature_palette  Soft, color-blind-aware paper figure colors.

P = struct();
P.proposed = hex2rgb_local('#5AB4AC');
P.proposedDark = hex2rgb_local('#2C7F79');
P.blue = hex2rgb_local('#92C5DE');
P.lightBlue = hex2rgb_local('#BBDDEB');
P.wrong = hex2rgb_local('#F4A582');
P.orange = hex2rgb_local('#FDD9A0');
P.gray = hex2rgb_local('#B8B8B8');
P.darkGray = hex2rgb_local('#4D4D4D');
P.lightGray = hex2rgb_local('#E2E2E2');
P.purple = hex2rgb_local('#B3A2C7');
P.green = hex2rgb_local('#A6DBA0');
P.risk = hex2rgb_local('#D8A1A1');
P.band = hex2rgb_local('#D9F0F0');
P.axis = hex2rgb_local('#333333');

if nargin < 1 || isempty(keys)
    colors = P;
    return;
end

keys = string(keys);
colors = zeros(numel(keys), 3);
for ik = 1:numel(keys)
    key = lower(char(keys(ik)));
    if contains(key, 'wrong') || contains(key, 'random') || contains(key, 'naive') || ...
            any(strcmp(key, {'w', 'wfg', 'dr', 'dn'}))
        colors(ik, :) = P.wrong;
    elseif contains(key, 'proposed') || contains(key, 'joint') || contains(key, 'adaptive') || ...
            contains(key, 'proj') || contains(key, 'p2.5') || any(strcmp(key, {'j', 'jm'}))
        colors(ik, :) = P.proposed;
    elseif contains(key, 'vp init') || contains(key, 'p2.0') || strcmp(key, 'dvp')
        colors(ik, :) = P.lightBlue;
    elseif contains(key, 'fixed') || contains(key, 'reference') || startsWith(key, 'f') || strcmp(key, 't')
        colors(ik, :) = P.gray;
    elseif contains(key, 'direct')
        colors(ik, :) = P.purple;
    else
        colors(ik, :) = P.blue;
    end
end
end

function rgb = hex2rgb_local(hex)
hex = char(hex);
if hex(1) == '#'
    hex = hex(2:end);
end
rgb = [hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))] / 255;
end
