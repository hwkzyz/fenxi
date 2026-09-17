function Fraw = r3_raw_response(R3,gMm)
%R3_RAW_RESPONSE Return one original COMSOL node curve for truth generation.

idx = find(abs(R3.ctx.gapList-gMm)<1e-12,1);
if isempty(idx)
    error('r3:MissingRawState','Raw FE gap state %.6g mm is unavailable.',gMm);
end
Fraw = griddedInterpolant(R3.ctx.xCell{idx},R3.ctx.yCell{idx},...
    'pchip','nearest');
end
