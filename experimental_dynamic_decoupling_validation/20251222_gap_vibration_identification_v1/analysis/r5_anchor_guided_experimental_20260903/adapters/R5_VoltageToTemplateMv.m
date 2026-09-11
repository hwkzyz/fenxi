function vMv = R5_VoltageToTemplateMv(vRaw, templateRow)
%R5_VOLTAGETOTEMPLATEMV Frozen Main10 voltage contract for R5 adapters.
% Raw Main07/Main08 voltage is converted to the Main05 mV baseline convention.
vRaw = double(vRaw(:));
assert(isfield(templateRow,'baseline') && isfinite(templateRow.baseline), ...
    'R5:VoltageContract','Template baseline is required for mV conversion.');
vMv = (vRaw - double(templateRow.baseline)) * 1000;
end
