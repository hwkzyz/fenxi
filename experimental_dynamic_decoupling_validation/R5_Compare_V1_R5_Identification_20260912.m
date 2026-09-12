function report = R5_Compare_V1_R5_Identification_20260912(caseFiles, r5Files, outputCsv)
% Compare window-wise V1 and R5 results without using strain or locking EO.
% caseFiles/r5Files are cell arrays of equal length; each item is a MAT file
% containing Result.WindowResult. This is a read-only audit/report function.
if nargin<3, outputCsv=''; end
assert(iscell(caseFiles)&&iscell(r5Files)&&numel(caseFiles)==numel(r5Files));
rows=struct('case_id',{},'window_id',{},'v1_EO',{},'r5_EO',{},'v1_f_Hz',{},'r5_f_Hz',{}, ...
    'delta_f_Hz',{},'v1_A_mm',{},'r5_A_mm',{},'delta_A_mm',{},'v1_rmse_mV',{},'r5_rmse_mV',{});
k=0;
for ic=1:numel(caseFiles)
    a=load(caseFiles{ic}); b=load(r5Files{ic});
    if isfield(a,'Result'), r1=a.Result; elseif isfield(a,'out'), r1=a.out; else, error('V1 result variable missing'); end
    if isfield(b,'Result'), r2=b.Result; elseif isfield(b,'out'), r2=b.out; else, error('R5 result variable missing'); end
    w1=r1.WindowResult; w2=r2.rows; n=min(numel(w1),numel(w2));
    for iw=1:n
        p1=pickfit(w1(iw)); p2=w2(iw); k=k+1;
        rows(k).case_id=string(ic); rows(k).window_id=double(getf(w1(iw),'windowId',getf(p2,'window_id',iw)));
        rows(k).v1_EO=getf(p1,'EO',getf(p1,'EO_id',NaN)); rows(k).r5_EO=getf(p2,'EO',getf(p2,'EO_id',NaN));
        rows(k).v1_f_Hz=getf(p1,'fHz',getf(p1,'freqHz',getf(p1,'frequency_hz',getf(p1,'f',NaN)))); rows(k).r5_f_Hz=getf(p2,'fHz',getf(p2,'freqHz',getf(p2,'frequency_hz',getf(p2,'f',NaN))));
        rows(k).delta_f_Hz=rows(k).r5_f_Hz-rows(k).v1_f_Hz;
        rows(k).v1_A_mm=getf(p1,'Amm',getf(p1,'amplitudeMm',getf(p1,'amplitude_mm',getf(p1,'A_mm',NaN)))); rows(k).r5_A_mm=getf(p2,'Amm',getf(p2,'amplitudeMm',getf(p2,'amplitude_mm',getf(p2,'A_mm',NaN))));
        rows(k).delta_A_mm=rows(k).r5_A_mm-rows(k).v1_A_mm;
        rows(k).v1_rmse_mV=getf(p1,'rmseMv',getf(p1,'plainRmseMv',getf(p1,'rmse_mv',getf(p1,'rmse_mV',NaN)))); rows(k).r5_rmse_mV=getf(p2,'rmseMv',getf(p2,'rmse_mv',getf(p2,'rmse_mV',NaN)));
    end
end
report=struct('schema','R5_V1_R5_IDENTIFICATION_COMPARISON_V1','rows',rows);
if ~isempty(outputCsv), writetable(struct2table(rows),outputCsv); end
end
function p=pickfit(w)
p=struct();
if isfield(w,'modelFits') && isfield(w.modelFits,'gap_only'), p=w.modelFits.gap_only;
elseif isfield(w,'fit'), p=w.fit;
elseif isfield(w,'bestFit'), p=w.bestFit; end
if isempty(fieldnames(p)), p=w; end
end
function v=getf(s,n,d)
v=d; if isstruct(s)&&isfield(s,n)&&isnumeric(s.(n))&&~isempty(s.(n))&&isfinite(s.(n)(1)), v=double(s.(n)(1)); end
end
