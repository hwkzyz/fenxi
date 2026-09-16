function audit = R5_AuditIndependentVibrationReference(resultFile, referenceCsv, outputDir, opts)
%R5_AUDITINDEPENDENTVIBRATIONREFERENCE Compare R5 vibration estimates with an
% independent strain/BTT reference.  Missing references are reported as
% unavailable; they are never silently treated as a pass.
if nargin < 3 || isempty(outputDir), outputDir = fileparts(resultFile); end
if nargin < 4 || isempty(opts), opts = struct(); end
freqTol = localOpt(opts, 'frequencyToleranceHz', 2.0);
ampRelTol = localOpt(opts, 'amplitudeRelativeTolerance', 0.25);
audit = struct('schema','R5_INDEPENDENT_VIBRATION_AUDIT_V1', ...
    'status','unavailable','pass',false,'resultFile',resultFile, ...
    'referenceFile',referenceCsv,'frequencyToleranceHz',freqTol, ...
    'amplitudeRelativeTolerance',ampRelTol,'matchedCount',0, ...
    'frequencyErrorMeanHz',NaN,'frequencyErrorMaxHz',NaN, ...
    'amplitudeRelativeErrorMean',NaN,'amplitudeRelativeErrorMax',NaN);
if exist(resultFile,'file') ~= 2 || exist(referenceCsv,'file') ~= 2
    localWrite(audit, outputDir); return;
end
S = load(resultFile);
if isfield(S,'out'), rows = S.out.rows; elseif isfield(S,'Result') && isfield(S.Result,'WindowResult')
    rows = S.Result.WindowResult;
else
    localWrite(audit, outputDir); return;
end
T = readtable(referenceCsv);
required = {'WindowID','WaveformFreq_Hz','WaveformAmp_mm','StrainFreq_Hz','StrainDerivedAmp_mm'};
if ~all(ismember(required,T.Properties.VariableNames))
    localWrite(audit, outputDir); return;
end
rid = nan(numel(rows),1); rf = rid; ra = rid;
for k=1:numel(rows)
    rid(k)=localNum(rows(k),{'window_id','WindowID'},k);
    rf(k)=localNum(rows(k),{'frequency_hz','freqHz','WaveformFreq_Hz'});
    ra(k)=localNum(rows(k),{'amplitude_mm','amplitudeMm','WaveformAmp_mm'});
end
[ok,ia,ib] = intersect(rid,double(T.WindowID)); %#ok<ASGLU>
if isempty(ia), localWrite(audit, outputDir); return; end
fRef = double(T.StrainFreq_Hz(ib)); aRef = double(T.StrainDerivedAmp_mm(ib));
valid = isfinite(rf(ia)) & isfinite(ra(ia)) & isfinite(fRef) & isfinite(aRef) & aRef > 0;
ef = rf(ia(valid))-fRef(valid); ea = abs(ra(ia(valid))-aRef(valid))./abs(aRef(valid));
audit.matchedCount = nnz(valid);
if audit.matchedCount > 0
    audit.frequencyErrorMeanHz=mean(abs(ef)); audit.frequencyErrorMaxHz=max(abs(ef));
    audit.amplitudeRelativeErrorMean=mean(ea); audit.amplitudeRelativeErrorMax=max(ea);
    audit.pass = audit.frequencyErrorMaxHz <= freqTol && audit.amplitudeRelativeErrorMax <= ampRelTol;
    audit.status = ternary(audit.pass,'pass','fail');
else
    audit.status='unavailable';
end
localWrite(audit, outputDir);
end

function v=localOpt(s,n,d), v=d; if isfield(s,n)&&isscalar(s.(n))&&isfinite(s.(n)), v=double(s.(n)); end, end
function v=localNum(s,names,d)
if nargin<3,d=NaN;end; v=d;
for i=1:numel(names), if isfield(s,names{i}) && isnumeric(s.(names{i})) && isscalar(s.(names{i})), v=double(s.(names{i})); return; end, end
end
function out=ternary(tf,a,b), if tf,out=a;else,out=b;end,end
function localWrite(audit,dirName)
if ~exist(dirName,'dir'), mkdir(dirName); end
fid=fopen(fullfile(dirName,'IndependentVibrationAudit.md'),'w'); if fid<0,return;end
c=onCleanup(@()fclose(fid)); fprintf(fid,'# Independent vibration reference audit\n\n- Status: **%s**\n- Matched windows: %d\n- Frequency max error: %.6g Hz\n- Amplitude relative max error: %.6g\n',audit.status,audit.matchedCount,audit.frequencyErrorMaxHz,audit.amplitudeRelativeErrorMax);
end
