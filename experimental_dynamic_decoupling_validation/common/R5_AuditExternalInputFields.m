function report = R5_AuditExternalInputFields(files, outputFile)
%R5_AUDITEXTERNALINPUTFIELDS Record MAT top-level fields before adaptation.
if ischar(files) || isstring(files), files=cellstr(files); end
assert(iscell(files),'R5:FilesRequired','files must be a cell array.');
report=struct('file',{},'exists',{},'variables',{},'bytes',{});
for k=1:numel(files)
    f=char(files{k}); e=isfile(f); report(k).file=f; report(k).exists=e;
    report(k).variables={}; report(k).bytes=0;
    if e
        w=whos('-file',f); report(k).variables={w.name}; report(k).bytes=sum([w.bytes]);
    end
end
if nargin>=2 && ~isempty(outputFile), save(outputFile,'report','-v7.3'); end
end
