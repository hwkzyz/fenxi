function manifest = R5_Build_ExternalInputManifest(caseRoot, outputFile)
%R5_BUILD_EXTERNALINPUTMANIFEST Index external gap_only artifacts without copying data.
% The manifest is deliberately read-only: it records candidate files so the
% case adapter can be reviewed before any dynamic identification is run.
if nargin < 1 || isempty(caseRoot), error('R5:CaseRootRequired','caseRoot is required.'); end
if nargin < 2 || isempty(outputFile), outputFile=fullfile(caseRoot,'external_input_manifest.mat'); end
assert(isfolder(caseRoot),'R5:MissingCaseRoot','Missing case root: %s',caseRoot);
patterns={'*Response*Surface*.mat','*Template*.mat','*GapLibrary*.mat','*Foundation*.mat','*GapAware*.mat','*sidecar*.mat'};
items=struct('role',{},'file',{},'bytes',{});
for k=1:numel(patterns)
    d=dir(fullfile(caseRoot,'**',patterns{k}));
    for j=1:numel(d)
        if d(j).isdir, continue; end
        role=classify_role(d(j).name);
        items(end+1)=struct('role',role,'file',fullfile(d(j).folder,d(j).name),'bytes',d(j).bytes); %#ok<AGROW>
    end
end
[~,ia]=unique({items.file},'stable'); items=items(ia);
manifest=struct('caseRoot',caseRoot,'created',datestr(now,30),'items',items);
save(outputFile,'manifest','-v7.3');
end
function role=classify_role(name)
n=lower(name);
if contains(n,'response')||contains(n,'surface'), role='response_surface';
elseif contains(n,'template'), role='low_speed_template';
elseif contains(n,'gaplibrary'), role='gap_library';
elseif contains(n,'foundation'), role='foundation';
elseif contains(n,'gapaware'), role='v1_gap_result';
elseif contains(n,'sidecar'), role='sensor_sidecar';
else, role='unclassified'; end
end
