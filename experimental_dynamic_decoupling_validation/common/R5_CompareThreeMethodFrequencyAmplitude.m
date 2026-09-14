function T = R5_CompareThreeMethodFrequencyAmplitude(foundationFile,v1File,r5File,outputDir)
%R5_COMPARETHREEMETHODFREQUENCYAMPLITUDE Compare Foundation, legacy V1 and R5.
% Only frequency and amplitude are compared here. All files must describe
% the same windows; this routine never selects or tunes a model.
if nargin<4 || isempty(outputDir), outputDir=fileparts(r5File); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
F=readRows(foundationFile,'foundation'); V=readRows(v1File,'v1'); R=readRows(r5File,'r5');
ids=intersect(intersect(F.window_id,V.window_id),R.window_id); ids=ids(:);
assert(~isempty(ids),'R5:NoCommonWindows');
T=table(ids,nan(size(ids)),nan(size(ids)),nan(size(ids)),nan(size(ids)),nan(size(ids)),nan(size(ids)), ...
 'VariableNames',{'window_id','foundation_frequency_hz','foundation_amplitude_mm','v1_frequency_hz','v1_amplitude_mm','r5_frequency_hz','r5_amplitude_mm'});
for k=1:height(T)
    id=ids(k); T{k,2:end}=[pick(F,id),pick(V,id),pick(R,id)];
end
T.v1_minus_foundation_frequency_hz=T.v1_frequency_hz-T.foundation_frequency_hz;
T.r5_minus_v1_frequency_hz=T.r5_frequency_hz-T.v1_frequency_hz;
T.r5_minus_foundation_frequency_hz=T.r5_frequency_hz-T.foundation_frequency_hz;
T.v1_minus_foundation_amplitude_mm=T.v1_amplitude_mm-T.foundation_amplitude_mm;
T.r5_minus_v1_amplitude_mm=T.r5_amplitude_mm-T.v1_amplitude_mm;
T.r5_minus_foundation_amplitude_mm=T.r5_amplitude_mm-T.foundation_amplitude_mm;
writetable(T,fullfile(outputDir,'ThreeMethod_FrequencyAmplitude_Comparison.csv'));
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 10.5]);
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax=nexttile; plot(ax,T.window_id,T.foundation_frequency_hz,'ko-',T.window_id,T.v1_frequency_hz,'b.-',T.window_id,T.r5_frequency_hz,'r.-','LineWidth',1); ylabel(ax,'Frequency (Hz)'); legend(ax,'Foundation','V1 gap-aware','R5','Location','best'); title(ax,'Frequency comparison'); styleAxes(ax);
ax=nexttile; plot(ax,T.window_id,T.foundation_amplitude_mm,'ko-',T.window_id,T.v1_amplitude_mm,'b.-',T.window_id,T.r5_amplitude_mm,'r.-','LineWidth',1); ylabel(ax,'Amplitude (mm)'); xlabel(ax,'Window'); legend(ax,'Foundation','V1 gap-aware','R5','Location','best'); title(ax,'Amplitude comparison'); styleAxes(ax);
exportgraphics(fig,fullfile(outputDir,'ThreeMethod_FrequencyAmplitude_Comparison.png'),'Resolution',300); close(fig);
end

function styleAxes(ax)
set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',0.75,'Box','on','TickDir','in','XGrid','off','YGrid','off');
end

function X=readRows(file,tag)
S=load(file); if isfield(S,'Result'), Q=S.Result; elseif isfield(S,'out'), Q=S.out; else, error('R5:UnknownResultSchema'); end
if isfield(Q,'WindowResult'), W=Q.WindowResult; elseif isfield(Q,'rows'), W=Q.rows; else, error('R5:NoWindowRows'); end
n=numel(W); X=table(nan(n,1),nan(n,1),nan(n,1),'VariableNames',{'window_id','frequency_hz','amplitude_mm'});
for i=1:n
    w=W(i); X.window_id(i)=num(w,{'window_id','windowId'},i);
    if isfield(w,'Result'), z=w.Result; elseif isfield(w,'CoreResult'), z=w.CoreResult; else, z=w; end
    % Legacy V1 stores the gap-aware estimate under modelFits.gap_only.
    % Keep this read-only extraction independent from the R5 result schema.
    if isfield(w,'modelFits') && isfield(w.modelFits,'gap_only') && any(strcmp(tag,{'v1','r5'}))
        z=w.modelFits.gap_only;
    end
    X.frequency_hz(i)=num(z,{'frequency_hz','fn_id','frequencyHz','freqHz','f_hz'},NaN);
    X.amplitude_mm(i)=num(z,{'amplitude_mm','A_id','amplitudeMm'},NaN);
    if isnan(X.frequency_hz(i)), X.frequency_hz(i)=num(w,{'frequency_hz','fn_id','freqHz','frequencyHz'},NaN); end
    if isnan(X.amplitude_mm(i)), X.amplitude_mm(i)=num(w,{'amplitude_mm','A_id'},NaN); end
end
X.Properties.Description=tag;
end
function y=num(s,names,d)
y=d; for j=1:numel(names), if isfield(s,names{j}) && isnumeric(s.(names{j})) && isscalar(s.(names{j})), y=double(s.(names{j})); return; end, end
end
function v=pick(X,id), q=X.window_id==id; assert(nnz(q)==1,'R5:DuplicateWindow'); v=[X.frequency_hz(q),X.amplitude_mm(q)]; end
