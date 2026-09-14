function report = Plot_B4_MethodAblationWaveform(resultFile,windowIndex,outputDir,bundleFile)
%PLOT_B4_METHODABLATIONWAVEFORM Compare stored waveform paths in one window.
% The plot is deliberately honest: a method is shown only when its result
% file stores a prediction on the same waveform bundle. Missing predictions
% are reported as unavailable instead of being replaced by another method.
if nargin<2||isempty(windowIndex), windowIndex=1; end
if nargin<3||isempty(outputDir), outputDir=fullfile(fileparts(resultFile),'method_ablation'); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
S=load(resultFile); R=pickResult(S); W=pickWindows(R); assert(windowIndex<=numel(W));
w=W(windowIndex);
if isfield(w,'bundle') || isfield(w,'CoreBundlePreview')
    B=pickBundle(w);
else
    assert(nargin>=4 && ~isempty(bundleFile),'A bundle file is required for row-only results.');
    Q=load(bundleFile); Rb=pickResult(Q); wb=pickWindows(Rb); assert(windowIndex<=numel(wb)); B=pickBundle(wb(windowIndex));
end
if ~isfield(B,'sensorIds')
    ids=unique(double(B.sensor_id(:)),'stable'); B.sensorIds=ids(:).';
    B.X=double(B.x(:)); B.V=1000*double(B.v(:));
    [~,B.sensorIndex]=ismember(double(B.sensor_id(:)),ids);
    if isfield(B,'template_v'), B.F0=1000*double(B.template_v(:)); end
end
fit=pickFit(w); ids=double(B.sensorIds(:).');
y=double(B.V(:)); x=double(B.X(:)); si=double(B.sensorIndex(:));
pred=struct('name',{},'value',{},'available',{});
pred(end+1)=struct('name','observed','value',y,'available',true);
pred(end+1)=struct('name','static baseline','value',getVec(B,'F0',numel(y)),'available',true);
pred(end+1)=struct('name','R5 no-gap','value',getVec(fit,'VPredNoGap',numel(y)),'available',false);
pred(end).available=~isempty(pred(end).value);
pred(end+1)=struct('name','R5 gap-aware','value',getVec(fit,'VPred',numel(y)),'available',false);
pred(end).available=~isempty(pred(end).value);
if isfield(w,'modelFits') && isfield(w.modelFits,'fixed') && isstruct(w.modelFits.fixed)
    z=getVec(w.modelFits.fixed,'VPred',numel(y));
    if ~isempty(z), pred(end+1)=struct('name','Foundation fixed-gap','value',z,'available',true); end
end
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 12]);
tl=tiledlayout(fig,numel(ids),1,'TileSpacing','compact','Padding','compact');
colors=lines(numel(pred)); rows=struct('sensorId',{},'method',{},'available',{},'rmse_mv',{});
for k=1:numel(ids)
    q=si==k; ax=nexttile(tl); hold(ax,'on'); box(ax,'on'); styleAxes(ax);
    for j=1:numel(pred)
        if ~pred(j).available, continue; end
        yy=pred(j).value; if strcmp(pred(j).name,'observed'), c=[0 0 0]; ls='-'; lw=0.8; else, c=colors(j,:); ls='-'; lw=1.0; end
        plotSegments(x(q),yy(q),c,ls,lw,pred(j).name);
        rows(end+1)=struct('sensorId',ids(k),'method',pred(j).name,'available',true,'rmse_mv',rmse(y(q),yy(q))); %#ok<AGROW>
    end
    ylabel(sprintf('S%d (mV)',ids(k))); if k==1, title(sprintf('B4 method ablation | window %d',windowIndex),'FontWeight','normal'); end
    if k==numel(ids), xlabel('x (mm)'); legend(ax,'Location','best','FontSize',7); end
end
exportgraphics(fig,fullfile(outputDir,sprintf('B4_method_ablation_W%02d.png',windowIndex)),'Resolution',300); close(fig);
T=struct2table(rows); writetable(T,fullfile(outputDir,sprintf('B4_method_ablation_W%02d.csv',windowIndex)));
report=struct('schema','B4_METHOD_ABLATION_WAVEFORM_V1','resultFile',resultFile,'windowIndex',windowIndex,'outputDir',outputDir,'methods',{ {pred.name} });
save(fullfile(outputDir,sprintf('B4_method_ablation_W%02d.mat',windowIndex)),'report','T','-v7.3');
end
function R=pickResult(S), if isfield(S,'Result'), R=S.Result; elseif isfield(S,'out'), R=S.out; else, error('Unknown result schema.'); end, end
function W=pickWindows(R), if isfield(R,'WindowResult'), W=R.WindowResult; elseif isfield(R,'rows'), W=R.rows; else, error('No windows.'); end, end
function B=pickBundle(w), if isfield(w,'bundle'), B=w.bundle; elseif isfield(w,'CoreBundlePreview'), B=w.CoreBundlePreview; else, error('No bundle.'); end, end
function f=pickFit(w), f=w; if isfield(w,'modelFits')&&isfield(w.modelFits,'gap_only')&&isstruct(w.modelFits.gap_only), f=w.modelFits.gap_only; end, end
function z=getVec(s,nm,n), z=[]; if isfield(s,nm)&&isnumeric(s.(nm))&&numel(s.(nm))==n, z=double(s.(nm)(:)); end, end
function e=rmse(a,b), q=isfinite(a)&isfinite(b); if any(q), e=sqrt(mean((a(q)-b(q)).^2)); else, e=NaN; end, end
function styleAxes(ax), set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',.75,'TickDir','in','Box','on','XGrid','off','YGrid','off'); end
function plotSegments(x,y,c,ls,lw,name)
x=double(x(:)); y=double(y(:)); q=isfinite(x)&isfinite(y); x=x(q); y=y(q); if isempty(x), return; end
span=max(x)-min(x); d=diff(x); cuts=[0;find(d<-.25*span|abs(d)>max(.25,.35*span));numel(x)]; first=true;
for j=1:numel(cuts)-1, ii=(cuts(j)+1):cuts(j+1); if numel(ii)<2, continue; end; [~,o]=sort(x(ii)); ii=ii(o); if first, plot(x(ii),y(ii),'Color',c,'LineStyle',ls,'LineWidth',lw,'DisplayName',name); first=false; else, plot(x(ii),y(ii),'Color',c,'LineStyle',ls,'LineWidth',lw,'HandleVisibility','off'); end; end
end
