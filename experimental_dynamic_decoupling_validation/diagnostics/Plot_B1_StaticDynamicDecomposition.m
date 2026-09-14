function report = Plot_B1_StaticDynamicDecomposition(resultFile,windowIndex,outputDir)
%PLOT_B1_STATICDYNAMICDECOMPOSITION Show static-to-dynamic voltage evidence.
% This is read-only: it uses the stored bundle and fits and never re-fits.
% Each sensor/turn is plotted as a separate segment, so no cross-turn line
% can create an artificial waveform.
if nargin<2 || isempty(windowIndex), windowIndex=1; end
if nargin<3 || isempty(outputDir), outputDir=fullfile(fileparts(resultFile),'static_dynamic_decomposition'); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
S=load(resultFile); R=pickResult(S); W=R.WindowResult; assert(windowIndex<=numel(W),'Window index out of range.');
w=W(windowIndex); assert(isfield(w,'bundle'),'Result window has no bundle.'); b=w.bundle;
ids=double(b.sensorIds(:).'); x=double(b.X(:)); y=double(b.V(:)); sid=double(b.sensorIndex(:)); f=pickFit(w);
full=vecField(f,'VPred',numel(y)); static=vecField(b,'F0',numel(y));
fixed=[]; if isfield(w,'modelFits') && isfield(w.modelFits,'fixed') && isstruct(w.modelFits.fixed)
    fixed=vecField(w.modelFits.fixed,'VPred',numel(y));
end
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 11]);
tl=tiledlayout(fig,numel(ids),1,'TileSpacing','compact','Padding','compact');
rows=struct('sensorId',{},'n',{},'rmseStaticMv',{},'rmseFullMv',{},'rmseNoGapMv',{});
for k=1:numel(ids)
    q=sid==k; nexttile(tl); hold on; box on; styleAxes(gca);
    plotSegments(x(q),y(q),'Color','k','LineStyle','-','DisplayName','observed');
    plotSegments(x(q),static(q),'Color',[.45 .45 .45],'LineStyle','--','DisplayName','static F_0');
    if ~isempty(fixed), plotSegments(x(q),fixed(q),'Color','b','LineStyle','-','DisplayName','no-gap dynamic'); end
    if ~isempty(full), plotSegments(x(q),full(q),'Color','r','LineStyle','-','DisplayName','gap-aware full-wave'); end
    ylabel(sprintf('S%d (mV)',ids(k))); if k==1, title(sprintf('B1 static to dynamic decomposition | window %d',windowIndex),'FontWeight','normal'); end
    if k==numel(ids), xlabel('x (mm)'); legend('Location','best','FontSize',7); end
    rows(k).sensorId=ids(k); rows(k).n=nnz(q);
    rows(k).rmseStaticMv=rmse(y(q),static(q)); rows(k).rmseFullMv=rmse(y(q),pickSlice(full,q)); rows(k).rmseNoGapMv=rmse(y(q),pickSlice(fixed,q));
end
exportgraphics(fig,fullfile(outputDir,sprintf('B1_static_dynamic_decomposition_W%02d.png',windowIndex)),'Resolution',300); close(fig);
T=struct2table(rows); writetable(T,fullfile(outputDir,sprintf('B1_static_dynamic_decomposition_W%02d.csv',windowIndex)));
report=struct('schema','B1_STATIC_DYNAMIC_DECOMPOSITION_V1','resultFile',resultFile,'windowIndex',windowIndex,'outputDir',outputDir,'hasStatic',~isempty(static),'hasNoGap',~isempty(fixed),'hasFull',~isempty(full));
save(fullfile(outputDir,sprintf('B1_static_dynamic_decomposition_W%02d.mat',windowIndex)),'report','T','-v7.3');
end
function R=pickResult(S), if isfield(S,'Result'), R=S.Result; elseif isfield(S,'out'), R=S.out; else, error('Unknown result schema.'); end, end
function f=pickFit(w), f=struct; if isfield(w,'modelFits') && isfield(w.modelFits,'gap_only') && isstruct(w.modelFits.gap_only), f=w.modelFits.gap_only; end, end
function z=vecField(s,nm,n), z=[]; if isfield(s,nm) && isnumeric(s.(nm)) && numel(s.(nm))==n, z=double(s.(nm)(:)); end, end
function e=rmse(a,b), q=isfinite(a)&isfinite(b); if any(q), e=sqrt(mean((a(q)-b(q)).^2)); else, e=NaN; end, end
function z=pickSlice(a,q), if isempty(a), z=nan(nnz(q),1); else, z=a(q); end, end
function styleAxes(ax), set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',.75,'TickDir','in','Box','on','XGrid','off','YGrid','off'); end
function plotSegments(x,y,varargin)
x=double(x(:)); y=double(y(:)); q=isfinite(x)&isfinite(y); x=x(q); y=y(q); if isempty(x), return; end
span=max(x)-min(x); d=diff(x); cuts=[0;find(d<-.25*span | abs(d)>max(.25,.35*span));numel(x)]; first=true;
for j=1:numel(cuts)-1, ii=(cuts(j)+1):cuts(j+1); if numel(ii)<2, continue; end; [~,o]=sort(x(ii)); ii=ii(o); a=varargin; if ~first, a(end+1:end+2)={'HandleVisibility','off'}; end; plot(x(ii),y(ii),a{:}); first=false; end
end
