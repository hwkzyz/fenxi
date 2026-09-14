function report = Plot_B_DynamicSingleWindowAudit(resultFile,windowIndex,outputDir,bundleFile)
%PLOT_B_DYNAMIC SINGLE-WINDOW audit. Reads stored bundle and fit only.
if nargin<2 || isempty(windowIndex), windowIndex=1; end
if nargin<3 || isempty(outputDir), outputDir=fullfile(fileparts(resultFile),'dynamic_window_audit'); end
if nargin<4, bundleFile=''; end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
S=load(resultFile); R=pick_result(S); W=pick_windows(R); assert(windowIndex<=numel(W),'Window index out of range.');
w=W(windowIndex); fit=pick_fit(w);
if isfield(w,'bundle') || isfield(w,'CoreBundlePreview')
    B=pick_bundle(w);
elseif ~isempty(bundleFile) && isfile(bundleFile)
    Q=load(bundleFile); Q=pick_result(Q); WB=pick_windows(Q); assert(windowIndex<=numel(WB),'Bundle window index out of range.'); B=pick_bundle(WB(windowIndex));
else
    error('Window has no waveform bundle.');
end
report=struct('window',windowIndex,'resultFile',resultFile,'bundleFile',bundleFile);
B = normalize_bundle(B);
fig=figure('Visible','off','Color','w','Position',[80 80 1500 1000]);
ids=double(B.sensorIds(:).'); tl=tiledlayout(fig,numel(ids),1,'TileSpacing','compact','Padding','compact');
for k=1:numel(ids)
    q=double(B.sensorIndex(:))==k; nexttile(tl); hold on; grid on; box on;
    plot_segmented(B.X(q),B.V(q),'k.','DisplayName','observed');
    if isfield(fit,'VPred') && numel(fit.VPred)==numel(B.V), plot_segmented(B.X(q),fit.VPred(q),'r-','DisplayName','full-wave'); end
    title(sprintf('S%d | dynamic waveform fit',ids(k))); xlabel('x (mm)'); ylabel('mV'); legend('Location','best');
end
exportgraphics(fig,fullfile(outputDir,'B1_fullwave_fit.png'),'Resolution',160); close(fig);
if isfield(fit,'uMm') && numel(fit.uMm)==numel(B.X)
    fig=figure('Visible','off','Color','w'); plot(B.T,fit.uMm,'b.'); grid on; xlabel('t (s)'); ylabel('u (mm)'); title(sprintf('Window %d gradient/VP displacement',windowIndex));
    exportgraphics(fig,fullfile(outputDir,'B2_displacement.png'),'Resolution',160); close(fig);
end
if isfield(fit,'CandidateTable') && istable(fit.CandidateTable)
    T=fit.CandidateTable; names=T.Properties.VariableNames; eo=find_name(names,{'EO','eo','EO_id'}); rm=find_name(names,{'plainRmseMv','rmse_mv','rmse'});
    if ~isempty(eo) && ~isempty(rm)
        fig=figure('Visible','off','Color','w'); plot(T{:,eo},T{:,rm},'ko-','LineWidth',1.1); grid on; xlabel('EO'); ylabel('candidate RMSE (mV)'); title(sprintf('Window %d EO candidate profile',windowIndex));
        exportgraphics(fig,fullfile(outputDir,'B3_eo_profile.png'),'Resolution',160); close(fig);
    end
end
report.fit=fit; report.outputDir=outputDir; save(fullfile(outputDir,'DynamicSingleWindowReport.mat'),'report','-v7.3');
end
function R=pick_result(S), if isfield(S,'Result'), R=S.Result; elseif isfield(S,'out'), R=S.out; else, error('Unknown result schema.'); end, end
function W=pick_windows(R), if isfield(R,'WindowResult'), W=R.WindowResult; elseif isfield(R,'rows'), W=R.rows; else, error('No windows found.'); end, end
function f=pick_fit(w), f=w; if isfield(w,'modelFits') && isfield(w.modelFits,'gap_only'), f=w.modelFits.gap_only; end, end
function B=pick_bundle(w), if isfield(w,'bundle'), B=w.bundle; elseif isfield(w,'CoreBundlePreview'), B=w.CoreBundlePreview; else, error('Window has no waveform bundle.'); end, end
function B=normalize_bundle(B)
% Accept both public R5 and legacy fixed-gap bundle field names.
if ~isfield(B,'sensorIds') && isfield(B,'sensor_id'), B.sensorIds=unique(double(B.sensor_id(:))).'; end
if ~isfield(B,'sensorIndex') && isfield(B,'sensor_id')
    [~,~,B.sensorIndex]=unique(double(B.sensor_id(:)),'stable');
end
if ~isfield(B,'X') && isfield(B,'x'), B.X=B.x; end
if ~isfield(B,'V') && isfield(B,'v'), B.V=B.v; end
if ~isfield(B,'T') && isfield(B,'t'), B.T=B.t; end
end
function n=find_name(names,c), n=''; for i=1:numel(c), if ismember(c{i},names), n=c{i}; return; end, end, end
function plot_segmented(x,y,varargin)
x=double(x(:)); y=double(y(:)); ok=isfinite(x)&isfinite(y); x=x(ok); y=y(ok); if isempty(x),return;end
span=max(x)-min(x); if span<=0, plot(x,y,varargin{:}); return; end
dx=diff(x); cuts=[0;find((dx < -0.25*span) | (abs(dx)>max(0.25,0.35*span)));numel(x)]; first=true;
for k=1:numel(cuts)-1
    ii=(cuts(k)+1):cuts(k+1); if numel(ii)<2,continue;end
    [~,ord]=sort(x(ii)); ii=ii(ord);
    a=varargin; if ~first,a(end+1:end+2)={'HandleVisibility','off'};end
    plot(x(ii),y(ii),a{:}); first=false;
end
end
