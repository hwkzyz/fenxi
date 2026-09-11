%% Strict V3 plots with one common, low-state-derived alignment per case.
% The same shift is applied to low, target, and prediction. No target-only
% realignment is used, so the plotted error remains a valid held-out error.
clear; clc;
thisDir=fileparts(mfilename('fullpath')); outDir=fullfile(thisDir,'results');
figDir=fullfile(outDir,'StrictV3_WaveformCommonPeakAlignment_20251222');
if ~exist(figDir,'dir'), mkdir(figDir); end
A=load(fullfile(outDir,'StrictV3_WaveformCases_20251222.mat')); x=A.x(:);
for im=1:numel(A.methodNames)
    p=nan(1,size(A.predCase,2));
    for k=1:size(A.predCase,2)
        ok=isfinite(A.predCase(:,k,im))&isfinite(A.targetCase(:,k));
        if nnz(ok)>20, p(k)=sqrt(mean((A.predCase(ok,k,im)-A.targetCase(ok,k)).^2)); end
    end
    [~,ord]=sort(p,'descend','MissingPlacement','last'); picks=unique([ord(1:min(6,numel(ord))),round(linspace(1,numel(ord),3))],'stable'); picks=picks(isfinite(p(picks))); picks=picks(1:min(9,numel(picks)));
    f=figure('Visible','off','Color','w','Position',[60 60 1600 1000]); tl=tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
    for jj=1:numel(picks)
        k=picks(jj); yp=A.predCase(:,k,im); yt=A.targetCase(:,k); yl=A.lowCase(:,k); ok=isfinite(yp)&isfinite(yt)&isfinite(yl); xx=x(ok); yy=yl(ok);
        ys=smoothdata(yy,'sgolay',31); [~,i0]=max(ys); x0=xx(i0); xa=xx-x0;
        nexttile; plot(xa,yt(ok),'k-','LineWidth',1.15); hold on; plot(xa,yp(ok),'r--','LineWidth',1.0); plot(xa,yl(ok),'b:','LineWidth',0.9); xline(0,'Color',[.5 .5 .5],'HandleVisibility','off'); grid on; box on;
        title(sprintf('B%d: %.2f -> %.2f mm | RMSE %.1f mV',A.caseMeta(k,1),A.caseMeta(k,2),A.caseMeta(k,3),p(k)),'FontSize',8);
        if jj>6, xlabel('x - low-state SG peak (mm)'); end; if mod(jj-1,3)==0, ylabel('mV'); end
    end
    ax=findall(f,'Type','axes'); lgd=legend(ax(1),{'target','prediction','low state'},'Location','southoutside','Orientation','horizontal'); lgd.Box='off';
    sgtitle([strrep(A.methodNames{im},'_',' ') ' | common low-state peak alignment'],'Interpreter','none');
    exportgraphics(f,fullfile(figDir,[A.methodNames{im} '_commonLowPeak.png']),'Resolution',180); close(f);
end
fprintf('Wrote common-low-peak aligned figures for %d methods to %s\n',numel(A.methodNames),figDir);
