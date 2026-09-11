%% Strict V3 waveform alignment audit.
% Fixed-coordinate panels are the only panels used for V3 errors. The
% peak-aligned panels are diagnostic only and must not be used for scoring.
clear; clc;
thisDir=fileparts(mfilename('fullpath')); outDir=fullfile(thisDir,'results');
figDir=fullfile(outDir,'StrictV3_WaveformAlignmentAudit_20251222');
if ~exist(figDir,'dir'), mkdir(figDir); end
A=load(fullfile(outDir,'StrictV3_WaveformCases_20251222.mat'));
x=A.x(:);
for im=1:numel(A.methodNames)
    p=nan(1,size(A.predCase,2));
    for k=1:size(A.predCase,2)
        ok=isfinite(A.predCase(:,k,im)) & isfinite(A.targetCase(:,k));
        if nnz(ok)>20, p(k)=sqrt(mean((A.predCase(ok,k,im)-A.targetCase(ok,k)).^2)); end
    end
    [~,ord]=sort(p,'descend','MissingPlacement','last');
    picks=unique([ord(1:min(6,numel(ord))),round(linspace(1,numel(ord),3))],'stable');
    picks=picks(isfinite(p(picks))); picks=picks(1:min(9,numel(picks)));
    f=figure('Visible','off','Color','w','Position',[60 60 1600 1000]);
    tl=tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
    for jj=1:numel(picks)
        k=picks(jj); yp=A.predCase(:,k,im); yt=A.targetCase(:,k); yl=A.lowCase(:,k);
        ok=isfinite(yp)&isfinite(yt)&isfinite(yl); xx=x(ok);
        % Common B2 coordinate: formal V3 comparison.
        nexttile; plot(xx,yt(ok),'k-','LineWidth',1.15); hold on;
        plot(xx,yp(ok),'r--','LineWidth',1.0); plot(xx,yl(ok),'b:','LineWidth',0.9);
        [~,it]=max(yt(ok)); [~,ip]=max(yp(ok)); [~,il]=max(yl(ok));
        plot(xx(it),yt(ok(it)),'ko','MarkerSize',4,'HandleVisibility','off');
        plot(xx(ip),yp(ok(ip)),'rs','MarkerSize',4,'HandleVisibility','off');
        plot(xx(il),yl(ok(il)),'bd','MarkerSize',4,'HandleVisibility','off');
        grid on; box on;
        title(sprintf('B%d: %.2f -> %.2f mm | RMSE %.1f mV',A.caseMeta(k,1),A.caseMeta(k,2),A.caseMeta(k,3),p(k)),'FontSize',8);
        if jj>6, xlabel('fixed B2 x (mm)'); end; if mod(jj-1,3)==0, ylabel('mV'); end
    end
    ax=findall(f,'Type','axes'); lgd=legend(ax(1),{'target (black)','prediction (red)','low state (blue)'},'Location','southoutside','Orientation','horizontal'); lgd.Box='off';
    sgtitle([strrep(A.methodNames{im},'_',' ') ' | fixed B2 coordinate'],'Interpreter','none');
    exportgraphics(f,fullfile(figDir,[A.methodNames{im} '_fixedB2.png']),'Resolution',180); close(f);

    % Diagnostic only: translate each curve to its own SG-like discrete peak.
    f=figure('Visible','off','Color','w','Position',[60 60 1600 1000]);
    tl=tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
    for jj=1:numel(picks)
        k=picks(jj); yp=A.predCase(:,k,im); yt=A.targetCase(:,k); yl=A.lowCase(:,k);
        ok=isfinite(yp)&isfinite(yt)&isfinite(yl); xx=x(ok);
        curves={yt(ok),yp(ok),yl(ok)}; names={'target','prediction','low'}; cols={'k-','r--','b:'};
        nexttile; hold on;
        for ic=1:3
            [~,ii]=max(curves{ic}); xq=xx-xx(ii); plot(xq,curves{ic},cols{ic},'LineWidth',1.0);
        end
        grid on; box on; title(sprintf('B%d: peak-aligned diagnostic | fixed RMSE %.1f mV',A.caseMeta(k,1),p(k)),'FontSize',8);
        if jj>6, xlabel('x - own peak (mm)'); end; if mod(jj-1,3)==0, ylabel('mV'); end
    end
    ax=findall(f,'Type','axes'); lgd=legend(ax(1),{'target','prediction','low state'},'Location','southoutside','Orientation','horizontal'); lgd.Box='off';
    sgtitle([strrep(A.methodNames{im},'_',' ') ' | diagnostic peak alignment (not scored)'],'Interpreter','none');
    exportgraphics(f,fullfile(figDir,[A.methodNames{im} '_peakDiagnostic.png']),'Resolution',180); close(f);
end
fprintf('Wrote fixed-coordinate and diagnostic alignment figures for %d methods to %s\n',numel(A.methodNames),figDir);
