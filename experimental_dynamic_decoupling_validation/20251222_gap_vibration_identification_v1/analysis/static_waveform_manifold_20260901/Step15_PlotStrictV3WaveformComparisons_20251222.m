%% Plot target versus prediction for every strict V3 method.
clear; clc;
thisDir=fileparts(mfilename('fullpath')); outDir=fullfile(thisDir,'results'); figDir=fullfile(outDir,'StrictV3_WaveformComparisons_20251222');
if ~exist(figDir,'dir'), mkdir(figDir); end
A=load(fullfile(outDir,'StrictV3_WaveformCases_20251222.mat'));
for im=1:numel(A.methodNames)
    f=figure('Visible','off','Color','w','Position',[80 80 1500 900]); tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
    % Show representative transitions plus the worst cases under this method.
    p=[];
    for k=1:size(A.predCase,2)
        y=A.predCase(:,k,im); t=A.targetCase(:,k); ok=isfinite(y)&isfinite(t);
        p(k)=sqrt(mean((y(ok)-t(ok)).^2));
    end
    [~,ord]=sort(p,'descend'); picks=unique([ord(1:min(6,numel(ord))), round(linspace(1,numel(ord),3))],'stable'); picks=picks(1:min(9,numel(picks)));
    for jj=1:numel(picks)
        k=picks(jj); nexttile; plot(A.x,A.targetCase(:,k),'k-','LineWidth',1.1); hold on; plot(A.x,A.predCase(:,k,im),'r--','LineWidth',1.0); plot(A.x,A.lowCase(:,k),'b:','LineWidth',0.9); grid on; box on;
        title(sprintf('B%d: %.2f -> %.2f mm | RMSE %.1f mV',A.caseMeta(k,1),A.caseMeta(k,2),A.caseMeta(k,3),p(k)),'FontSize',8);
        if jj>6, xlabel('x (mm)'); end; if mod(jj-1,3)==0, ylabel('mV'); end
    end
    lgd=legend('target','prediction','low state','Location','southoutside','Orientation','horizontal'); lgd.Box='off';
    sgtitle(strrep(A.methodNames{im},'_',' '),'Interpreter','none');
    exportgraphics(f,fullfile(figDir,[A.methodNames{im} '.png']),'Resolution',180); close(f);
end
fprintf('Wrote %d strict V3 method comparison figures to %s\n',numel(A.methodNames),figDir);
