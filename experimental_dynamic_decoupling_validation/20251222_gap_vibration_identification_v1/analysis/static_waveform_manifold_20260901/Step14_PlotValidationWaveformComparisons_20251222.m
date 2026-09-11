%% Representative waveform comparisons for V1-V4
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath')); outDir = fullfile(thisDir, 'results');
U = load(fullfile(outDir, 'StaticWaveformFamily_20251222.mat')); F = U.StaticWaveformFamily;
Y = F.waveformsB2FixedPeakAlignedMv; x = F.xB2OprAxisMm(:); gaps = F.nominalGapMm(:); blades = F.bladeIds(:).';
v1 = readtable(fullfile(outDir, 'V1_SameBladeTransition_20251222.csv'));
v2 = readtable(fullfile(outDir, 'V2_LeaveOneGap_20251222.csv'));
v3 = readtable(fullfile(outDir, 'V3_LeaveOneBlade_IncrementTransfer_20251222.csv'));
v4 = readtable(fullfile(outDir, 'V4_ExperimentalClosure_PerWindow_20251222.csv'));

% Select a representative case closest to the median gate metric.
i1 = closest_to(v1.relativeRMSEpct, median(v1.relativeRMSEpct));
ib = find(blades == v1.blade(i1), 1); il = find(abs(gaps-v1.lowGapMm(i1)) < 1e-9, 1); it = find(abs(gaps-v1.targetGapMm(i1)) < 1e-9, 1);
yt1 = squeeze(Y(:,ib,it)); C1 = fit_gap(squeeze(Y(:,ib,:)), gaps, it); yh1 = squeeze(Y(:,ib,il)) + surface(C1, gaps(it)) - surface(C1, gaps(il));

i2 = closest_to(v2.relativeRMSEpct, median(v2.relativeRMSEpct));
ib2 = find(blades == v2.blade(i2), 1); ih = find(abs(gaps-v2.heldOutGapMm(i2)) < 1e-9, 1);
yt2 = squeeze(Y(:,ib2,ih)); yh2 = surface(fit_gap(squeeze(Y(:,ib2,:)), gaps, ih), gaps(ih));

i3 = closest_to(v3.relativeToBaseline, median(v3.relativeToBaseline));
ih3 = find(blades == v3.heldOutBlade(i3), 1); il3 = find(abs(gaps-v3.lowGapMm(i3)) < 1e-9, 1); it3 = find(abs(gaps-v3.targetGapMm(i3)) < 1e-9, 1);
others = setdiff(1:numel(blades), ih3);
C3 = fit_increment_model(Y(:,others,:), gaps);
d3 = predict_increment_model(C3, gaps(il3), gaps(it3)); y03 = squeeze(Y(:,ih3,il3)); yt3 = squeeze(Y(:,ih3,it3)); yh3 = y03+d3;

fig = figure('Color','w','Units','centimeters','Position',[2 2 17 12.5]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
wave_panel(x, yt1, yh1, sprintf('V1 | B%d: %.2f -> %.2f mm', v1.blade(i1), v1.lowGapMm(i1), v1.targetGapMm(i1)), 'Measured target','Increment prediction');
wave_panel(x, yt2, yh2, sprintf('V2 | B%d | held-out %.2f mm', v2.blade(i2), v2.heldOutGapMm(i2)), 'Measured target','Leave-one-gap prediction');
wave_panel(x, yt3, yh3, sprintf('V3 | held-out B%d: %.2f -> %.2f mm', v3.heldOutBlade(i3), v3.lowGapMm(i3), v3.targetGapMm(i3)), 'Measured target','Cross-blade increment transfer');
nexttile; yyaxis left; plot(v4.windowId, abs(v4.frequencyDiffHz_freeMinusNoB), 'o-','Color',[0.12 0.35 0.62],'MarkerFaceColor',[0.12 0.35 0.62]); yline(0.05,'--','Color',[0.70 0.14 0.12]); ylabel('|df| (Hz)');
yyaxis right; plot(v4.windowId, abs(v4.RMSEDiffmV_freeMinusNoB), 's-','Color',[0.85 0.38 0.10],'MarkerFaceColor',[0.85 0.38 0.10]); yline(0.10,'--','Color',[0.70 0.14 0.12]); ylabel('|dRMSE| (mV)'); xlabel('High-speed window'); title('V4 | closure across 18 windows','FontWeight','normal'); box on; grid off; set(gca,'TickDir','in');
export_figure(fig, fullfile(outDir,'ValidationWaveformComparisons_V1V4_20251222'));
fprintf('Saved representative waveform comparison figure under %s\n', outDir);

function wave_panel(x, target, predicted, titleText, l1, l2)
nexttile; plot(x,target,'k-','LineWidth',1.15); hold on; plot(x,predicted,'Color',[0.12 0.35 0.62],'LineWidth',1.0); hold off; box on; grid off; set(gca,'TickDir','in'); xlabel('Peak-aligned position (mm)'); ylabel('Voltage (mV)'); title(titleText,'FontWeight','normal'); legend({l1,l2},'Location','best','Box','off','FontSize',7.5);
end
function C = fit_gap(Y,gaps,held), tr=setdiff(1:numel(gaps),held); C=nan(size(Y,1),3); for ix=1:size(Y,1), y=Y(ix,tr); ok=isfinite(y); if nnz(ok)>=4, A=[ones(nnz(ok),1),1./gaps(tr(ok)),log(gaps(tr(ok)))]; C(ix,:)=A\y(ok).'; end, end, end
function y=surface(C,g), y=C(:,1)+C(:,2)./g+C(:,3).*log(g); end
function i=closest_to(v,x), [~,i]=min(abs(v-x)); end
function C=fit_increment_model(Ytrain,gaps)
nX=size(Ytrain,1); nB=size(Ytrain,2); nG=numel(gaps); C=nan(nX,3); A=[];
for il=1:nG-1, it=il+1; A(end+1,:)=[1,1/gaps(it)-1/gaps(il),log(gaps(it)/gaps(il))]; end
AA=repmat(A,nB,1);
for ix=1:nX, y=[]; for il=1:nG-1, y=[y; squeeze(Ytrain(ix,:,il+1)-Ytrain(ix,:,il)).']; end; ok=isfinite(y); if nnz(ok)>=3, C(ix,:)=AA(ok,:)\y(ok); end, end
end
function d=predict_increment_model(C,g1,g2), d=C*[1;1/g2-1/g1;log(g2/g1)]; end
function export_figure(fig,baseName), if exist('exportgraphics','file'), exportgraphics(fig,[baseName '.png'],'Resolution',300); exportgraphics(fig,[baseName '.pdf'],'ContentType','vector'); else, print(fig,[baseName '.png'],'-dpng','-r300'); print(fig,[baseName '.pdf'],'-dpdf','-painters'); end, end
