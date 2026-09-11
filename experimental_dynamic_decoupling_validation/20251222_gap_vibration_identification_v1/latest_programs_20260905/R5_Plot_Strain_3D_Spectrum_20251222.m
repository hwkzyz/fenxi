function outputFiles = R5_Plot_Strain_3D_Spectrum_20251222(regionTag)
if nargin < 1 || isempty(regionTag), regionTag = 'R04'; end
thisDir = fileparts(mfilename('fullpath')); addpath(genpath(thisDir));
P = Setup_Paths_20251222(); outDir = fullfile(P.paths.strainValidation,'strain_3d_spectrum');
if ~isfolder(outDir), mkdir(outDir); end
step04File=fullfile(P.paths.prepared,'Step04_BTT_STE_Resonance_Regions_20251222.mat'); A4=load(step04File,'alignResult','strainSource'); bestTauSec=A4.alignResult.best_tau_sec; rawStrainOffsetSec=A4.strainSource.strain_time_offset_total_sec;
switch upper(string(regionTag))
    case "R04", resultFile=fullfile(P.paths.gapResults,'Main_GapAware_VPTopK_FullWave_20251222_B5_S123_r04_gapaware.mat'); timePad=[31.4 33.2]; expectedHz=629.86; bandHz=[600 650]; tag='B5_R04_BTTAligned';
    case "R01", resultFile=fullfile(P.paths.gapResults,'Main_GapAware_VPTopK_FullWave_20251222_B1_S123_r01_gapaware.mat'); timePad=[49.7 51.6]; expectedHz=582.22; bandHz=[560 610]; tag='B1_R01_BTTAligned';
    otherwise, error('Use R01 or R04.');
end
R=load(resultFile,'Result'); W=R.Result.WindowResult(:); ws=arrayfun(@(w)w.timeWindow(1),W); we=arrayfun(@(w)w.timeWindow(2),W);
sourceFile='E:\试验数据\20251222\应变片数据\AI1-03_20251222204737.mat'; outputFiles=strings(0,1);
for channel=[1 3]
    fileNow=regexprep(sourceFile,'AI1-03',sprintf('AI1-%02d',channel)); S=load(fileNow,'Datas'); t=double(S.Datas(:,1))+rawStrainOffsetSec-bestTauSec; y=double(S.Datas(:,2)); use=t>=timePad(1)&t<=timePad(2); t=t(use); y=detrend(y(use)); Fs=1/median(diff(t));
    winLen=max(256,round(0.20*Fs)); noverlap=round(0.90*winLen); nfft=max(4096,2^nextpow2(winLen)); [Z,F,T]=spectrogram(y,hann(winLen),noverlap,nfft,Fs); T=T+timePad(1); A=abs(Z)*2/sum(hann(winLen)); keep=F>=500&F<=700; Fp=F(keep); Ap=A(keep,:); [TT,FF]=meshgrid(T,Fp);
    fig=figure('Visible','off','Color','w','Position',[80 80 1250 760]); surf(TT,FF,Ap,'EdgeColor','none'); shading interp; view(43,58); axis tight; box on; grid on; xlabel('Time (s)'); ylabel('Frequency (Hz)'); zlabel('Strain amplitude (\mu\epsilon)'); title(sprintf('%s | AI1-%02d | 3-D strain spectrum | tau=%.3f s, raw offset=%.3f s',tag,channel,bestTauSec,rawStrainOffsetSec),'Interpreter','none'); colormap(turbo); colorbar; hold on; zmax=max(Ap(:),[],'omitnan'); plot3([timePad(1) timePad(2)],[expectedHz expectedHz],[zmax zmax],'k--','LineWidth',1.4); plot3([timePad(1) timePad(2)],[bandHz(1) bandHz(1)],[zmax*.98 zmax*.98],'w:','LineWidth',1); plot3([timePad(1) timePad(2)],[bandHz(2) bandHz(2)],[zmax*.98 zmax*.98],'w:','LineWidth',1); for k=1:numel(ws), if we(k)>=timePad(1)&&ws(k)<=timePad(2), plot3([ws(k) ws(k)],[500 700],[zmax zmax],'Color',[.15 .15 .15 .35],'LineWidth',.5); end, end; legend({'spectrogram','registered resonance','band limits','BTT windows'},'Location','northeast','Box','off'); base=fullfile(outDir,sprintf('Strain3D_%s_AI1-%02d',tag,channel)); exportgraphics(fig,[base '.png'],'Resolution',220); savefig(fig,[base '.fig']); close(fig); outputFiles(end+1)=string([base '.png']);
end
fprintf('3-D strain spectra saved under %s\n',outDir);
end
