%% Open-set transfer: the target tilt surface is absent from the candidate library.
% One anchor-gap trace from the hidden target tilt selects the closest
% candidate surface. The selected candidate's within-surface gap law then
% predicts the target tilt's remaining gaps. No target tilt label is used.

thisDir = fileparts(mfilename('fullpath'));
dataPath = fullfile(thisDir,'output','waveform_families','TiltGap_WaveformSurfaceData.mat');
outDir = fullfile(thisDir,'output','surface_identification_gap_law');
if ~exist(outDir,'dir'), mkdir(outDir); end
D = load(dataPath,'Data'); D=D.Data;
x=D.x_mm(:); gap=D.gap_mm(:); tilt=D.tilt_deg(:); V=D.response_V;
[nX,nGap,nTilt]=size(V); gref=0.8;
X=[ones(nGap,1),1./gap,log(gap/gref)];
targetTilt=zeros(nTilt*nGap,1); anchorGap=targetTilt; selectedTilt=targetTilt;
candidateDistance=targetTilt; transferRmse=targetTilt; transferNrmsePct=targetTilt;
transferCorr=targetTilt; peakXMae=targetTilt; peakAmpMae=targetTilt;
oracleTransferNrmsePct=targetTilt; row=0;
incrementRmse=targetTilt; incrementNrmsePct=targetTilt; incrementCorr=targetTilt;
incrementPeakXMae=targetTilt; incrementPeakAmpMae=targetTilt;

for it=1:nTilt
    candidates=setdiff(1:nTilt,it);
    for ia=1:nGap
        row=row+1; targetTilt(row)=tilt(it); anchorGap(row)=gap(ia);
        trainGap=setdiff(1:nGap,ia);
        predAnchor=nan(nX,numel(candidates)); B=cell(numel(candidates),1);
        for k=1:numel(candidates)
            ic=candidates(k);
            % Candidate surfaces are calibration data and may use all gaps;
            % only the hidden target tilt is excluded from the library.
            B{k}=(X \ V(:,:,ic).').';
            predAnchor(:,k)=B{k}*X(ia,:).';
        end
        y=V(:,ia,it); d=sqrt(mean((predAnchor-y).^2,1));
        [candidateDistance(row),kbest]=min(d); selectedTilt(row)=tilt(candidates(kbest));
        pred=B{kbest}*X.'; truth=V(:,trainGap,it); predOther=pred(:,trainGap);
        % Incremental transfer: retain the observed anchor trace and add the
        % candidate surface's predicted gap increment.
        predIncrement = y + pred(:,trainGap) - pred(:,ia);
        ei=predIncrement-truth; incrementRmse(row)=sqrt(mean(ei(:).^2));
        incrementNrmsePct(row)=100*incrementRmse(row)/range(truth(:));
        incrementCorr(row)=corr(truth(:),predIncrement(:));
        [incrementPeakXMae(row),incrementPeakAmpMae(row)]=peak_errors(x,truth,predIncrement);
        e=predOther-truth; transferRmse(row)=sqrt(mean(e(:).^2));
        transferNrmsePct(row)=100*transferRmse(row)/range(truth(:));
        transferCorr(row)=corr(truth(:),predOther(:));
        [peakXMae(row),peakAmpMae(row)]=peak_errors(x,truth,predOther);
        % Oracle uses the hidden target surface's own gap law, only as a
        % decomposition reference and never for the blind selection.
        Btrue=(X(trainGap,:) \ V(:,trainGap,it).').'; ptrue=Btrue*X.';
        eo=ptrue(:,trainGap)-truth;
        oracleTransferNrmsePct(row)=100*sqrt(mean(eo(:).^2))/range(truth(:));
    end
end

T=table(targetTilt,anchorGap,selectedTilt,candidateDistance,transferRmse, ...
    transferNrmsePct,transferCorr,peakXMae,peakAmpMae,incrementRmse, ...
    incrementNrmsePct,incrementCorr,incrementPeakXMae,incrementPeakAmpMae, ...
    oracleTransferNrmsePct);
writetable(T,fullfile(outDir,'open_set_per_anchor.csv'));
summary=table(tilt,zeros(nTilt,1),zeros(nTilt,1),zeros(nTilt,1), ...
    zeros(nTilt,1),zeros(nTilt,1),zeros(nTilt,1), ...
    'VariableNames',{'target_tilt_deg','mean_candidate_distance_V', ...
    'mean_transfer_NRMSE_pct','worst_transfer_NRMSE_pct', ...
    'mean_transfer_correlation','mean_SG_peak_x_MAE_mm', ...
    'mean_SG_peak_amp_MAE_V'});
for it=1:nTilt
    r=T.targetTilt==tilt(it);
    summary.mean_candidate_distance_V(it)=mean(T.candidateDistance(r));
    summary.mean_transfer_NRMSE_pct(it)=mean(T.transferNrmsePct(r));
    summary.worst_transfer_NRMSE_pct(it)=max(T.transferNrmsePct(r));
    summary.mean_transfer_correlation(it)=mean(T.transferCorr(r));
    summary.mean_SG_peak_x_MAE_mm(it)=mean(T.peakXMae(r));
    summary.mean_SG_peak_amp_MAE_V(it)=mean(T.peakAmpMae(r));
end
writetable(summary,fullfile(outDir,'open_set_tilt_summary.csv'));
overall=table(mean(T.candidateDistance),mean(T.transferNrmsePct), ...
    max(T.transferNrmsePct),mean(T.transferCorr),mean(T.peakXMae), ...
    mean(T.peakAmpMae),mean(T.incrementNrmsePct),max(T.incrementNrmsePct), ...
    mean(T.incrementCorr),mean(T.incrementPeakXMae),mean(T.incrementPeakAmpMae), ...
    mean(T.oracleTransferNrmsePct), ...
    'VariableNames',{'mean_candidate_distance_V','mean_transfer_NRMSE_pct', ...
    'worst_transfer_NRMSE_pct','mean_transfer_correlation', ...
    'mean_SG_peak_x_MAE_mm','mean_SG_peak_amp_MAE_V', ...
    'mean_increment_NRMSE_pct','worst_increment_NRMSE_pct', ...
    'mean_increment_correlation','mean_increment_peak_x_MAE_mm', ...
    'mean_increment_peak_amp_MAE_V', ...
    'mean_oracle_transfer_NRMSE_pct'});
writetable(overall,fullfile(outDir,'open_set_summary.csv'));
save(fullfile(outDir,'OpenSetSurfaceTransferGapLawResults.mat'),'T','summary','overall', ...
    'x','gap','tilt','gref','dataPath');
fprintf('\nOpen-set surface transfer (target tilt excluded)\n'); disp(overall); disp(summary);
fprintf('Outputs written to:\n%s\n',outDir);

function [xm,am]=peak_errors(x,truth,pred)
n=size(truth,2); xt=zeros(n,1); at=xt; xp=xt; ap=xt;
for j=1:n
    [at(j),q]=sg_peak(truth(:,j)); xt(j)=x(q);
    [ap(j),q]=sg_peak(pred(:,j)); xp(j)=x(q);
end
xm=mean(abs(xp-xt)); am=mean(abs(ap-at));
end
function [a,idx]=sg_peak(y)
y=sgolayfilt(y(:),3,11); [a,idx]=max(y);
end
