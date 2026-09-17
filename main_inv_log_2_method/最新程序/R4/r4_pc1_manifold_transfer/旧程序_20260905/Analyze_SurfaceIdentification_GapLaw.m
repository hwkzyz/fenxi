%% Surface identification from one gap waveform, followed by gap transfer.
% The tilt label is hidden during identification.  For each candidate tilt,
% the within-tilt gap law is fitted as B0(x)+B1(x)/g+B2(x)log(g/gref).
% The anchor gap is left out of every candidate fit to prevent leakage.

thisDir = fileparts(mfilename('fullpath'));
dataPath = fullfile(thisDir, 'output', 'waveform_families', ...
    'TiltGap_WaveformSurfaceData.mat');
outDir = fullfile(thisDir, 'output', 'surface_identification_gap_law');
if ~exist(outDir, 'dir'), mkdir(outDir); end

S = load(dataPath, 'Data');
Data = S.Data;
x = Data.x_mm(:); gap = Data.gap_mm(:); tilt = Data.tilt_deg(:);
V = Data.response_V;
[nX,nGap,nTilt] = size(V);
gref = 0.8;
Xgap = [ones(nGap,1), 1./gap, log(gap/gref)];

% Fit each candidate surface from the non-anchor gaps.  B is x-by-3.
anchorGap = zeros(nTilt*nGap,1); targetTilt = anchorGap;
selectedTilt = anchorGap; selectedRank = anchorGap;
anchorRmse = anchorGap; oracleAnchorRmse = anchorGap;
transferRmse = anchorGap; transferNrmsePct = anchorGap;
transferCorr = anchorGap; peakXMae = anchorGap; peakAmpMae = anchorGap;
oracleTransferRmse = anchorGap; oracleTransferNrmsePct = anchorGap;
row = 0;
allDistances = nan(nTilt*nGap,nTilt);

for iTarget = 1:nTilt
    for iAnchor = 1:nGap
        row = row + 1;
        targetTilt(row) = tilt(iTarget); anchorGap(row) = gap(iAnchor);
        trainGap = setdiff(1:nGap,iAnchor);
        B = cell(nTilt,1);
        predAnchor = nan(nX,nTilt);
        for iCandidate = 1:nTilt
            % Fit the candidate surface without the observed anchor trace.
            B{iCandidate} = (Xgap(trainGap,:) \ V(:,trainGap,iCandidate).').';
            predAnchor(:,iCandidate) = B{iCandidate} * Xgap(iAnchor,:).';
        end
        yAnchor = V(:,iAnchor,iTarget);
        d = sqrt(mean((predAnchor-yAnchor).^2,1));
        allDistances(row,:) = d;
        [~,order] = sort(d,'ascend');
        iSelected = order(1); selectedTilt(row) = tilt(iSelected);
        selectedRank(row) = find(order == iTarget,1);
        anchorRmse(row) = d(iSelected); oracleAnchorRmse(row) = d(iTarget);

        % Predict the five non-anchor gaps using the selected surface.
        predSelected = B{iSelected} * Xgap.';
        truthOther = V(:,trainGap,iTarget);
        err = predSelected(:,trainGap) - truthOther;
        transferRmse(row) = sqrt(mean(err(:).^2));
        transferNrmsePct(row) = 100*transferRmse(row)/range(truthOther(:));
        predOther = predSelected(:,trainGap);
        transferCorr(row) = corr(truthOther(:),predOther(:));
        [peakXMae(row),peakAmpMae(row)] = peak_errors(x,truthOther,predSelected(:,trainGap));

        % Oracle: correct candidate selected, with the same leave-anchor fit.
        oraclePred = B{iTarget} * Xgap.';
        eo = oraclePred(:,trainGap) - truthOther;
        oracleTransferRmse(row) = sqrt(mean(eo(:).^2));
        oracleTransferNrmsePct(row) = 100*oracleTransferRmse(row)/range(truthOther(:));
    end
end

correct = selectedTilt == targetTilt;
T = table(targetTilt,anchorGap,selectedTilt,selectedRank,correct,anchorRmse, ...
    oracleAnchorRmse,transferRmse,transferNrmsePct,transferCorr,peakXMae, ...
    peakAmpMae,oracleTransferRmse,oracleTransferNrmsePct);
writetable(T,fullfile(outDir,'surface_identification_per_anchor.csv'));

tiltSummary = table(tilt,zeros(nTilt,1),zeros(nTilt,1),zeros(nTilt,1), ...
    zeros(nTilt,1),zeros(nTilt,1),zeros(nTilt,1), ...
    'VariableNames',{'tilt_deg','identification_accuracy','mean_rank', ...
    'mean_transfer_NRMSE_pct','worst_transfer_NRMSE_pct', ...
    'mean_peak_x_MAE_mm','mean_peak_amp_MAE_V'});
for k=1:nTilt
    r = T.targetTilt == tilt(k);
    tiltSummary.identification_accuracy(k) = mean(T.correct(r));
    tiltSummary.mean_rank(k) = mean(T.selectedRank(r));
    tiltSummary.mean_transfer_NRMSE_pct(k) = mean(T.transferNrmsePct(r));
    tiltSummary.worst_transfer_NRMSE_pct(k) = max(T.transferNrmsePct(r));
    tiltSummary.mean_peak_x_MAE_mm(k) = mean(T.peakXMae(r));
    tiltSummary.mean_peak_amp_MAE_V(k) = mean(T.peakAmpMae(r));
end
writetable(tiltSummary,fullfile(outDir,'surface_identification_tilt_summary.csv'));

methodSummary = table(mean(T.correct),mean(T.transferNrmsePct), ...
    max(T.transferNrmsePct),mean(T.transferCorr),mean(T.peakXMae), ...
    mean(T.peakAmpMae),mean(T.oracleTransferNrmsePct), ...
    'VariableNames',{'overall_identification_accuracy','mean_transfer_NRMSE_pct', ...
    'worst_transfer_NRMSE_pct','mean_transfer_correlation', ...
    'mean_SG_peak_x_MAE_mm','mean_SG_peak_amp_MAE_V', ...
    'mean_oracle_transfer_NRMSE_pct'});
writetable(methodSummary,fullfile(outDir,'surface_identification_summary.csv'));

% Distance-margin diagnostic: positive margin means the true surface wins.
trueD = nan(height(T),1); bestWrongD = trueD;
for r=1:height(T)
    iTrue = find(tilt == T.targetTilt(r),1);
    trueD(r) = allDistances(r,iTrue);
    z = allDistances(r,:); z(iTrue)=Inf; bestWrongD(r)=min(z);
end
margin = bestWrongD-trueD;
writetable(table(T.targetTilt,T.anchorGap,trueD,bestWrongD,margin), ...
    fullfile(outDir,'surface_identification_distance_margin.csv'));

save(fullfile(outDir,'SurfaceIdentificationGapLawResults.mat'), ...
    'T','tiltSummary','methodSummary','allDistances','x','gap','tilt','gref','dataPath');

fprintf('\nOne-gap surface identification with leave-anchor gap-law fits\n');
disp(methodSummary);
fprintf('Overall identification accuracy: %.2f%% (%d/%d)\n', ...
    100*mean(T.correct),sum(T.correct),height(T));
disp(tiltSummary);
fprintf('Outputs written to:\n%s\n',outDir);

function [xMae,ampMae] = peak_errors(x,truth,pred)
nt = size(truth,2); xp=zeros(nt,1); ap=xp; xq=xp; aq=xp;
for j=1:nt
    [ap(j),it]=sg_peak(truth(:,j)); xp(j)=x(it);
    [aq(j),iq]=sg_peak(pred(:,j)); xq(j)=x(iq);
end
xMae=mean(abs(xq-xp)); ampMae=mean(abs(aq-ap));
end

function [a,idx] = sg_peak(y)
y=sgolayfilt(y(:),3,11); [a,idx]=max(y);
end
