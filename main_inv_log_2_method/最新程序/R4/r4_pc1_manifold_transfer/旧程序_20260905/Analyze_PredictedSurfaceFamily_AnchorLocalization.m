%% Predicted surface family + one-gap anchor localization (tilt hidden).
% A complete target tilt surface is removed.  Training surfaces are fit by
% the within-surface gap law, compressed by SVD, and linearly parameterized
% along their first latent coordinate.  The target anchor trace alone then
% localizes a continuous predicted surface; the remaining gaps are predicted
% from that surface.  Tilt labels are never used by the algorithm.

thisDir=fileparts(mfilename('fullpath'));
dataPath=fullfile(thisDir,'output','waveform_families','TiltGap_WaveformSurfaceData.mat');
outDir=fullfile(thisDir,'output','surface_identification_gap_law');
if ~exist(outDir,'dir'), mkdir(outDir); end
D=load(dataPath,'Data'); D=D.Data;
x=D.x_mm(:); gap=D.gap_mm(:); tilt=D.tilt_deg(:); V=D.response_V;
[nX,nGap,nTilt]=size(V); gref=0.8;
X=[ones(nGap,1),1./gap,log(gap/gref)]; nCoeff=3*nX;
anchorGap=zeros(nTilt*nGap,1); targetTilt=anchorGap; latentTarget=anchorGap;
latentHat=anchorGap; anchorRmse=anchorGap; transferRmse=anchorGap;
transferNrmsePct=anchorGap; transferCorr=anchorGap; peakXMae=anchorGap;
peakAmpMae=anchorGap; oracleTransferNrmsePct=anchorGap; row=0;

for it=1:nTilt
    trainTilt=setdiff(1:nTilt,it);
    % q(:,k) is the fitted (B0,B1,B2) vector for training surface k.
    Q=zeros(nCoeff,numel(trainTilt));
    for k=1:numel(trainTilt)
        B=(X \ V(:,:,trainTilt(k)).').';
        Q(:,k)=B(:);
    end
    mu=mean(Q,2); [U,Sv,~]=svd(Q-mu,'econ');
    r=min(2,size(U,2)); Ur=U(:,1:r); Z=Ur'*(Q-mu);
    % Use the first latent coordinate only to define the ordered surface path.
    [z1,ord]=sort(Z(1,:),'ascend'); Qord=Q(:,ord);
    targetQ=(X \ V(:,:,it).').'; targetQ=targetQ(:);
    % Target latent coordinate is used only after the blind prediction, for
    % diagnostics; it is never supplied to candidate generation/localization.
    latentTarget(it)=Ur(:,1)'*(targetQ-mu);
    for ia=1:nGap
        row=row+1; targetTilt(row)=tilt(it); anchorGap(row)=gap(ia);
        y=V(:,ia,it); trainGap=setdiff(1:nGap,ia);
        % Candidate latent coordinate: dense interval plus endpoint extension.
        dz=max(diff(z1)); tmin=z1(1)-0.5*dz; tmax=z1(end)+0.5*dz;
        tgrid=linspace(tmin,tmax,401); nT=numel(tgrid);
        predAnchor=zeros(nX,nT);
        for kt=1:nT
            q=interp1(z1,Qord.',tgrid(kt),'linear','extrap').';
            B=reshape(q,nX,3);
            predAnchor(:,kt)=B*X(ia,:).';
        end
        d=sqrt(mean((predAnchor-y).^2,1)); [~,ih]=min(d);
        latentHat(row)=tgrid(ih); anchorRmse(row)=d(ih);
        qhat=interp1(z1,Qord.',latentHat(row),'linear','extrap').';
        Bhat=reshape(qhat,nX,3); pred=Bhat*X.';
        truth=V(:,trainGap,it); predOther=pred(:,trainGap); e=predOther-truth;
        transferRmse(row)=sqrt(mean(e(:).^2));
        transferNrmsePct(row)=100*transferRmse(row)/range(truth(:));
        transferCorr(row)=corr(truth(:),predOther(:));
        [peakXMae(row),peakAmpMae(row)]=peak_errors(x,truth,predOther);
        Btrue=(X(trainGap,:) \ V(:,trainGap,it).').'; ptrue=Btrue*X.';
        eo=ptrue(:,trainGap)-truth;
        oracleTransferNrmsePct(row)=100*sqrt(mean(eo(:).^2))/range(truth(:));
    end
end

T=table(targetTilt,anchorGap,latentHat,anchorRmse,transferRmse, ...
    transferNrmsePct,transferCorr,peakXMae,peakAmpMae,oracleTransferNrmsePct);
writetable(T,fullfile(outDir,'predicted_family_localization_per_anchor.csv'));
summary=table(tilt,zeros(nTilt,1),zeros(nTilt,1),zeros(nTilt,1), ...
    zeros(nTilt,1),zeros(nTilt,1),zeros(nTilt,1), ...
    'VariableNames',{'target_tilt_deg','mean_anchor_RMSE_V', ...
    'mean_transfer_NRMSE_pct','worst_transfer_NRMSE_pct', ...
    'mean_transfer_correlation','mean_SG_peak_x_MAE_mm', ...
    'mean_SG_peak_amp_MAE_V'});
for it=1:nTilt
    rr=T.targetTilt==tilt(it);
    summary.mean_anchor_RMSE_V(it)=mean(T.anchorRmse(rr));
    summary.mean_transfer_NRMSE_pct(it)=mean(T.transferNrmsePct(rr));
    summary.worst_transfer_NRMSE_pct(it)=max(T.transferNrmsePct(rr));
    summary.mean_transfer_correlation(it)=mean(T.transferCorr(rr));
    summary.mean_SG_peak_x_MAE_mm(it)=mean(T.peakXMae(rr));
    summary.mean_SG_peak_amp_MAE_V(it)=mean(T.peakAmpMae(rr));
end
writetable(summary,fullfile(outDir,'predicted_family_localization_tilt_summary.csv'));
overall=table(mean(T.anchorRmse),mean(T.transferNrmsePct),max(T.transferNrmsePct), ...
    mean(T.transferCorr),mean(T.peakXMae),mean(T.peakAmpMae), ...
    mean(T.oracleTransferNrmsePct), ...
    'VariableNames',{'mean_anchor_RMSE_V','mean_transfer_NRMSE_pct', ...
    'worst_transfer_NRMSE_pct','mean_transfer_correlation', ...
    'mean_SG_peak_x_MAE_mm','mean_SG_peak_amp_MAE_V', ...
    'mean_oracle_transfer_NRMSE_pct'});
writetable(overall,fullfile(outDir,'predicted_family_localization_summary.csv'));
save(fullfile(outDir,'PredictedSurfaceFamilyLocalizationResults.mat'), ...
    'T','summary','overall','x','gap','tilt','gref','dataPath');
fprintf('\nPredicted surface family + blind anchor localization\n'); disp(overall); disp(summary);
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
