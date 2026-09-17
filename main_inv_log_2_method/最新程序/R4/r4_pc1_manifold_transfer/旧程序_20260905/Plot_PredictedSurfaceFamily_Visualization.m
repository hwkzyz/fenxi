%% Visualization of predicted surface-family localization and gap transfer.
% Uses the same blind open-set protocol as the analysis script.  The target
% tilt is excluded from the surface library; only one anchor-gap waveform is
% used for localization.

thisDir=fileparts(mfilename('fullpath'));
dataPath=fullfile(thisDir,'output','waveform_families','TiltGap_WaveformSurfaceData.mat');
outDir=fullfile(thisDir,'output','surface_identification_gap_law');
D=load(dataPath,'Data'); D=D.Data;
x=D.x_mm(:); gap=D.gap_mm(:); tilt=D.tilt_deg(:); V=D.response_V;
[nX,nGap,nTilt]=size(V); X=[ones(nGap,1),1./gap,log(gap/0.8)];

% Representative open-set case: hidden 2 deg surface, anchor at 0.9 mm.
it=find(tilt==2,1); ia=find(abs(gap-0.9)<1e-12,1);
trainTilt=setdiff(1:nTilt,it); Q=zeros(3*nX,numel(trainTilt));
for k=1:numel(trainTilt)
    B=(X \ V(:,:,trainTilt(k)).').'; Q(:,k)=B(:);
end
mu=mean(Q,2); [U,~,~]=svd(Q-mu,'econ'); U=U(:,1:min(2,size(U,2)));
Z=U'*(Q-mu); [z1,ord]=sort(Z(1,:)); Qord=Q(:,ord);
dz=max(diff(z1)); tgrid=linspace(z1(1)-0.5*dz,z1(end)+0.5*dz,401);
y=V(:,ia,it); anchorPred=zeros(nX,numel(tgrid));
for k=1:numel(tgrid)
    q=interp1(z1,Qord.',tgrid(k),'linear','extrap').';
    anchorPred(:,k)=reshape(q,nX,3)*X(ia,:).';
end
d=sqrt(mean((anchorPred-y).^2,1)); [~,ib]=min(d); zhat=tgrid(ib);
qhat=interp1(z1,Qord.',zhat,'linear','extrap').'; Bhat=reshape(qhat,nX,3);
Vpred=Bhat*X.'; Vtrue=V(:,:,it);

% Error matrix for all target tilts and anchor gaps, regenerated consistently.
errMat=nan(nTilt,nGap); peakMat=nan(nTilt,nGap);
for jt=1:nTilt
    tr=setdiff(1:nTilt,jt); QQ=zeros(3*nX,numel(tr));
    for k=1:numel(tr), BB=(X \ V(:,:,tr(k)).').'; QQ(:,k)=BB(:); end
    mm=mean(QQ,2); [UU,~,~]=svd(QQ-mm,'econ'); UU=UU(:,1:min(2,size(UU,2)));
    zz=UU'*(QQ-mm); [zz1,oo]=sort(zz(1,:)); QQ=QQ(:,oo); dd=max(diff(zz1));
    tg=linspace(zz1(1)-0.5*dd,zz1(end)+0.5*dd,401);
    for ja=1:nGap
        pa=zeros(nX,numel(tg));
        for k=1:numel(tg)
            qq=interp1(zz1,QQ.',tg(k),'linear','extrap').'; pa(:,k)=reshape(qq,nX,3)*X(ja,:).';
        end
        [~,kk]=min(sqrt(mean((pa-V(:,ja,jt)).^2,1)));
        qq=interp1(zz1,QQ.',tg(kk),'linear','extrap').'; pp=reshape(qq,nX,3)*X.';
        eo=pp(:,setdiff(1:nGap,ja))-V(:,setdiff(1:nGap,ja),jt);
        truthOther=V(:,setdiff(1:nGap,ja),jt);
        errMat(jt,ja)=100*sqrt(mean(eo(:).^2))/range(truthOther(:));
        peakMat(jt,ja)=mean(abs(sg_peak_x(x,pp(:,setdiff(1:nGap,ja)))-sg_peak_x(x,V(:,setdiff(1:nGap,ja),jt))));
    end
end

fig=figure('Name','Predicted surface family visualization','Color','w', ...
    'Units','centimeters','Position',[2 2 18 14]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
% (a) Anchor localization.
nexttile; hold on;
plot(x,y,'k','LineWidth',1.5,'DisplayName','Hidden target anchor');
plot(x,anchorPred(:,ib),'Color',[0.85 0.20 0.08],'LineWidth',1.25,'DisplayName','Best predicted surface');
plot(x,anchorPred(:,max(1,ib-40)),'--','Color',[0.45 0.45 0.45],'LineWidth',0.8,'DisplayName','Neighbour candidate');
xlabel('Circumferential coordinate, x (mm)'); ylabel('Response, V (V)');
title('(a) Anchor-waveform localization'); legend('Location','best','Box','off'); style_axes(gca);
% (b) True and predicted surfaces.
nexttile; imagesc(gap,x,Vpred); axis xy tight; hold on; plot(gap(ia),x(round(end/2)),'wo','MarkerFaceColor','k','MarkerSize',4);
xlabel('Clearance, g (mm)'); ylabel('x (mm)'); title('(b) Predicted surface (hidden tilt = 2.0 deg)'); cb=colorbar; cb.Label.String='V (V)'; style_axes(gca);
% (c) Other-gap waveform transfer.
nexttile; hold on; others=setdiff(1:nGap,ia); cols=lines(numel(others));
for k=1:numel(others)
    j=others(k); plot(x,Vtrue(:,j),'Color',cols(k,:),'LineWidth',1.2,'DisplayName',sprintf('Truth, g=%.1f',gap(j)));
    plot(x,Vpred(:,j),'--','Color',cols(k,:),'LineWidth',0.9,'HandleVisibility','off');
end
xlabel('Circumferential coordinate, x (mm)'); ylabel('Response, V (V)'); title('(c) Other-gap transfer'); legend('Location','best','Box','off'); style_axes(gca);
% (d) Error across all hidden tilts and anchors.
nexttile; imagesc(gap,tilt,errMat); axis xy tight; xlabel('Anchor gap, g_a (mm)'); ylabel('Hidden tilt (deg)');
title('(d) Transfer NRMSE (%)'); cb=colorbar; cb.Label.String='NRMSE (%)'; style_axes(gca);

exportgraphics(fig,fullfile(outDir,'Fig_PredictedSurfaceFamily_Visualization.png'),'Resolution',300);
exportgraphics(fig,fullfile(outDir,'Fig_PredictedSurfaceFamily_Visualization.pdf'),'ContentType','vector');
fprintf('Representative target: tilt %.1f deg, anchor gap %.1f mm; best anchor RMSE %.6g V.\n', ...
    tilt(it),gap(ia),d(ib));
fprintf('Figures written to:\n%s\n',outDir);

function style_axes(ax)
set(ax,'FontName','Times New Roman','FontSize',8,'TickDir','in','Box','on','LineWidth',0.7);
end
function xp=sg_peak_x(x,Y)
n=size(Y,2); xp=zeros(1,n);
for q=1:n, z=sgolayfilt(Y(:,q),3,11); [~,ii]=max(z); xp(q)=x(ii); end
end
