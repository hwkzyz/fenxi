%% End-to-end visualization: gap library -> surface family -> localization -> transfer.
% The target tilt is hidden from the algorithm.  This figure visualizes the
% data and operations rather than only the final prediction error.

thisDir=fileparts(mfilename('fullpath'));
dataPath=fullfile(thisDir,'output','waveform_families','TiltGap_WaveformSurfaceData.mat');
outDir=fullfile(thisDir,'output','surface_identification_gap_law');
if ~exist(outDir,'dir'), mkdir(outDir); end
D=load(dataPath,'Data'); D=D.Data;
x=D.x_mm(:); gap=D.gap_mm(:); tilt=D.tilt_deg(:); V=D.response_V;
[nX,nGap,nTilt]=size(V); X=[ones(nGap,1),1./gap,log(gap/0.8)];

% Representative leave-one-surface-out case.
it=find(tilt==2,1); ia=find(abs(gap-0.9)<1e-12,1); trainTilt=setdiff(1:nTilt,it);
Q=zeros(3*nX,numel(trainTilt)); Btrain=cell(numel(trainTilt),1);
for k=1:numel(trainTilt)
    Btrain{k}=(X \ V(:,:,trainTilt(k)).').'; Q(:,k)=Btrain{k}(:);
end
mu=mean(Q,2); [U,Sv,~]=svd(Q-mu,'econ'); r=min(2,size(U,2)); Ur=U(:,1:r); Z=Ur'*(Q-mu);
[z1,ord]=sort(Z(1,:)); Qord=Q(:,ord); trainTiltOrd=trainTilt(ord);
dz=max(diff(z1)); zgrid=linspace(z1(1)-0.5*dz,z1(end)+0.5*dz,401);
y=V(:,ia,it); anchorPred=zeros(nX,numel(zgrid));
for k=1:numel(zgrid)
    q=interp1(z1,Qord.',zgrid(k),'linear','extrap').'; anchorPred(:,k)=reshape(q,nX,3)*X(ia,:).';
end
dist=sqrt(mean((anchorPred-y).^2,1)); [~,ib]=min(dist); zhat=zgrid(ib);
qhat=interp1(z1,Qord.',zhat,'linear','extrap').'; Bhat=reshape(qhat,nX,3); Vhat=Bhat*X.';

% Fit quality of the within-tilt gap law for one representative training face.
itr=find(tilt==1.5,1); Bdemo=(X \ V(:,:,itr).').'; Vdemo=Bdemo*X.';

fig=figure('Name','End-to-end surface-transfer story','Color','w', ...
    'Units','centimeters','Position',[2 2 20 20]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

% (a) Discrete gap waveform library.
nexttile; hold on; cmap=parula(nGap);
for j=1:nGap, plot(x,Vdemo(:,j),'Color',cmap(j,:),'LineWidth',1.0); end
xlabel('x (mm)'); ylabel('V (V)'); title('(a) Gap waveform library');
text(0.03,0.94,sprintf('one calibration face, tilt = %.1f deg',tilt(itr)), ...
    'Units','normalized','FontSize',7,'Color',[0.2 0.2 0.2]); style_axes(gca);

% (b) Within-face gap law reconstruction.
nexttile; hold on;
for j=1:nGap
    plot(x,Vdemo(:,j),'Color',cmap(j,:),'LineWidth',0.8);
    plot(x,V(:,j,itr),'k:','LineWidth',0.45,'HandleVisibility','off');
end
xlabel('x (mm)'); ylabel('V (V)'); title('(b) Fitted within-face gap law');
text(0.03,0.94,'solid: fitted law; dotted: raw COMSOL', ...
    'Units','normalized','FontSize',7,'Color',[0.2 0.2 0.2]); style_axes(gca);

% (c) A fitted face as an x-by-gap surface.
nexttile; imagesc(gap,x,Vdemo); axis xy tight; xlabel('g (mm)'); ylabel('x (mm)');
title('(c) One parameterized waveform face'); cb=colorbar; cb.Label.String='V (V)'; style_axes(gca);

% (d) Surface-family latent coordinates and hidden target.
nexttile; hold on; scatter(Z(1,:),Z(min(2,r),:),38,tilt(trainTilt),'filled');
qtarget=(X \ V(:,:,it).').'; qtarget=qtarget(:); ztarget=Ur'*(qtarget-mu);
plot(ztarget(1),ztarget(min(2,r)),'kp','MarkerSize',10,'MarkerFaceColor','w','LineWidth',1.1);
xlabel('Latent coordinate z_1'); ylabel('z_2'); title('(d) SVD surface manifold');
text(0.04,0.07,'star: hidden target surface', 'Units','normalized','FontSize',7); style_axes(gca);

% (e) Blind localization uses only the anchor waveform mismatch.
nexttile; hold on; plot(zgrid,dist,'k','LineWidth',1.1); plot(zhat,dist(ib),'ro','MarkerFaceColor','r','MarkerSize',4);
xlabel('Latent coordinate z'); ylabel('Anchor mismatch RMSE (V)'); title('(e) Anchor localizes the face');
text(0.04,0.92,sprintf('anchor: g_a = %.1f mm',gap(ia)), ...
    'Units','normalized','FontSize',7,'Color',[0.2 0.2 0.2]);
text(0.04,0.84,sprintf('selected z = %.3f',zhat), ...
    'Units','normalized','FontSize',7,'Color',[0.2 0.2 0.2]); style_axes(gca);

% (f) Predicted versus true hidden surface.
nexttile; imagesc(gap,x,Vhat); axis xy tight; xlabel('g (mm)'); ylabel('x (mm)');
title('(f) Predicted hidden surface'); cb=colorbar; cb.Label.String='V (V)'; style_axes(gca);

% (g) Other-gap waveform transfer; anchor is intentionally omitted.
nexttile; hold on; others=setdiff(1:nGap,ia); cols=lines(numel(others));
for k=1:numel(others)
    j=others(k); plot(x,V(:,j,it),'Color',cols(k,:),'LineWidth',1.05,'DisplayName',sprintf('g=%.1f mm',gap(j)));
    plot(x,Vhat(:,j),'--','Color',cols(k,:),'LineWidth',0.8,'HandleVisibility','off');
end
xlabel('x (mm)'); ylabel('V (V)'); title('(g) Other-gap transfer');
legend('Location','southoutside','Orientation','horizontal','Box','off','FontSize',6); style_axes(gca);

% (h) Full protocol error map.
errMap=nan(nTilt,nGap);
for jt=1:nTilt
    tr=setdiff(1:nTilt,jt); QQ=zeros(3*nX,numel(tr));
    for k=1:numel(tr), BB=(X \ V(:,:,tr(k)).').'; QQ(:,k)=BB(:); end
    mm=mean(QQ,2); [UU,~,~]=svd(QQ-mm,'econ'); UU=UU(:,1:min(2,size(UU,2)));
    zz=UU'*(QQ-mm); [zz1,oo]=sort(zz(1,:)); QQ=QQ(:,oo); dd=max(diff(zz1)); tg=linspace(zz1(1)-0.5*dd,zz1(end)+0.5*dd,401);
    for ja=1:nGap
        pa=zeros(nX,numel(tg));
        for k=1:numel(tg), qq=interp1(zz1,QQ.',tg(k),'linear','extrap').'; pa(:,k)=reshape(qq,nX,3)*X(ja,:).'; end
        [~,kk]=min(sqrt(mean((pa-V(:,ja,jt)).^2,1))); qq=interp1(zz1,QQ.',tg(kk),'linear','extrap').'; pp=reshape(qq,nX,3)*X.';
        oth=setdiff(1:nGap,ja); e=pp(:,oth)-V(:,oth,jt); errMap(jt,ja)=100*sqrt(mean(e(:).^2))/range(V(:,oth,jt), 'all');
    end
end
nexttile; imagesc(gap,tilt,errMap); axis xy tight; xlabel('Anchor gap g_a (mm)'); ylabel('Hidden tilt (deg)');
title('(h) Transfer error across protocol'); cb=colorbar; cb.Label.String='NRMSE (%)'; style_axes(gca);

% (i) Protocol summary: tilt is not supplied to the localization step.
nexttile; axis off;
text(0.02,0.90,'(i) Blind surface-transfer protocol','FontWeight','bold','FontSize',8);
text(0.04,0.70,'1. Fit the within-face gap law','FontSize',7);
text(0.04,0.53,'2. Build a low-rank surface family','FontSize',7);
text(0.04,0.36,'3. Localize with one anchor waveform','FontSize',7);
text(0.04,0.19,'4. Transfer to the remaining gaps','FontSize',7);
text(0.04,0.04,'Tilt labels are used only for evaluation.', ...
    'FontSize',6.5,'Color',[0.25 0.25 0.25]);

exportgraphics(fig,fullfile(outDir,'Fig_EndToEnd_SurfaceTransfer_Story.png'),'Resolution',300);
exportgraphics(fig,fullfile(outDir,'Fig_EndToEnd_SurfaceTransfer_Story.pdf'),'ContentType','vector');
fprintf('Saved end-to-end figure for hidden tilt %.1f deg, anchor gap %.1f mm.\n%s\n', ...
    tilt(it),gap(ia),outDir);

function style_axes(ax)
set(ax,'FontName','Times New Roman','FontSize',7.5,'TickDir','in','Box','on','LineWidth',0.65);
end
