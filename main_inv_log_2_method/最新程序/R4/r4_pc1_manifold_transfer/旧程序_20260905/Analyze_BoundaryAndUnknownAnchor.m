%% Boundary-surface extrapolation and unknown-anchor identifiability.
% Ideal noiseless simulation: tilt labels are never supplied to localization.

thisDir=fileparts(mfilename('fullpath'));
dataPath=fullfile(thisDir,'output','waveform_families','TiltGap_WaveformSurfaceData.mat');
outDir=fullfile(thisDir,'output','boundary_unknown_anchor');
if ~exist(outDir,'dir'), mkdir(outDir); end
D=load(dataPath,'Data'); D=D.Data;
x=D.x_mm(:); gap=D.gap_mm(:); tilt=D.tilt_deg(:); V=D.response_V;
[nX,nGap,nTilt]=size(V); gref=0.8;
X=[ones(nGap,1),1./gap,log(gap/gref)];
gSearch=linspace(min(gap),max(gap),401).';

% Two endpoint faces are hidden; each endpoint is tested with every anchor gap.
targetFaces=[1 nTilt];
Known=struct('target',{},'anchor',{},'rmse_anchor',{},'nrmse_other',{}, ...
    'corr_other',{},'peak_x_mae',{},'peak_amp_mae',{},'z_error',{},'g_hat',{},'g_error',{});
Joint=Known; caseNo=0; jointNo=0;
for target=targetFaces
    train=setdiff(1:nTilt,target);
    [z1,Qord,mu,Ur]=fit_family(V,X,train);
    dz=max(diff(z1));
    zgrid=linspace(z1(1)-1.25*dz,z1(end)+1.25*dz,601);
    for ia=1:nGap
        % Candidate family evaluated at all gaps for this latent coordinate.
        Cand=build_candidates(zgrid,Qord,z1,nX,X);
        yAnchor=V(:,ia,target);
        anchorPred=squeeze(Cand(:,ia,:));
        dist=sqrt(mean((anchorPred-yAnchor).^2,1));
        [~,ib]=min(dist); zhat=zgrid(ib);
        qhat=interp1(z1,Qord.',zhat,'linear','extrap').'; Vhat=reshape(qhat,nX,3)*X.';
        oth=setdiff(1:nGap,ia); E=Vhat(:,oth)-V(:,oth,target);
        Toth=V(:,oth,target); Poth=Vhat(:,oth);
        [px,pa]=peak_metrics(V(:,oth,target),Vhat(:,oth),x);
        caseNo=caseNo+1;
        Known(caseNo)=struct('target',target,'anchor',ia,'rmse_anchor',dist(ib), ...
            'nrmse_other',100*sqrt(mean(E(:).^2))/range(Toth,'all'), ...
            'corr_other',corr(Toth(:),Poth(:)), ...
            'peak_x_mae',px,'peak_amp_mae',pa,'z_error',NaN,'g_hat',gap(ia),'g_error',0);
        % Joint search: anchor gap is not supplied to the locator.
        J=zeros(numel(zgrid),numel(gSearch));
        for iz=1:numel(zgrid)
            B=reshape(interp1(z1,Qord.',zgrid(iz),'linear','extrap').',nX,3);
            P=B*[ones(1,numel(gSearch));1./gSearch.';log(gSearch.'/gref)];
            J(iz,:)=sqrt(mean((P-yAnchor).^2,1));
        end
        [~,lin]=min(J(:)); [izg,igg]=ind2sub(size(J),lin);
        zjoint=zgrid(izg); gjoint=gSearch(igg);
        qjoint=interp1(z1,Qord.',zjoint,'linear','extrap').'; Bj=reshape(qjoint,nX,3);
        Vjoint=Bj*X.'; Ej=Vjoint(:,oth)-V(:,oth,target);
        Pjoth=Vjoint(:,oth);
        [pxj,paj]=peak_metrics(V(:,oth,target),Vjoint(:,oth),x);
        jointNo=jointNo+1;
        Joint(jointNo)=struct('target',target,'anchor',ia,'rmse_anchor',J(izg,igg), ...
            'nrmse_other',100*sqrt(mean(Ej(:).^2))/range(Toth,'all'), ...
            'corr_other',corr(Toth(:),Pjoth(:)), ...
            'peak_x_mae',pxj,'peak_amp_mae',paj,'z_error',zjoint,'g_hat',gjoint,'g_error',gjoint-gap(ia));
        if ia==1
            save_case_figure(outDir,x,gap,tilt,target,ia,zgrid,dist,J,gSearch, ...
                V(:,:,target),Vhat,Vjoint,Joint(jointNo));
        end
    end
end

Summary=struct(); Summary.description='Ideal boundary-face extrapolation and unknown-anchor joint search';
Summary.target_faces=targetFaces; Summary.tilt_deg=tilt; Summary.gap_mm=gap;
Summary.known_anchor=Known; Summary.unknown_anchor=Joint;
Summary.known_anchor_mean_nrmse=mean([Known.nrmse_other]);
Summary.unknown_anchor_mean_nrmse=mean([Joint.nrmse_other]);
save(fullfile(outDir,'Boundary_UnknownAnchor_Results.mat'),'Summary','-v7.3');
fprintf('Boundary/unknown-anchor analysis complete: %d known-anchor and %d joint-search cases.\n',caseNo,jointNo);
fprintf('Known-anchor mean other-gap NRMSE = %.4f%%; unknown-anchor mean = %.4f%%.\n', ...
    Summary.known_anchor_mean_nrmse,Summary.unknown_anchor_mean_nrmse);

function [z1,Qord,mu,Ur]=fit_family(V,X,train)
nX=size(V,1); Q=zeros(3*nX,numel(train));
for k=1:numel(train), B=(X\V(:,:,train(k)).').'; Q(:,k)=B(:); end
mu=mean(Q,2); [U,~,~]=svd(Q-mu,'econ'); Ur=U(:,1:min(2,size(U,2))); Z=Ur'*(Q-mu);
[z1,ord]=sort(Z(1,:)); Qord=Q(:,ord);
end

function Cand=build_candidates(zgrid,Qord,z1,nX,X)
Cand=zeros(nX,size(X,1),numel(zgrid));
for k=1:numel(zgrid)
q=interp1(z1,Qord.',zgrid(k),'linear','extrap').'; Cand(:,:,k)=reshape(q,nX,3)*X.';
end
end

function [xMae,aMae]=peak_metrics(T,P,x)
n=size(T,2); xt=zeros(n,1); xp=zeros(n,1); at=zeros(n,1); ap=zeros(n,1);
for j=1:n
 [xt(j),at(j)]=sg_peak(x,T(:,j)); [xp(j),ap(j)]=sg_peak(x,P(:,j));
end
xMae=mean(abs(xp-xt)); aMae=mean(abs(ap-at));
end

function [xp,yp]=sg_peak(x,y)
if exist('sgolayfilt','file')==2
    ys=sgolayfilt(y,3,min(11,2*floor((numel(y)-1)/2)+1));
else
    ys=y;
end
[~,k]=max(ys); k=max(2,min(numel(y)-1,k));
xx=x(k-1:k+1); yy=ys(k-1:k+1); c=polyfit(xx,yy,2);
xp=-c(2)/(2*c(1)); yp=polyval(c,xp);
end

function save_case_figure(outDir,x,gap,tilt,target,ia,zgrid,dist,J,gSearch,Vtrue,Vknown,Vjoint,joint)
fig=figure('Name','Boundary and unknown-anchor validation','Color','w','Units','centimeters','Position',[2 2 21 14]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile; imagesc(gap,x,Vtrue); axis xy tight; xlabel('g (mm)'); ylabel('x (mm)'); title('True hidden boundary surface'); colorbar; style_axes(gca);
nexttile; imagesc(gap,x,Vknown); axis xy tight; xlabel('g (mm)'); ylabel('x (mm)'); title('Known-anchor prediction'); colorbar; style_axes(gca);
nexttile; imagesc(gap,x,Vjoint); axis xy tight; xlabel('g (mm)'); ylabel('x (mm)'); title('Unknown-anchor prediction'); colorbar; style_axes(gca);
nexttile; plot(zgrid,dist,'k','LineWidth',1.2); xlabel('Latent coordinate z'); ylabel('Anchor RMSE (V)'); title(sprintf('Known g_a = %.1f mm',gap(ia))); style_axes(gca);
nexttile; imagesc(gSearch,zgrid,J); axis xy tight; hold on;
plot(joint.g_hat,joint.z_error,'wp','MarkerFaceColor','k','MarkerSize',8);
xlabel('Candidate g (mm)'); ylabel('Latent coordinate z'); title('Joint cost J(z,g)');
text(0.04,0.92,sprintf('minimum: g = %.3f mm',joint.g_hat),'Units','normalized','FontSize',7,'Color','w');
colorbar; style_axes(gca);
nexttile; hold on; oth=setdiff(1:numel(gap),ia); cols=lines(numel(oth));
for k=1:numel(oth), j=oth(k); plot(x,Vtrue(:,j),'Color',cols(k,:),'LineWidth',1); plot(x,Vjoint(:,j),'--','Color',cols(k,:),'LineWidth',0.8); end
xlabel('x (mm)'); ylabel('V (V)'); title(sprintf('Other-gap transfer, hidden tilt %.1f deg',tilt(target))); style_axes(gca);
exportgraphics(fig,fullfile(outDir,sprintf('BoundaryCase_tilt%.1f_anchor%.1f.png',tilt(target),gap(ia))),'Resolution',300);
exportgraphics(fig,fullfile(outDir,sprintf('BoundaryCase_tilt%.1f_anchor%.1f.pdf',tilt(target),gap(ia))),'ContentType','vector');
end

function style_axes(ax), set(ax,'FontName','Times New Roman','FontSize',8,'TickDir','in','Box','on','LineWidth',0.65); end
