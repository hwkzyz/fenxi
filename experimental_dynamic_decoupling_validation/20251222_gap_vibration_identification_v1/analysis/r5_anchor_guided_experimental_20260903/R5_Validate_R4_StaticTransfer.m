function R = R5_Validate_R4_StaticTransfer(familyFile, outputDir)
% Independent static validation of the R4 waveform-surface transfer idea.
% It does not alter any Main05/Main07/Main10 production artifact.
if nargin < 1 || isempty(familyFile)
    familyFile=fullfile(fileparts(mfilename('fullpath')),'results', ...
        'r4_frontend_only','r4_experimental_surface_family.mat');
end
if nargin < 2 || isempty(outputDir)
    outputDir=fullfile(fileparts(mfilename('fullpath')),'results', ...
        'r4_frontend_only','static_validation');
end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
load(familyFile,'F');
x=F.x_mm(:); g=F.gap_mm(:); Y=F.waveforms_mv; support=F.support_mask(:);
nB=numel(F.blade_id); nG=numel(g);

% Same-blade / leave-gap: a held gap is absent from every gap-law fit;
% another measured gap of that blade is the single localization anchor.
sameRows={}; sameDetails={}; ks=0;
for ih=1:nG
    keep=setdiff(1:nG,ih); ia=nearest_index_local(g(keep),F.gref_mm); ia=keep(ia);
    M=build_family_local(Y,g,F.gref_mm,keep,1:nB);
    for ib=1:nB
        [z,B]=localize_local(M,Y(:,ib,ia),g(ia),support,F.gref_mm);
        yp=evaluate_local(B,g(ih),F.gref_mm);
        yt=Y(:,ib,ih);
        ks=ks+1; sameDetails{ks,1}=detail_local(ib,ih,ia,z,x,yt,yp,support);
        sameRows{ks,1}=row_local(sameDetails{ks},'same_blade_leave_gap');
    end
end
same=vertcat(sameRows{:});

% Leave-blade: no waveform of the held blade is used to fit its surface.
% The reference-gap waveform is the only observed anchor for localization.
leaveRows={}; leaveDetails={}; kl=0; ia=nearest_index_local(g,F.gref_mm);
for ib=1:nB
    train=setdiff(1:nB,ib);
    M=build_family_local(Y,g,F.gref_mm,1:nG,train);
    [z,B]=localize_local(M,Y(:,ib,ia),g(ia),support,F.gref_mm);
    for ig=setdiff(1:nG,ia)
        yp=evaluate_local(B,g(ig),F.gref_mm); yt=Y(:,ib,ig);
        kl=kl+1; leaveDetails{kl,1}=detail_local(ib,ig,ia,z,x,yt,yp,support);
        leaveRows{kl,1}=row_local(leaveDetails{kl},'leave_blade');
    end
end
leave=vertcat(leaveRows{:});
writetable(same,fullfile(outputDir,'R5_same_blade_leave_gap.csv'));
writetable(leave,fullfile(outputDir,'R5_leave_blade.csv'));
plot_representative_local(same,sameDetails, ...
    fullfile(outputDir,'R5_same_blade_leave_gap_waveform.png'), ...
    'Same blade, held gap');
plot_representative_local(leave,leaveDetails, ...
    fullfile(outputDir,'R5_leave_blade_waveform.png'), ...
    'Held blade, reference-gap anchor');
R=struct('schema','R5_R4_STATIC_TRANSFER_VALIDATION_V1', ...
    'familyFile',familyFile,'same_blade_leave_gap',same, ...
    'leave_blade',leave,'support_point_count',nnz(support), ...
    'same_blade_mean_rmse_mv',mean(same.rmse_mv), ...
    'leave_blade_mean_rmse_mv',mean(leave.rmse_mv));
save(fullfile(outputDir,'R5_static_transfer_validation.mat'),'R','sameDetails','leaveDetails','-v7.3');
fprintf('R5 static transfer validation saved: %s\n',outputDir);
end

function M=build_family_local(Y,g,gref,keep,blades)
nX=size(Y,1); Q=nan(3*nX,numel(blades));
X=[ones(numel(keep),1),1./g(keep),log(g(keep)/gref)];
for k=1:numel(blades)
    B=(X\squeeze(Y(:,blades(k),keep)).').'; Q(:,k)=B(:);
end
mu=mean(Q,2); [U,~,~]=svd(Q-mu,'econ'); r=min(2,size(U,2));
Z=U(:,1:r).'*(Q-mu); [z,ord]=sort(Z(1,:));
M=struct('q',Q(:,ord),'z',z,'x_count',nX,'mean_q',mu);
end

function [zBest,Bbest]=localize_local(M,yAnchor,gAnchor,support,gref)
zGrid=linspace(M.z(1),M.z(end),301); best=inf; zBest=NaN; Bbest=[];
for iz=1:numel(zGrid)
    q=interp1(M.z,M.q.',zGrid(iz),'linear').'; B=reshape(q,M.x_count,3);
    yp=evaluate_local(B,gAnchor,gref); e=sqrt(mean((yp(support)-yAnchor(support)).^2));
    if e<best, best=e; zBest=zGrid(iz); Bbest=B; end
end
end

function y=evaluate_local(B,g,gref)
y=B*[1;1/g;log(g/gref)];
end

function k=nearest_index_local(v,target)
[~,k]=min(abs(v-target));
end

function d=detail_local(blade,held,anchor,z,x,yTrue,yPred,support)
e=yPred(support)-yTrue(support); scale=max(yTrue(support))-min(yTrue(support));
d=struct('blade_id',blade,'held_gap_index',held,'anchor_gap_index',anchor, ...
    'latent_coordinate',z,'x_mm',x,'true_mv',yTrue,'pred_mv',yPred, ...
    'support',support,'rmse_mv',sqrt(mean(e.^2)), ...
    'nrmse_percent',100*sqrt(mean(e.^2))/max(scale,eps));
end

function T=row_local(d,method)
T=table(string(method),d.blade_id,d.held_gap_index,d.anchor_gap_index, ...
    d.latent_coordinate,d.rmse_mv,d.nrmse_percent,nnz(d.support), ...
    'VariableNames',{'method','blade_id','held_gap_index','anchor_gap_index', ...
    'latent_coordinate','rmse_mv','nrmse_percent','support_point_count'});
end

function plot_representative_local(T,D,file,titleText)
[~,k]=min(abs(T.rmse_mv-median(T.rmse_mv))); d=D{k};
fig=figure('Visible','off','Color','w','Position',[100 100 820 420]);
plot(d.x_mm(d.support),d.true_mv(d.support),'k-','LineWidth',1.6); hold on;
plot(d.x_mm(d.support),d.pred_mv(d.support),'r--','LineWidth',1.6);
xlabel('x (mm)'); ylabel('Voltage (mV)'); grid on;
title(sprintf('%s: B%d, held gap index %d, RMSE %.2f mV', ...
    titleText,d.blade_id,d.held_gap_index,d.rmse_mv));
legend({'held-out measured waveform','R4 anchor-localized prediction'}, ...
    'Location','best');
exportgraphics(fig,file,'Resolution',180); close(fig);
end
