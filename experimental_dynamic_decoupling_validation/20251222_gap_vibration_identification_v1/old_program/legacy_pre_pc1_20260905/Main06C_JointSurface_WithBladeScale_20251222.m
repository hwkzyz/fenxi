%% Step06C: frozen paths plus static blade response-scale correction
% Static-only model: Y_b = a_b F(g_b(x),x) + b_b.
% The blade scale/intercept are removed before the common surface is fitted;
% neither parameter is released in the downstream high-speed identification.

clear; clc; close all;
cfg = struct('referenceBladeId',2,'deltaBounds',[-0.25 0.25], ...
    'qBounds',[-0.06 0.06],'aBounds',[0.30 3.00],'bBounds',[-1500 1500], ...
    'minStates',8,'minPoints',300,'saveFigures',true);
pc = Setup_Paths_20251222(); outDir = pc.paths.calibrationWork;
step05File = fullfile(outDir,'Step05_Response_Surface_20251222_RefBladeAnchor.mat');
if ~isfile(step05File), error('Missing Step05 B2 anchor output.'); end
S=load(step05File,'responseSurface'); base=S.responseSurface;
x=base.xGrid(:); Y=base.waveforms; blades=base.bladeIds(:).'; gaps=base.trueGapMm(:).';
ibRef=find(blades==cfg.referenceBladeId,1); fitMask=base.effectiveWindow(:)&all(isfinite(base.coeff),2);
glo=min(base.gTrainMm); ghi=max(base.gTrainMm); nB=numel(blades); nG=numel(gaps);
delta=zeros(1,nB); q=zeros(1,nB); a=ones(1,nB); b=zeros(1,nB); rmse=nan(1,nB); nstate=zeros(1,nB);
for ib=1:nB
    if ib==ibRef
        [rmse(ib),nstate(ib)]=objective_local([0 0 1 0],Y(:,ib,:),x,gaps,base,fitMask,glo,ghi,cfg);
        continue
    end
    fun=@(p)objective_local(p,Y(:,ib,:),x,gaps,base,fitMask,glo,ghi,cfg);
    opt=optimset('Display','off','MaxIter',800,'MaxFunEvals',1800,'TolX',1e-6,'TolFun',1e-3);
    p=fminsearch(@(z)bounded_local(z,fun,cfg),[0 0 1 0],opt);
    p(1)=min(max(p(1),cfg.deltaBounds(1)),cfg.deltaBounds(2)); p(2)=min(max(p(2),cfg.qBounds(1)),cfg.qBounds(2));
    p(3)=min(max(p(3),cfg.aBounds(1)),cfg.aBounds(2)); p(4)=min(max(p(4),cfg.bBounds(1)),cfg.bBounds(2));
    [rmse(ib),nstate(ib)]=fun(p); delta(ib)=p(1); q(ib)=p(2); a(ib)=p(3); b(ib)=p(4);
end

path=nan(numel(x),nB,nG); valid=false(size(path)); Yn=nan(size(Y)); Ypred=nan(size(Y));
for ib=1:nB
    for ig=1:nG
        path(:,ib,ig)=gaps(ig)+delta(ib)+q(ib).*x;
        valid(:,ib,ig)=fitMask & path(:,ib,ig)>=glo & path(:,ib,ig)<=ghi & isfinite(Y(:,ib,ig));
        f=evalF_local(base.coeff,base.g0Mm,path(:,ib,ig));
        Ypred(:,ib,ig)=a(ib).*f+b(ib);
        Yn(:,ib,ig)=(Y(:,ib,ig)-b(ib))./a(ib);
    end
end

% Joint regression on normalized all-blade samples.
coeff=nan(numel(x),3); Yfit=nan(size(Yn)); count=zeros(numel(x),1); condA=nan(numel(x),1);
for ix=1:numel(x)
    if ~fitMask(ix), continue; end
    gg=squeeze(path(ix,:,:)); yy=squeeze(Yn(ix,:,:)); vv=squeeze(valid(ix,:,:)); gg=gg(vv); yy=yy(vv);
    if numel(gg)<12 || rank([ones(numel(gg),1),1./gg,log(gg./base.g0Mm)])<3, continue; end
    A=[ones(numel(gg),1),1./gg,log(gg./base.g0Mm)]; coeff(ix,:)=(A\yy).'; count(ix)=numel(gg); condA(ix)=cond(A);
    for ib=1:nB, for ig=1:nG, if valid(ix,ib,ig), gg0=path(ix,ib,ig); Yfit(ix,ib,ig)=coeff(ix,1)+coeff(ix,2)/gg0+coeff(ix,3)*log(gg0/base.g0Mm); end, end, end
end

rows=cell(nB*nG,1); k=0;
for ib=1:nB, for ig=1:nG, k=k+1; vv=fitMask&valid(:,ib,ig)&isfinite(Yfit(:,ib,ig)); r=Yn(vv,ib,ig)-Yfit(vv,ib,ig); rows{k}=table(blades(ib),gaps(ig),a(ib),b(ib),delta(ib),q(ib),nnz(vv),sqrt(mean(r.^2,'omitnan')),'VariableNames',{'bladeId','nominalGapMm','aBlade','bBladeMv','deltaGapMm','qGapPerXMm','pointCount','rmseMv'}); end, end
fitTable=vertcat(rows{:}); bladeSummary=groupsummary(fitTable,'bladeId',{'mean','std'},{'rmseMv'});
pathTable=table(blades(:),delta(:),q(:),a(:),b(:),atan(q(:))*180/pi,rmse(:),nstate(:),'VariableNames',{'bladeId','deltaGapMm','qGapPerXMm','aBlade','bBladeMv','tiltAngleDeg','scaledPathRmseMv','validGapStateCount'});

JointResponseSurface=base; JointResponseSurface.method='B2_anchor_frozen_paths_with_static_blade_scale';
JointResponseSurface.description='Static blade scale/intercept are fitted and removed; high-speed model uses only frozen F(g,x) increments.';
JointResponseSurface.pathTable=pathTable; JointResponseSurface.gapPathMm=path; JointResponseSurface.pathValid=valid;
JointResponseSurface.bladeNormalizedWaveforms=Yn; JointResponseSurface.coeff=coeff; JointResponseSurface.YfitAllBlades=Yfit;
JointResponseSurface.fitTable=fitTable; JointResponseSurface.bladeFitSummary=bladeSummary; JointResponseSurface.sampleCountByX=count; JointResponseSurface.conditionNumberByX=condA;
matFile=fullfile(outDir,'Step06C_Joint_Response_Surface_WithBladeScale_20251222.mat');
save(matFile,'JointResponseSurface','-v7.3'); writetable(pathTable,fullfile(outDir,'Step06C_BladeScale_PathSummary_20251222.csv')); writetable(fitTable,fullfile(outDir,'Step06C_JointSurface_Fit_20251222.csv')); writetable(bladeSummary,fullfile(outDir,'Step06C_JointSurface_BladeSummary_20251222.csv'));

if cfg.saveFigures
    f=figure('Color','w','Position',[80 80 1300 760]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    nexttile; hold on; grid on; for ib=1:nB, plot(x,path(:,ib,round(nG/2)),'LineWidth',1.2,'DisplayName',sprintf('B%d',blades(ib))); end; yline(glo,'k--'); yline(ghi,'k--'); xlabel('x (mm)'); ylabel('g path (mm)'); title('Frozen paths'); legend('Location','best');
    nexttile; yyaxis left; bar(blades,a); ylabel('a_b'); yyaxis right; plot(blades,b,'ko-'); ylabel('b_b (mV)'); xlabel('Blade ID'); title('Static blade response mapping');
    nexttile; hold on; grid on; for ib=1:nB, z=fitTable.bladeId==blades(ib); plot(fitTable.nominalGapMm(z),fitTable.rmseMv(z),'o-','DisplayName',sprintf('B%d',blades(ib))); end; xlabel('Nominal gap (mm)'); ylabel('Normalized joint RMSE (mV)'); legend('Location','best');
    nexttile; semilogy(x,condA); grid on; xlabel('x (mm)'); ylabel('cond(A)'); title('Joint regression conditioning');
    saveas(f,fullfile(outDir,'Step06C_JointSurface_WithBladeScale_Diagnostics_20251222.png'));
end
disp(pathTable); fprintf('Step06C saved: %s\n',matFile);

function v=bounded_local(p,fun,cfg)
if p(1)<cfg.deltaBounds(1)||p(1)>cfg.deltaBounds(2)||p(2)<cfg.qBounds(1)||p(2)>cfg.qBounds(2)||p(3)<cfg.aBounds(1)||p(3)>cfg.aBounds(2)||p(4)<cfg.bBounds(1)||p(4)>cfg.bBounds(2), v=1e10; else, v=fun(p); end
end
function [rmse,nstate]=objective_local(p,Yraw,x,gaps,base,fitMask,glo,ghi,cfg)
Y=squeeze(Yraw); rr=[]; nstate=0;
for ig=1:numel(gaps), g=gaps(ig)+p(1)+p(2).*x; vv=fitMask&g>=glo&g<=ghi&isfinite(Y(:,ig)); if nnz(vv)>=cfg.minPoints/numel(gaps), pred=p(3).*evalF_local(base.coeff,base.g0Mm,g)+p(4); rr=[rr;Y(vv,ig)-pred(vv)]; nstate=nstate+1; end, end
if nstate<cfg.minStates||numel(rr)<cfg.minPoints, rmse=1e8; else, rmse=sqrt(mean(rr.^2)); end
end
function f=evalF_local(c,g0,g), f=nan(size(g)); vv=isfinite(g)&all(isfinite(c),2); f(vv)=c(vv,1)+c(vv,2)./g(vv)+c(vv,3).*log(g(vv)./g0); end
