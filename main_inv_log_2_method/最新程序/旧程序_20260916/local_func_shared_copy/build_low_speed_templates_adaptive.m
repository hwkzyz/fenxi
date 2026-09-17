function low=build_low_speed_templates_adaptive(dataLow,alpha_k,R_tip,domain,xGrid,method,opts)
%BUILD_LOW_SPEED_TEMPLATES_ADAPTIVE Low-speed-only regularized template.
% Hyperparameters are selected by grouped revolution cross-validation.  A
% complete revolution is always assigned to one fold, so samples from the
% same passage can never occur in both training and validation data.
if nargin<7||isempty(opts),opts=struct();end
dx=get_opt(opts,'gridSpacingMm',.02);minCount=get_opt(opts,'minBinCount',5);
minCoverage=get_opt(opts,'minCoverage',.90);nFold=get_opt(opts,'numFolds',5);
mapRatio=get_opt(opts,'mapWindowRatio',.48);
mapped=map_highspeed_to_space(dataLow,alpha_k,R_tip,domain,0,mapRatio);
xBin=(domain(1):dx:domain(2)).';if xBin(end)<domain(2)-eps,xBin(end+1)=domain(2);end
xOut=xGrid(:);sensorIds=unique(mapped.S_v(:),'stable');revs=unique(mapped.rev_v(:));
foldId=mod((1:numel(revs))-1,nFold)+1;
switch lower(char(method))
    case {'cv_sg','split_sg'}
        width=get_opt(opts,'candidateWindowMm',[.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
        span=max(5,2*floor((width/dx)/2)+1);width=span*dx;
        fitfun=@(m,s)fit_sg(m,s,xBin,minCount,minCoverage,span);
        [~,cv]=grouped_cv(mapped,revs,foldId,sensorIds,xBin,fitfun,numel(span));
        % Prediction-risk minimization is used for the forward template.
        % The deliberately smoother one-SE derivative is handled separately
        % by method D in the benchmark.
        [~,choice]=min(cv.mean);
        [Tb,Cb,Sb]=raw_bins(mapped,sensorIds,xBin,minCount,minCoverage);
        Tfit=apply_sg(Tb,span(choice));Dfit=gradient_matrix(Tfit,dx);
        low=package_low(mapped,xOut,xBin,Tfit,Dfit,Cb,Sb,sensorIds,method,dx);
        low.selected_span=span(choice);low.selected_window_mm=width(choice);
        low.cv_candidates=table(span(:),width(:),cv.mean(:),cv.se(:),cv.eligible(:),...
            'VariableNames',{'span_bins','window_mm','cv_rmse_V','cv_se_V','one_se_eligible'});
    case 'cv_spline'
        mult=get_opt(opts,'candidateToleranceMultiplier',[.001 .002 .005 .01 .02 .05 .1 .2]);
        fitfun=@(m,s)fit_spline(m,s,xBin,minCount,minCoverage,mult);
        [choice,cv]=grouped_cv(mapped,revs,foldId,sensorIds,xBin,fitfun,numel(mult));
        [Tb,Cb,Sb]=raw_bins(mapped,sensorIds,xBin,minCount,minCoverage);
        [Tfit,Dfit,pp,tol]=apply_spline(Tb,Cb,Sb,xBin,mult(choice));
        low=package_low(mapped,xOut,xBin,Tfit,Dfit,Cb,Sb,sensorIds,method,dx);
        low.selected_tolerance_multiplier=mult(choice);low.spline_pp=pp;low.spline_tolerance=tol;
        low.cv_candidates=table(mult(:),cv.mean(:),cv.se(:),cv.eligible(:),...
            'VariableNames',{'tolerance_multiplier','cv_rmse_V','cv_se_V','one_se_eligible'});
    otherwise,error('Unknown adaptive template method: %s',method);
end
if strcmpi(method,'cv_sg'),low.cv_rule='minimum grouped-revolution prediction error';else,low.cv_rule='largest regularization within one standard error of minimum grouped-CV error';end
low.cv_fold_count=nFold;low.cv_grouping='complete revolutions';
low.map_window_ratio=mapRatio;
end

function [choice,out]=grouped_cv(mapped,revs,foldId,sensorIds,xBin,fitfun,nCand)
loss=nan(max(foldId),numel(sensorIds),nCand);
for f=1:max(foldId)
    testRev=revs(foldId==f);train=~ismember(mapped.rev_v,testRev);test=~train;
    for is=1:numel(sensorIds)
        tr=train&mapped.S_v==sensorIds(is);te=test&mapped.S_v==sensorIds(is);
        fits=fitfun(submap(mapped,tr),1);
        for ic=1:nCand
            pred=interp1(xBin,fits(:,ic),mapped.x_v(te),'pchip','extrap');
            loss(f,is,ic)=sqrt(mean((mapped.V_a(te)-pred).^2,'omitnan'));
        end
    end
end
z=reshape(loss,[],nCand);mu=mean(z,1,'omitnan');se=std(z,0,1,'omitnan')./sqrt(sum(isfinite(z),1));
[~,imin]=min(mu);eligible=mu<=mu(imin)+se(imin);choice=find(eligible,1,'last');
out=struct('mean',mu,'se',se,'eligible',eligible,'raw',loss);
end

function fits=fit_sg(mapped,~,xBin,minCount,minCoverage,span)
[Tb,~,~]=raw_bins(mapped,unique(mapped.S_v(:),'stable'),xBin,minCount,minCoverage);
fits=nan(numel(xBin),numel(span));for i=1:numel(span),q=apply_sg(Tb,span(i));fits(:,i)=q(:,1);end
end

function fits=fit_spline(mapped,~,xBin,minCount,minCoverage,mult)
[Tb,Cb,Sb]=raw_bins(mapped,unique(mapped.S_v(:),'stable'),xBin,minCount,minCoverage);
fits=nan(numel(xBin),numel(mult));for i=1:numel(mult),q=apply_spline(Tb,Cb,Sb,xBin,mult(i));fits(:,i)=q(:,1);end
end

function [Tb,Cb,Sb]=raw_bins(mapped,sensorIds,xBin,minCount,minCoverage)
nB=numel(xBin);nS=numel(sensorIds);Tb=nan(nB,nS);Cb=zeros(nB,nS);Sb=nan(nB,nS);dx=median(diff(xBin));
for is=1:nS
    idx=mapped.S_v==sensorIds(is);x=mapped.x_v(idx);v=mapped.V_a(idx);
    ib=round((x-xBin(1))/dx)+1;ok=isfinite(x)&isfinite(v)&ib>=1&ib<=nB;ib=ib(ok);v=v(ok);
    for j=1:nB,q=v(ib==j);Cb(j,is)=numel(q);if numel(q)>=minCount,Tb(j,is)=mean(q);Sb(j,is)=std(q);end,end
    valid=isfinite(Tb(:,is));if nnz(valid)<minCoverage*nB,error('Sensor %d has insufficient bin coverage.',sensorIds(is));end
    Tb(:,is)=fillmissing(Tb(:,is),'linear','EndValues','nearest');
    Sb(:,is)=fillmissing(Sb(:,is),'nearest');
end
end

function T=apply_sg(T,span)
for is=1:size(T,2),s=min(span,size(T,1));if mod(s,2)==0,s=s-1;end;if s>=5,T(:,is)=smoothdata(T(:,is),'sgolay',s,'omitnan');end,end
end

function [T,D,pp,tols]=apply_spline(Tb,Cb,Sb,x,mult)
T=nan(size(Tb));D=T;pp=cell(1,size(Tb,2));tols=zeros(1,size(Tb,2));
for is=1:size(Tb,2)
    varMean=(Sb(:,is).^2)./max(Cb(:,is),1);noiseEnergy=sum(varMean(isfinite(varMean)));
    tol=max(eps,mult*noiseEnergy);sp=spaps(x.',Tb(:,is).',tol);pp{is}=sp;tols(is)=tol;
    T(:,is)=fnval(sp,x).';D(:,is)=fnval(fnder(sp),x).';
end
end

function low=package_low(mapped,xOut,xBin,Tb,Db,Cb,Sb,sensorIds,method,dx)
T=interp1(xBin,Tb,xOut,'pchip','extrap');D=interp1(xBin,Db,xOut,'pchip','extrap');
C=interp1(xBin,Cb,xOut,'nearest','extrap');S=interp1(xBin,Sb,xOut,'nearest','extrap');
low=struct('xGrid',xOut,'templateBySensor',T,'derivativeBySensor',D,...
    'templateLow',mean(T,2,'omitnan'),'derivativeLow',mean(D,2,'omitnan'),...
    'sensorIds',sensorIds,'countBySensor',C,'stdBySensor',S,'countLow',sum(C,2),...
    'stdLow',mean(S,2,'omitnan'),'mapped',mapped,'method',char(method),'grid_spacing_mm',dx);
end

function D=gradient_matrix(T,dx)
D=zeros(size(T));for i=1:size(T,2),D(:,i)=gradient(T(:,i),dx);end
end
function m=submap(m,idx)
names={'t_v','V_a','x_v','rev_v','S_v','theta_v'};for i=1:numel(names),if isfield(m,names{i}),m.(names{i})=m.(names{i})(idx);end,end
end
function v=get_opt(s,n,d)
if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
