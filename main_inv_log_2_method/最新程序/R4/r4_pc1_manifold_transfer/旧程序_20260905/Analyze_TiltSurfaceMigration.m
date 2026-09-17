function Results = Analyze_TiltSurfaceMigration()
%ANALYZE_TILTSURFACEMIGRATION Validate complete response-surface migration.
thisDir = fileparts(mfilename('fullpath'));
dataPath = fullfile(thisDir,'output','waveform_families','TiltGap_WaveformSurfaceData.mat');
load(dataPath,'Data');
x = Data.x_mm(:); g = Data.gap_mm(:); a = Data.tilt_deg(:); V = Data.response_V;
nA = numel(a); kList = (2:nA-1).';
names = ["nearest_surface","local_linear_surface","pod_linear_surface","global_quadratic_surface"];
rows = repmat(metric_row(),numel(kList)*numel(names),1); ir = 0;
for k = kList.'
    train = setdiff(1:nA,k); [~,j] = min(abs(a(train)-a(k)));
    P = {V(:,:,train(j)), local_linear(V,a,k), pod_linear(V,a,k,train,3), global_quadratic(V,a,k,train)};
    for im = 1:numel(names)
        ir = ir+1; rows(ir) = evaluate_surface(P{im},V(:,:,k),x,g,a(k),names(im));
    end
end
Metrics = struct2table(rows);
outDir = fullfile(thisDir,'output','surface_migration');
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(Metrics,fullfile(outDir,'tilt_surface_migration_metrics.csv'));
save(fullfile(outDir,'tilt_surface_migration_results.mat'),'Metrics','x','g','a','kList','Data','-v7.3');
fig = figure('Color','w','Units','centimeters','Position',[2 2 17 8], ...
    'Name','Tilt surface migration audit','NumberTitle','off');
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
nexttile; plot_metric(Metrics,'surface_relative_rmse','Full-surface relative RMSE');
nexttile; plot_metric(Metrics,'max_peak_x_error_mm','Maximum SG peak-position error (mm)');
drawnow;
pngPath = fullfile(outDir,'tilt_surface_migration_metrics.png');
pdfPath = fullfile(outDir,'tilt_surface_migration_metrics.pdf');
exportgraphics(fig,pngPath,'Resolution',600); exportgraphics(fig,pdfPath,'ContentType','vector');
Results = struct('metrics',Metrics,'output_dir',string(outDir),'figure',fig, ...
    'png_path',string(pngPath),'pdf_path',string(pdfPath));
disp(Metrics); fprintf('Saved tilt-surface migration audit to:\n%s\n',outDir);
end

function S = local_linear(V,a,k)
w = (a(k)-a(k-1))/(a(k+1)-a(k-1)); S = (1-w)*V(:,:,k-1)+w*V(:,:,k+1);
end

function S = global_quadratic(V,a,k,train)
S = zeros(size(V,1),size(V,2));
for ix=1:size(V,1)
    for ig=1:size(V,2)
        p=polyfit(a(train),squeeze(V(ix,ig,train)),2);
        S(ix,ig)=polyval(p,a(k));
    end
end
end

function S = pod_linear(V,a,k,train,r)
% Interpolate low-dimensional surface coordinates, retaining only r POD modes.
M = reshape(V(:,:,train),[],numel(train));
mu = mean(M,2);
[U,~,~] = svd(M-mu,'econ');
r = min([r,size(U,2)]);
U = U(:,1:r);
C = U'*(M-mu);
iLo = find(train==k-1,1); iHi = find(train==k+1,1);
w = (a(k)-a(k-1))/(a(k+1)-a(k-1));
cPred = (1-w)*C(:,iLo)+w*C(:,iHi);
S = reshape(mu+U*cPred,size(V,1),size(V,2));
end

function r = evaluate_surface(pred,truth,x,g,tiltDeg,method)
xp=nan(1,numel(g)); xt=xp; yp=xp; yt=xp;
for ig=1:numel(g), [xp(ig),yp(ig)]=sg_peak(x,pred(:,ig)); [xt(ig),yt(ig)]=sg_peak(x,truth(:,ig)); end
d=pred(:)-truth(:); scale=norm(truth(:)-mean(truth(:)));
r=struct('tilt_deg',tiltDeg,'method',method,'surface_relative_rmse',norm(d)/max(scale,eps), ...
    'surface_correlation',corr(pred(:),truth(:)),'max_peak_x_error_mm',max(abs(xp-xt)), ...
    'mean_peak_x_error_mm',mean(abs(xp-xt)),'max_peak_amplitude_error_pct',100*max(abs(yp-yt)./yt), ...
    'mean_peak_amplitude_error_pct',100*mean(abs(yp-yt)./yt));
end

function [xp,yp]=sg_peak(x,y)
ys=smoothdata(y(:),'sgolay',11); [~,im]=max(ys); q=max(1,im-2):min(numel(x),im+2); p=polyfit(x(q),ys(q),2);
xp=-p(2)/(2*p(1)); if ~(p(1)<0 && xp>=x(q(1)) && xp<=x(q(end))), xp=x(im); end; yp=polyval(p,xp);
end

function r=metric_row()
r=struct('tilt_deg',NaN,'method',"",'surface_relative_rmse',NaN,'surface_correlation',NaN, ...
    'max_peak_x_error_mm',NaN,'mean_peak_x_error_mm',NaN,'max_peak_amplitude_error_pct',NaN, ...
    'mean_peak_amplitude_error_pct',NaN);
end

function plot_metric(T,field,yLabel)
methods=unique(T.method,'stable'); hold on;
for i=1:numel(methods), q=T.method==methods(i); plot(T.tilt_deg(q),T.(field)(q),'o-','LineWidth',1.0,'DisplayName',strrep(methods(i),'_',' ')); end
hold off; set(gca,'FontName','Times New Roman','FontSize',8.5,'TickDir','in','Box','on','LineWidth',0.65);
xlabel('Held-out FE tilt, \alpha (deg)','FontName','Times New Roman','FontSize',9); ylabel(yLabel,'FontName','Times New Roman','FontSize',9);
legend('Location','northwest','FontSize',7.5); grid off;
end
