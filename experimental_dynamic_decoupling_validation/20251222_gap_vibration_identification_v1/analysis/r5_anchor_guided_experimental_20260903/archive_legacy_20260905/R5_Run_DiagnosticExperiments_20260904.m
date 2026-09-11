function R5_Run_DiagnosticExperiments_20260904()
% Independent R5 diagnostics. Does not alter Main05/Main07/Main10 or inputs.
thisDir = fileparts(mfilename('fullpath'));
pkg = fileparts(fileparts(thisDir));
outDir = fullfile(thisDir, 'results', 'r5_diagnostics_20260904');
figDir = fullfile(outDir, 'figures');
if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(figDir), mkdir(figDir); end

familyFile = fullfile(thisDir, 'results', 'r4_frontend_only', 'r4_experimental_surface_family.mat');
registrationFile = fullfile(thisDir, 'results', 'r4_frontend_only', 'r4_b2_anchor_registration.mat');
templateFile = fullfile(pkg, 'inputs', 'prepared', 'foundation', 'step04_low_speed_template', ...
    'Template_OPRCenterStd_AdaptiveSG_LowSpeed_AllBlades_S123_20251222.mat');
baseLibraryFile = fullfile(pkg, 'inputs', 'calibration', ...
    'Step06I_OffsetTiltShared_GapLibrary_20251222_B1_S123.mat');
r5LibraryFile = fullfile(thisDir, 'results', 'r4_frontend_only', 'R5_Main10_CompatibleLibrary_B1_S123.mat');

load(familyFile, 'F');
load(registrationFile, 'G');
load(templateFile, 'Template');

% A/B/C isolates only the waveform-family representation. All cases use an
% anchor at a known static gap, so this is not a claim about gap identifiability.
staticCases = [run_static_validation_local(F, 'same_blade_leave_gap'); ...
    run_static_validation_local(F, 'leave_one_blade')];
writetable(staticCases, fullfile(outDir, 'R5_ModelOrder_StaticCases.csv'));
staticSummary = summarize_static_local(staticCases);
writetable(staticSummary, fullfile(outDir, 'R5_ModelOrder_StaticSummary.csv'));
plot_static_representatives_local(staticCases, figDir);

% Joint (gap, latent) surfaces use the experimental low-speed templates.
% They expose ridges; they do not turn an unobserved local gap into ground truth.
[jointRows, jointDetail, commonRows, bootstrapRows] = joint_localization_diagnostics_local(F, G, Template);
writetable(jointRows, fullfile(outDir, 'R5_JointGapLatent_Localization.csv'));
writetable(commonRows, fullfile(outDir, 'R5_IndependentVsCommonLatent.csv'));
writetable(bootstrapRows, fullfile(outDir, 'R5_B2_ResidualBootstrap.csv'));
save(fullfile(outDir, 'R5_JointGapLatentSurfaces.mat'), 'jointRows', 'jointDetail', 'commonRows', ...
    'bootstrapRows', 'familyFile', 'registrationFile', 'templateFile', '-v7.3');
plot_joint_representatives_local(jointDetail, figDir);

% Contract checks are deliberately structural plus an exact coefficient-forward
% identity test. The latter checks the R5 export, not the identification result.
contract = library_contract_local(baseLibraryFile, r5LibraryFile, F);
writetable(contract, fullfile(outDir, 'R5_LibraryContract.csv'));
registrationComparison = registration_comparison_local(baseLibraryFile, G);
writetable(registrationComparison, fullfile(outDir, 'R5_B2Registration_vs_Main10.csv'));

write_report_local(outDir, staticSummary, jointRows, commonRows, bootstrapRows, contract, registrationComparison);
fprintf('R5 diagnostic outputs written to %s\n', outDir);
end

function rows = run_static_validation_local(F, mode)
x = F.x_mm(:); g = F.gap_mm(:).'; phi = [ones(size(g)); 1 ./ g; log(g ./ F.gref_mm)];
Y = double(F.waveforms_mv); blade = 1:size(Y,2);  % x, blade, gap in Step05
rows = table();
for ib = 1:numel(blade)
    for ig = 1:numel(g)
        switch mode
            case 'same_blade_leave_gap'
                testBlade = ib; testGap = ig;
                trainG = g; trainG(ig) = [];
                keepG = setdiff(1:numel(g), ig);
                trainY = Y(:, :, keepG);
            case 'leave_one_blade'
                testBlade = ib; testGap = ig;
                keepB = setdiff(1:numel(blade), ib);
                trainY = Y(:, keepB, :); trainG = g;
            otherwise
                error('Unknown mode.');
        end
        [Q, mu, U, score] = fit_family_local(permute(trainY,[1 3 2]), trainG, F.gref_mm);
        yTrue = Y(:, testBlade, testGap);
        for method = ["A_PC1_interpolated", "B_strict_POD1", "C_strict_POD2"]
            [qHat, anchorRmse] = estimate_q_known_gap_local(Q, mu, U, score, trainG, F.gref_mm, yTrue, g(testGap), char(method));
            B = reshape(qHat, numel(x), 3);
            yHat = B * phi(:, testGap);
            rmse = sqrt(mean((yHat-yTrue).^2));
            nrmse = 100 * rmse / max(range(yTrue), eps);
            row = table(string(mode), blade(testBlade), g(testGap), method, rmse, nrmse, anchorRmse, ...
                {x}, {yTrue}, {yHat}, 'VariableNames', {'validation','blade_id','held_gap_mm','method', ...
                'rmse_mv','nrmse_percent','anchor_rmse_mv','x_mm','y_true_mv','y_hat_mv'});
            rows = [rows; row]; %#ok<AGROW>
        end
    end
end
end

function [Q, mu, U, score] = fit_family_local(Y, g, gref)
nx = size(Y, 1); nb = size(Y, 3); X = [ones(numel(g),1), 1./g(:), log(g(:)./gref)];
Q = nan(nx*3, nb);
for ib = 1:nb
    % Keep the same coefficient memory layout as the production sidecar:
    % q = [B0(x); B1(x); B2(x)] = B(:), with B sized [nx,3].
    B = (X \ squeeze(Y(:, :, ib)).').';
    Q(:, ib) = B(:);
end
mu = mean(Q, 2);
[U,~,~] = svd(Q-mu, 'econ');
score = U.' * (Q-mu);
end

function [qHat, bestRmse] = estimate_q_known_gap_local(Q, mu, U, score, g, gref, y, ga, method)
nx = numel(y); h = [1; 1/ga; log(ga/gref)];
switch method
    case 'A_PC1_interpolated'
        z = score(1,:); [zs, ord] = sort(z);
        zGrid = linspace(zs(1), zs(end), 601);
        bestRmse = inf; qHat = Q(:,ord(1));
        for k = 1:numel(zGrid)
            q = interp1(zs, Q(:,ord).', zGrid(k), 'linear', 'extrap').';
            yh = reshape(q,nx,3) * h; r = sqrt(mean((yh-y).^2));
            if r < bestRmse, bestRmse = r; qHat = q; end
        end
    otherwise
        if strcmp(method, 'B_strict_POD1'), r = 1; else, r = min(2, size(U,2)); end
        A = zeros(nx, r);
        for j = 1:r, A(:,j) = reshape(U(:,j), nx, 3) * h; end
        b = y - reshape(mu,nx,3)*h;
        s = A \ b;
        lo = min(score(1:r,:), [], 2); hi = max(score(1:r,:), [], 2);
        s = min(max(s, lo), hi);
        qHat = mu + U(:,1:r)*s;
        bestRmse = sqrt(mean((reshape(qHat,nx,3)*h-y).^2));
end
end

function S = summarize_static_local(T)
methods = unique(T.method); validations = unique(T.validation); S = table();
for iv=1:numel(validations), for im=1:numel(methods)
    q=T.nrmse_percent(T.validation==validations(iv) & T.method==methods(im));
    r=T.rmse_mv(T.validation==validations(iv) & T.method==methods(im));
    S=[S; table(validations(iv), methods(im), numel(q), median(r), prctile(r,90), prctile(r,95), max(r), ...
        median(q), prctile(q,90), prctile(q,95), max(q), 'VariableNames', {'validation','method','case_count', ...
        'rmse_p50_mv','rmse_p90_mv','rmse_p95_mv','rmse_max_mv','nrmse_p50_pct','nrmse_p90_pct','nrmse_p95_pct','nrmse_max_pct'})]; %#ok<AGROW>
end, end
end

function plot_static_representatives_local(T, figDir)
for iv = unique(T.validation).'
    tt = T(T.validation==iv,:);
    spread = groupsummary(tt, {'blade_id','held_gap_mm'}, 'range', 'rmse_mv');
    [~, k] = max(spread.range_rmse_mv);
    b=spread.blade_id(k); g=spread.held_gap_mm(k); sel=tt(tt.blade_id==b & tt.held_gap_mm==g,:);
    f=figure('Visible','off','Color','w'); hold on; grid on; box on;
    plot(sel.x_mm{1},sel.y_true_mv{1},'k','LineWidth',1.6,'DisplayName','measured');
    for j=1:height(sel), plot(sel.x_mm{j},sel.y_hat_mv{j},'LineWidth',1.1,'DisplayName',char(sel.method(j))); end
    xlabel('x / mm'); ylabel('baseline-removed voltage / mV');
    title(sprintf('%s: blade %g, held gap %.3f mm', char(iv), b, g),'Interpreter','none'); legend('Location','best');
    exportgraphics(f,fullfile(figDir,sprintf('R5_ModelOrder_%s_B%d_g%.3f.png',char(iv),b,g)),'Resolution',180); close(f);
end
end

function [rows, detail, commonRows, bootstrapRows] = joint_localization_diagnostics_local(F, G, Template)
zGrid = linspace(min(F.latent_sorted), max(F.latent_sorted), 301);
gGrid = linspace(min(F.gap_mm), max(F.gap_mm), 281);
phi = [ones(1,numel(gGrid)); 1./gGrid; log(gGrid./F.gref_mm)];
rows=table(); detail=struct([]); commonRows=table(); bootstrapRows=table(); nD=0;
[~,familyOrder]=ismember(F.blade_order_sorted,F.blade_id);
qSorted=F.q(:,familyOrder);
bladeIds=unique([Template.SensorBlade.blade_id]);
for ib=1:numel(bladeIds)
    bid=bladeIds(ib); perSensor=struct([]);
    items=Template.SensorBlade([Template.SensorBlade.blade_id]==bid);
    for ia=1:numel(items)
        a=items(ia); sid=double(a.sensor_id); reg=G.registration([G.registration.sensor_id]==sid);
        % Template x_grid is already in mm and v_grid is in V.  Apply the
        % calibrated coordinate map before interpolation onto the family grid.
        xRaw = double(a.x_grid(:));
        yRaw = (double(a.v_grid(:))-double(a.baseline))*1000;
        mm=isfinite(xRaw)&isfinite(yRaw);
        if isfield(a,'valid_grid_mask'), mm=mm&logical(a.valid_grid_mask(:)); end
        if isfield(a,'domain_effective_mask'), mm=mm&logical(a.domain_effective_mask(:)); end
        xRegistered=reg.x_scale*(xRaw-reg.tau_mm);
        xAnchor=interp1(xRegistered(mm),yRaw(mm),F.x_mm,'pchip',NaN);
        use=F.support_mask & isfinite(xAnchor);
        y=xAnchor(use); x=F.x_mm(use);
        J=nan(numel(zGrid),numel(gGrid));
        for iz=1:numel(zGrid)
            q=interp1(F.latent_sorted, qSorted.', zGrid(iz), 'linear').';
            B=reshape(q,numel(F.x_mm),3); Bx=B(use,:);
            yh=reg.voltage_gain*(Bx*phi)+reg.voltage_offset_mv;
            if size(yh,1)~=numel(y), error('R5:DiagnosticSize','Anchor/model size mismatch CH%d.',sid); end
            J(iz,:)=sqrt(mean((yh-y).^2,1));
        end
        [jmin,lin]=min(J(:)); [iz,ig]=ind2sub(size(J),lin);
        profG=min(J,[],1); profZ=min(J,[],2);
        nearG=gGrid(profG<=jmin+1); nearZ=zGrid(profZ<=jmin+1);
        nD=nD+1; detail(nD).blade_id=bid; detail(nD).sensor_id=sid; detail(nD).z_grid=zGrid; detail(nD).g_grid=gGrid; detail(nD).J_mv=single(J);
        detail(nD).x_mm=x; detail(nD).y_mv=y; detail(nD).support_index=find(use); detail(nD).best_z=zGrid(iz); detail(nD).best_g=gGrid(ig);
        row=table(bid,sid,zGrid(iz),gGrid(ig),jmin,min(nearG),max(nearG),min(nearZ),max(nearZ),numel(x), ...
          'VariableNames',{'blade_id','sensor_id','best_z','best_gap_mm','best_rmse_mv','gap_1mv_min_mm','gap_1mv_max_mm','z_1mv_min','z_1mv_max','point_count'}); rows=[rows;row]; %#ok<AGROW>
        perSensor(ia).J=J; perSensor(ia).n=numel(x); perSensor(ia).sid=sid;
        if bid==2
            bootstrapRows=[bootstrapRows; bootstrap_b2_local(F,reg,x,y,zGrid,gGrid,J,bid,sid)]; %#ok<AGROW>
        end
    end
    if ~isempty(perSensor)
        [indR, commonR, zc, gs]=compare_common_z_local(perSensor,zGrid,gGrid);
        commonRows=[commonRows; table(bid,indR,commonR,commonR-indR,100*(commonR/indR-1),zc,{gs}, ...
          'VariableNames',{'blade_id','independent_rmse_mv','common_z_rmse_mv','common_minus_independent_mv','relative_penalty_percent','common_z','common_gap_mm_by_sensor'})]; %#ok<AGROW>
    end
end
end

function B = bootstrap_b2_local(F,reg,x,y,zGrid,gGrid,J,bid,sid)
% Residual bootstrap tests template residual sensitivity. It is not a multi-lap CI.
[~,lin]=min(J(:)); [iz,ig]=ind2sub(size(J),lin);
[~,familyOrder]=ismember(F.blade_order_sorted,F.blade_id); qSorted=F.q(:,familyOrder);
q=interp1(F.latent_sorted,qSorted.',zGrid(iz),'linear').';
supportIndex=arrayfun(@(v)find(abs(F.x_mm-v)==min(abs(F.x_mm-v)),1),x);
Bx=reshape(q,numel(F.x_mm),3); Bx=Bx(supportIndex,:);
yhat=reg.voltage_gain*(Bx*[1;1/gGrid(ig);log(gGrid(ig)/F.gref_mm)])+reg.voltage_offset_mv; res=y-yhat;
rng(20260904+sid); nBoot=100; zg=nan(nBoot,1); gg=nan(nBoot,1);
for k=1:nBoot
    ys=yhat+res(randi(numel(res),numel(res),1)); best=inf;
    for iz0=1:numel(zGrid)
        q0=interp1(F.latent_sorted,qSorted.',zGrid(iz0),'linear').';
        B0=reshape(q0,numel(F.x_mm),3); B0=B0(supportIndex,:);
        yh=reg.voltage_gain*(B0*[ones(1,numel(gGrid));1./gGrid;log(gGrid./F.gref_mm)])+reg.voltage_offset_mv;
        [v,j]=min(sqrt(mean((yh-ys).^2,1))); if v<best, best=v; zg(k)=zGrid(iz0); gg(k)=gGrid(j); end
    end
end
B=table(repmat(bid,nBoot,1),repmat(sid,nBoot,1),zg,gg,'VariableNames',{'blade_id','sensor_id','bootstrap_z','bootstrap_gap_mm'});
end

function [ri, rc, zc, gs] = compare_common_z_local(S,zGrid,gGrid)
den=sum([S.n]); indSum=0; for j=1:numel(S), indSum=indSum+S(j).n*min(S(j).J(:)).^2; end; ri=sqrt(indSum/den);
obj=zeros(numel(zGrid),1); gaps=zeros(numel(S),numel(zGrid));
for iz=1:numel(zGrid), for j=1:numel(S), [v,ig]=min(S(j).J(iz,:)); obj(iz)=obj(iz)+S(j).n*v^2; gaps(j,iz)=gGrid(ig); end, end
[v,k]=min(obj); rc=sqrt(v/den); zc=zGrid(k); gs=gaps(:,k).';
end

function plot_joint_representatives_local(D, figDir)
for k=1:numel(D)
    if D(k).blade_id~=2, continue; end
    f=figure('Visible','off','Color','w');
    imagesc(D(k).g_grid,D(k).z_grid,D(k).J_mv); axis xy; colorbar; hold on;
    plot(D(k).best_g,D(k).best_z,'wp','MarkerFaceColor','w','MarkerSize',10);
    xlabel('gap / mm'); ylabel('latent coordinate'); title(sprintf('B2 CH%d joint RMSE surface',D(k).sensor_id));
    exportgraphics(f,fullfile(figDir,sprintf('R5_JointSurface_B2_CH%d.png',D(k).sensor_id)),'Resolution',180); close(f);
end
end

function T = library_contract_local(baseFile,r5File,F)
A=load(baseFile,'CorrectedGapLibrary'); B=load(r5File,'CorrectedGapLibrary'); a=A.CorrectedGapLibrary; b=B.CorrectedGapLibrary;
checks=strings(0,1); pass=false(0,1); value=strings(0,1);
checks(end+1,1)="responseSurface_xGrid_equal"; pass(end+1,1)=isequaln(a.responseSurface.xGrid,b.responseSurface.xGrid); value(end+1,1)=string(numel(a.responseSurface.xGrid));
checks(end+1,1)="responseSurface_gref_equal"; pass(end+1,1)=isequaln(a.responseSurface.g0Mm,b.responseSurface.g0Mm); value(end+1,1)=string(a.responseSurface.g0Mm);
checks(end+1,1)="responseSurface_coeff_shape"; pass(end+1,1)=isequal(size(b.responseSurface.coeff),[numel(F.x_mm),3]); value(end+1,1)=string(mat2str(size(b.responseSurface.coeff)));
checks(end+1,1)="responseSurface_coeff_finite"; pass(end+1,1)=all(isfinite(b.responseSurface.coeff(:))); value(end+1,1)="finite";
preserved={'waveforms','effectiveWindow','dFdgGrid','dFdxGrid'};
for k=1:numel(preserved)
    name=preserved{k};
    checks(end+1,1)="preserved_"+string(name);
    pass(end+1,1)=isequaln(a.responseSurface.(name),b.responseSurface.(name));
    value(end+1,1)=string(mat2str(size(b.responseSurface.(name))));
end
[~,ord]=ismember(F.blade_order_sorted,F.blade_id);
z=b.r5.localized_latent_coordinate;
qExpected=interp1(F.latent_sorted,F.q(:,ord).',z,'linear').';
coeffExpected=reshape(qExpected,numel(F.x_mm),3);
maxCoeffError=max(abs(coeffExpected-b.responseSurface.coeff),[],'all');
checks(end+1,1)="exported_coeff_matches_localized_family"; pass(end+1,1)=maxCoeffError<1e-9; value(end+1,1)=string(maxCoeffError);
checks(end+1,1)="sensor_latent_consistency_gate"; pass(end+1,1)=b.r5.sensor_latent_range_fraction<=b.r5.sensor_consistency_limit_fraction; value(end+1,1)=string(b.r5.sensor_latent_range_fraction);
T=table(checks,pass,value,'VariableNames',{'check','pass','observed'});
end

function T = registration_comparison_local(baseFile,G)
A=load(baseFile,'CorrectedGapLibrary'); sensor=A.CorrectedGapLibrary.sensor;
rows=table();
for k=1:numel(G.registration)
    r=G.registration(k); j=find([sensor.sensorId]==r.sensor_id,1);
    assert(~isempty(j),'R5:Main10SensorMissing','Main10 sensor CH%d is missing.',r.sensor_id);
    s=sensor(j);
    rows=[rows; table(r.sensor_id,r.tau_mm,s.tauMm,r.tau_mm-s.tauMm, ...
        r.x_scale,s.xScale,r.x_scale-s.xScale,r.voltage_gain,s.voltageGain, ...
        r.voltage_gain-s.voltageGain,r.voltage_offset_mv,s.voltageOffsetMv, ...
        r.voltage_offset_mv-s.voltageOffsetMv,r.rmse_mv,s.lowFitRmseMv, ...
        'VariableNames',{'sensor_id','r5_tau_mm','main10_tau_mm','delta_tau_mm', ...
        'r5_x_scale','main10_x_scale','delta_x_scale','r5_voltage_gain', ...
        'main10_voltage_gain','delta_voltage_gain','r5_voltage_offset_mv', ...
        'main10_voltage_offset_mv','delta_voltage_offset_mv','r5_anchor_rmse_mv', ...
        'main10_lowfit_rmse_mv'})]; %#ok<AGROW>
end
T=rows;
end

function write_report_local(outDir,S,J,C,B,K,R)
fid=fopen(fullfile(outDir,'R5_Diagnostic_Report.md'),'w');
fprintf(fid,'# R5 diagnostic results\n\n');
fprintf(fid,'This is a sidecar analysis. No production Main05/Main07/Main10 file was changed.\n\n');
fprintf(fid,'## Model-order comparison\n\n');
fprintf(fid,'%s\n\n',evalc('disp(S)'));
fprintf(fid,'A is the current PC1 ordered/interpolated family. B and C are strict POD-1/POD-2 reconstructions. These tests use a known static anchor gap; they diagnose representation error only.\n\n');
fprintf(fid,'## Joint gap-latent localization\n\n');
fprintf(fid,'The exported CSV contains the minimum and the operational +1 mV profile spans. A broad span is a ridge diagnostic, not a confidence interval or a measured local clearance.\n\n');
fprintf(fid,'## Three-sensor latent test\n\n');
fprintf(fid,'%s\n\n',evalc('disp(C)'));
fprintf(fid,'Independent versus shared latent coordinates are compared with point-count weighted RMS. The residual bootstrap is spatial-residual only because this template artifact stores aggregate grids rather than all raw laps.\n\n');
fprintf(fid,'## B2 registration versus frozen Main10 correction\n\n');
fprintf(fid,'%s\n\n',evalc('disp(R)'));
fprintf(fid,'These parameter sets are numerically different coordinate/voltage maps. Interface compatibility does not establish physical equivalence.\n\n');
fprintf(fid,'## Contract\n\n');
fprintf(fid,'%s\n',evalc('disp(K)')); fclose(fid);
end
