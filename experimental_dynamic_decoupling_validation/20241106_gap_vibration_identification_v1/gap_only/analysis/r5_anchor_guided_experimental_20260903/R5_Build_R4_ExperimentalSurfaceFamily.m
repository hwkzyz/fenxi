function F = R5_Build_R4_ExperimentalSurfaceFamily(responseFile, outputFile, excludeBlade)
% Build an R4-compatible experimental waveform-surface family.
% Production files are read only. The hidden surface coordinate is blade
% response geometry, not a physical tilt angle.
if nargin < 1 || isempty(responseFile)
    responseFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'inputs','calibration','Step05_Response_Surface_20251222.mat');
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(fileparts(mfilename('fullpath')),'results', ...
        'r4_experimental_surface_family.mat');
end
if nargin < 3, excludeBlade = []; end
sourceResponseFile = responseFile;
responseFile = R5_StageMatForMatlab(responseFile,'r5_response_surface');
S = load(responseFile,'responseSurface');
R = S.responseSurface;
x = double(R.xGrid(:));
g = double(R.trueGapMm(:));
Y = double(R.waveforms);                 % x, blade, gap
% Main05 stores the static calibration waveform in mV, while Main07's
% formal template is baseline-removed mV. Keep one voltage contract before
% any anchor distance is computed.
if isfield(R,'voltageUnit') && strcmpi(char(string(R.voltageUnit)),'V')
    Y = 1000*Y;
end
blade = double(R.bladeIds(:).');
if ~isempty(excludeBlade)
    keepBlade = blade ~= double(excludeBlade);
    assert(any(keepBlade),'R5:EmptyLOBOTraining','No training blades remain after excluding B%d.',excludeBlade);
    blade = blade(keepBlade);
    Y = Y(:,keepBlade,:);
else
    keepBlade = true(size(blade));
end
support = true(size(x));
if isfield(R,'effectiveWindow'), support = logical(R.effectiveWindow(:)); end
gref = double(R.g0Mm);
X = [ones(numel(g),1), 1./g, log(g/gref)];
nX = numel(x); nB = numel(blade);
Q = nan(3*nX,nB); Yfit = nan(size(Y)); fitRmse = nan(nB,1);
for ib=1:nB
    B = (X \ squeeze(Y(:,ib,:)).').';
    Q(:,ib) = B(:);
    Yfit(:,ib,:) = reshape((X*B.').',nX,1,numel(g));
    e = squeeze(Yfit(:,ib,:))-squeeze(Y(:,ib,:));
    fitRmse(ib) = sqrt(mean(e(support,:).^2,'all','omitnan'));
end
mu = mean(Q,2,'omitnan'); Qc = Q-mu;
[U,Sv,~] = svd(Qc,'econ');
r = min(2,size(U,2)); Ur = U(:,1:r); Z = Ur.'*Qc;
[z1,ord] = sort(Z(1,:));
F = struct();
F.schema = 'R5_R4_EXPERIMENTAL_SURFACE_FAMILY_V1';
F.method = 'within_blade_gap_law_plus_PC1_ordered_coefficient_manifold';
F.sourceFile = sourceResponseFile; F.stagedSourceFile = responseFile;
F.x_mm=x; F.gap_mm=g; F.gref_mm=gref;
F.trainingBladeIds = blade; F.excludedTargetBlade = excludeBlade;
F.blade_id=blade; F.waveforms_mv=Y; F.waveforms_fit_mv=Yfit;
F.B0 = Q(1:nX,:); F.B1=Q(nX+1:2*nX,:); F.B2=Q(2*nX+1:3*nX,:);
F.q=Q; F.mean_q=mu; F.basis=Ur; F.latent=Z; F.latent_sorted=z1;
F.blade_order_sorted=blade(ord); F.support_mask=support; F.fit_rmse_mv=fitRmse;
F.units.voltage='mV'; F.units.x='mm'; F.units.gap='mm';
F.identifiability_note = ['latent coordinate is a waveform-manifold coordinate; ' ...
    'it is not a measurable tilt angle or an absolute local clearance.'];
if ~exist(fileparts(outputFile),'dir'), mkdir(fileparts(outputFile)); end
save(outputFile,'F','-v7.3');
fprintf('R5 R4 experimental surface family saved: %s\n',outputFile);
end
