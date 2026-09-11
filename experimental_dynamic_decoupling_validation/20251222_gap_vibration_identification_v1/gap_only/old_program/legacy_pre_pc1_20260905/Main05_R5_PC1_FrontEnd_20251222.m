%% Frozen R5 experimental front-end
% Main05 is the only calibration/transfer entry in the current route.
% It builds the PC1 surface family, performs B2 registration, localizes
% the target surface from the low-speed anchor, and exports Main10 sidecars.
clc
rootDir = fileparts(mfilename('fullpath'));
r5Dir = fullfile(rootDir, 'analysis', 'r5_anchor_guided_experimental_20260903');
addpath(r5Dir);
responseFile = fullfile(rootDir,'inputs','calibration','Step05_Response_Surface_20251222.mat');
% Use the package-local short alias to keep MATLAB v7.3/HDF5 paths below the
% Windows MAX_PATH boundary.  It is a byte-for-byte copy of the frozen
% AdaptiveSG template with the long historical filename.
templateFile = fullfile(rootDir,'inputs','prepared','foundation','step04_low_speed_template', ...
    'Template_AdaptiveSG.mat');
outRoot = fullfile(r5Dir,'results','r4_frontend_only');
if ~exist(outRoot,'dir'), mkdir(outRoot); end
manifest = struct();
manifest.methodFreezeId = 'R4R5_PC1_ANCHOR_SIDECAR_20260905';
manifest.generatedAt = datestr(now,30);
manifest.packageRoot = rootDir;
manifest.frontEnd = 'target-excluded PC1 surface family -> B2 affine registration -> shared-z profile localization';
manifest.targetBlades = [1 5];
manifest.entries = struct([]);
for targetBlade = [1 5]
    entryIndex = find([1 5]==targetBlade,1);
    targetRoot = fullfile(outRoot,sprintf('B%d',targetBlade));
    if ~exist(targetRoot,'dir'), mkdir(targetRoot); end
    familyFile = fullfile(targetRoot,'r4_experimental_surface_family_lobo.mat');
    registrationFile = fullfile(targetRoot,'r4_b2_anchor_registration.mat');
    R5_Build_R4_ExperimentalSurfaceFamily(responseFile,familyFile,targetBlade);
    FFcheck = load(familyFile,'F');
    assert(double(FFcheck.F.excludedTargetBlade)==double(targetBlade) && ...
        ~ismember(double(targetBlade),double(FFcheck.F.trainingBladeIds)), ...
        'R5:LOBOContract','Family for B%d does not satisfy target-excluded training contract.',targetBlade);
    R5_Calibrate_B2_AnchorRegistration(familyFile,templateFile,registrationFile,1.0);
    localizationDir = fullfile(targetRoot,'surface_localization');
    R5_Localize_R4_Surface_FromLowSpeed(familyFile,templateFile, ...
        localizationDir,1.0,registrationFile);
    localizationFile = fullfile(localizationDir,'loc.mat');
    baseLibrary = fullfile(rootDir,'inputs','calibration',sprintf( ...
        'Step06I_OffsetTiltShared_GapLibrary_20251222_B%d_S123.mat',targetBlade));
    outputFile = fullfile(outRoot,sprintf('R5_Main10_CompatibleLibrary_B%d_S123.mat',targetBlade));
    R5_Export_Main10_CompatibleLibrary(familyFile,localizationFile, ...
        baseLibrary,outputFile,targetBlade);
    FF = load(familyFile,'F');
    manifest.entries(entryIndex).targetBlade = targetBlade;
    manifest.entries(entryIndex).excludedTargetBlade = targetBlade;
    manifest.entries(entryIndex).trainingBladeIds = FF.F.trainingBladeIds;
    manifest.entries(entryIndex).familyFile = strrep(familyFile,[rootDir filesep],'');
    manifest.entries(entryIndex).registrationFile = strrep(registrationFile,[rootDir filesep],'');
    manifest.entries(entryIndex).localizationFile = strrep(localizationFile,[rootDir filesep],'');
    manifest.entries(entryIndex).sidecarFile = strrep(outputFile,[rootDir filesep],'');
end
% B2 is a target-independent experimental-domain adapter. The two LOBO
% branches must therefore produce the same registration within numerical
% precision; otherwise registration has become target-dependent.
G1 = load(fullfile(rootDir,manifest.entries(1).registrationFile),'G');
G5 = load(fullfile(rootDir,manifest.entries(2).registrationFile),'G');
assert(isequaln(G1.G.registration,G5.G.registration), ...
    'R5:RegistrationNotShared','B1/B5 LOBO runs produced different B2 registrations.');
% Write the formal package contract next to the generated sidecars.  Main06
% deliberately fails closed when this manifest is absent or from another
% frozen method version.
baseContract = fullfile(rootDir,'inputs','calibration', ...
    'Step06I_OffsetTiltShared_GapLibrary_20251222_B1_S123.mat');
Bc = load(baseContract,'CorrectedGapLibrary');
manifest.sidecarFiles = { ...
    strrep(fullfile(outRoot,'R5_Main10_CompatibleLibrary_B1_S123.mat'),[rootDir filesep],''), ...
    strrep(fullfile(outRoot,'R5_Main10_CompatibleLibrary_B5_S123.mat'),[rootDir filesep],'')};
manifest.sensor = Bc.CorrectedGapLibrary.sensor;
save(fullfile(outRoot,'R5_Formal_Manifest.mat'),'manifest','-v7');
fprintf('R5 front-end complete. Main10 sidecars are under: %s\n', outRoot);
