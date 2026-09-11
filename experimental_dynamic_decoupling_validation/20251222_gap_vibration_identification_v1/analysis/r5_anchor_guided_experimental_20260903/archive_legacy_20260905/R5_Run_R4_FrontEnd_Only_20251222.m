function out = R5_Run_R4_FrontEnd_Only_20251222()
% Front-end-only R5 run. Main08/Main09/Main10 remain untouched.
root=fileparts(mfilename('fullpath')); pkg=fileparts(fileparts(root));
responseFile=fullfile(pkg,'inputs','calibration','Step05_Response_Surface_20251222.mat');
templateFile=fullfile(pkg,'inputs','prepared','foundation','step04_low_speed_template', ...
 'Template_OPRCenterStd_AdaptiveSG_LowSpeed_AllBlades_S123_20251222.mat');
outRoot=fullfile(root,'results','r4_frontend_only');
familyFile=fullfile(outRoot,'r4_experimental_surface_family.mat');
registrationFile=fullfile(outRoot,'r4_b2_anchor_registration.mat');
R5_Build_R4_ExperimentalSurfaceFamily(responseFile,familyFile);
G=R5_Calibrate_B2_AnchorRegistration(familyFile,templateFile,registrationFile,1.0);
T=R5_Localize_R4_Surface_FromLowSpeed(familyFile,templateFile,fullfile(outRoot,'surface_localization'),1.0,registrationFile);
localizationFile=fullfile(outRoot,'surface_localization','R5_R4_surface_localization.mat');
main10Libraries=struct();
for targetBlade=[1 5]
    baseLibrary=fullfile(pkg,'inputs','calibration',sprintf( ...
        'Step06I_OffsetTiltShared_GapLibrary_20251222_B%d_S123.mat',targetBlade));
    main10Library=fullfile(outRoot,sprintf('R5_Main10_CompatibleLibrary_B%d_S123.mat',targetBlade));
    R5_Export_Main10_CompatibleLibrary(familyFile,localizationFile,baseLibrary,main10Library,targetBlade);
    main10Libraries.(sprintf('B%d',targetBlade))=main10Library;
end
out=struct('familyFile',familyFile,'registration',G,'localization',T, ...
    'main10Library',main10Libraries.B1,'main10Libraries',main10Libraries,'outputRoot',outRoot);
save(fullfile(outRoot,'R5_R4_frontend_result.mat'),'out','-v7.3');
end
