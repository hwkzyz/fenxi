function M = R5_Build_RidgeSidecarLibraries(diagnosticFile,familyFile,baseDir,outputDir,targets,injectSensorG0)
% Build copy-on-write libraries at representative near-optimal (g,z) ridge points.
% This prepares the subsequent invariant-to-dynamics Main10 propagation test.
if nargin<1||isempty(diagnosticFile), diagnosticFile=fullfile(fileparts(mfilename('fullpath')),'results','r5_diagnostics_20260904','R5_JointGapLatentSurfaces.mat'); end
if nargin<2||isempty(familyFile), familyFile=fullfile(fileparts(mfilename('fullpath')),'results','r4_frontend_only','r4_experimental_surface_family.mat'); end
if nargin<3||isempty(baseDir), baseDir=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'inputs','calibration'); end
if nargin<4||isempty(outputDir), outputDir=fullfile(fileparts(mfilename('fullpath')),'results','ridge_sidecars'); end
if nargin<5||isempty(targets), targets=[1 5]; end
if nargin<6||isempty(injectSensorG0), injectSensorG0=false; end
if ~isfolder(outputDir), mkdir(outputDir); end
D=load(diagnosticFile,'jointDetail'); load(familyFile,'F');
manifest=table(); M=struct();
for target=targets
    dd=D.jointDetail([D.jointDetail.blade_id]==target);
    assert(numel(dd)>=2,'R5:RidgeSensors','Need at least two sensors for B%d.',target);
    z=dd(1).z_grid; g=dd(1).g_grid; obj=zeros(numel(z),1); bestGap=zeros(numel(dd),numel(z));
    totalN=0;
    for k=1:numel(dd)
        [v,ig]=min(double(dd(k).J_mv),[],2); n=numel(dd(k).x_mm); obj=obj+n*v.^2; totalN=totalN+n; bestGap(k,:)=g(ig);
    end
    obj=sqrt(obj/max(totalN,1)); jmin=min(obj); admissible=obj<=jmin+1;
    admissibleIdx=find(admissible);
    izPick=admissibleIdx(round(linspace(1,numel(admissibleIdx),5)));
    [~,bestIdx]=min(obj);
    if ~ismember(bestIdx,izPick)
        [~,nearestInterior]=min(abs(izPick(2:4)-bestIdx));
        izPick(nearestInterior+1)=bestIdx;
        izPick=sort(izPick);
    end
    for p=1:numel(izPick)
        iz0=izPick(p); zz=z(iz0); gg=bestGap(:,iz0).'; Jk=obj(iz0);
        [~,ord]=ismember(F.blade_order_sorted,F.blade_id);
        q=interp1(F.latent_sorted,F.q(:,ord).',zz,'linear').'; coeff=reshape(q,numel(F.x_mm),3);
        base=fullfile(baseDir,sprintf('Step06I_OffsetTiltShared_GapLibrary_20251222_B%d_S123.mat',target)); S=load(base,'CorrectedGapLibrary'); R=S.CorrectedGapLibrary.responseSurface;
        assert(isequal(size(R.coeff),size(coeff)) && isequaln(R.xGrid(:),F.x_mm(:)),'R5:RidgeContract','Base library contract mismatch.');
        R.coeff=coeff; R.method=sprintf('R5 ridge sidecar B%d point %d',target,p); S.CorrectedGapLibrary.responseSurface=R;
        % The unchanged Main10 contract keeps vLow fixed at its original
        % calibration baseline. Therefore surface-only propagation is the
        % valid default; changing sensor.g0Mm is an explicit diagnostic stress
        % test and is expected to require a corresponding baseline rebuild.
        if injectSensorG0
            for sidx=1:numel(S.CorrectedGapLibrary.sensor)
                sid=S.CorrectedGapLibrary.sensor(sidx).sensorId; ir=find([dd.sensor_id]==sid,1);
                if ~isempty(ir), S.CorrectedGapLibrary.sensor(sidx).g0Mm=gg(ir); end
            end
        end
        S.CorrectedGapLibrary.r5_ridge=struct('target_blade',target,'ridge_index',p,'effective_gap_mm_by_sensor',gg,'latent_coordinate',zz,'joint_rmse_mv',Jk,'joint_rmse_min_mv',jmin,'profile_delta_mv',Jk-jmin,'profile_threshold_mv',1,'sensor_g0_injected',injectSensorG0);
        file=fullfile(outputDir,sprintf('R5_RidgeSidecar_B%d_P%d.mat',target,p)); CorrectedGapLibrary=S.CorrectedGapLibrary; save(file,'CorrectedGapLibrary','-v7.3');
        manifest=[manifest; table(target,p,zz,Jk,Jk-jmin,gg(1),gg(2),gg(3),string(fullfile(outputDir,sprintf('R5_RidgeSidecar_B%d_P%d.mat',target,p))), ...
            'VariableNames',{'target_blade','ridge_index','latent_coordinate','joint_rmse_mv','delta_j_mv','gap_s1_mm','gap_s2_mm','gap_s3_mm','library_file'})]; %#ok<AGROW>
    end
end
writetable(manifest,fullfile(outputDir,'R5_RidgeSidecar_Manifest.csv')); M=manifest;
end
