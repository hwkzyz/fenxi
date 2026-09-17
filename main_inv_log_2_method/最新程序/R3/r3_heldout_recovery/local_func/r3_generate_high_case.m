function [high,highMap,truth] = r3_generate_high_case(R3,gTruth,Atrue,phiInjected,sigmaV,seed)
%R3_GENERATE_HIGH_CASE Generate raw-node FE truth and its mapped phase label.

cfg=R3.cfg;Fraw=r3_raw_response(R3,gTruth);
deltaF=0;if isfield(R3.P,'deltaFTrueHz'),deltaF=R3.P.deltaFTrueHz;end
fTrue=R3.P.eoTrue*cfg.RPM_high/60+deltaF;
if Atrue==0,Agen=[];fGen=[];phiGen=[];else,Agen=Atrue;fGen=fTrue;phiGen=phiInjected;end
rng(seed,'twister');
high=simulate_rotating_waveform_from_template(Fraw,cfg.RPM_high,...
    cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,sigmaV,R3.libCal.domain,...
    Agen,fGen,phiGen,'fixed_std');
highMap=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,R3.libCal.domain,.02);

if Atrue>0
    u=Atrue*sin(2*pi*fTrue*highMap.t_v(:)+phiInjected);
    B=[sin(R3.P.eoTrue*highMap.theta_v(:)),cos(R3.P.eoTrue*highMap.theta_v(:))];
    coef=B\u;
    amplitudeMapped=hypot(coef(1),coef(2));
    phaseMapped=atan2(coef(2),coef(1));
    mappingResidualRms=rms(u-B*coef);
else
    amplitudeMapped=0;phaseMapped=NaN;mappingResidualRms=0;
end
truth=struct('g_reference_mm',R3.P.gReferenceMm,'g_truth_mm',gTruth,...
    'delta_gap_truth_mm',gTruth-R3.P.gReferenceMm,'A_true_mm',Atrue,...
    'A_truth_mapped_mm',amplitudeMapped,'eo_true',R3.P.eoTrue,...
    'delta_f_true_hz',deltaF,'f_true_hz',fTrue,'phi_injected_rad',phiInjected,...
    'phi_mapped_rad',phaseMapped,'phase_mapping_residual_rms_mm',mappingResidualRms,...
    'high_seed',seed,'noise_std_V',sigmaV);
end
