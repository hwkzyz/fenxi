function Result = Run_LowHighSingleVPBridge(runMode)
%RUN_LOWHIGHSINGLEVPBRIDGE Shared-SNR low/high single-frequency validation.
if nargin < 1, runMode = "smoke"; end
runMode = lower(string(runMode));
root = fileparts(mfilename('fullpath')); mainDir = fileparts(root);
addpath(mainDir, '-begin'); addpath(fullfile(mainDir, 'local_func'), '-begin');
ctx = load_inv_log_2_project_context(mainDir); trust = load_fixed_trust_domain(mainDir);
lib = make_fixed_trust_template_library(ctx.gapList, ctx.xCell, ctx.yCell, NaN, ...
    ctx.cfgAna.xGridN, get_inv_log_2_model_def(), trust);
base = make_inv_log_2_demo_case(ctx, .2); rf = base.cfgCase.RPM_high / 60;
if runMode == "smoke"
    eoList = [1 2.5 4.5]; gapPairs = [.2 .2; .8 .5]; snrList = 20; N = 600;
else
    eoList = [1:5 2.5 3.5 4.5 5.7]; gapPairs = [.2 .2; .5 .5; .8 .8; .5 .2; .8 .5]; snrList = [Inf 20]; N = 900;
end
rows = repmat(row0(), 0, 1); tAll = tic;
for ig = 1:size(gapPairs,1)
    gLow = gapPairs(ig,1); gHigh = gapPairs(ig,2);
    for isnr = 1:numel(snrList)
        snrDb = snrList(isnr); c = base.cfgCase; c.RPM_low = min(c.RPM_high,300); c.NumRevs_low = 8;
        lowData = simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),c,lib.domain,snrDb,710000+100*ig+isnr);
        low = aggregate_low_speed_template(lowData,c.alpha_k,c.R_tip,lib.domain,lib.xGrid);
        st = estimate_highspeed_static_gap_raw(low.mapped,lib,c);
        for ie = 1:numel(eoList)
            fTrue = rf * eoList(ie);
            d = simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gHigh,z),c.RPM_high,c.NumRevs_high,c.fs,c.R_tip,c.alpha_k,0,lib.domain,.25,fTrue,pi/4,'noise_ratio');
            sig = rms(d.V_clean-min(d.V_clean)); noiseStd = ternary(isinf(snrDb),0,sig/10^(snrDb/20));
            rng(810000+1000*ig+100*isnr+ie,'twister'); d.V_cap=d.V_clean+noiseStd*randn(size(d.V_clean));
            m=map_highspeed_to_space(d,c.alpha_k,c.R_tip,lib.domain,.02); m=subset_map(m,N);
            c.f1Grid=5:5:1500; c.f2Grid=1500; c.numVarproCandidates=12;
            c.mainVpGapHalfWidth=.45; c.mainVpGapN=61; c.mainVpCoarseCount=31; c.mainVpMaxKeep=9;
            c.mainUseProjectedJointMuList=0;
            tc=tic; fit2=run_main_vp_joint(m,lib,c,st);
            fit1=refit_single_model(m,lib,fit2);
            nObs=numel(m.V_a); bic1=nObs*log(max(fit1.rmse^2,eps))+5*log(nObs); bic2=nObs*log(max(fit2.rmse^2,eps))+8*log(nObs);
            chooseSingle=bic1<=bic2; elapsed=toc(tc);
            fEst=fit1.f; r=row0(); r.class=ternary(abs(eoList(ie)-round(eoList(ie)))<1e-9,"sync","async"); r.eo=eoList(ie); r.snr_db=snrDb; r.g_low_true=gLow; r.g_low_hat=st.gHat; r.g_high_true=gHigh; r.g_high_hat=fit1.g; r.g_error_mm=fit1.g-gHigh; r.frequency_true_hz=fTrue; r.frequency_est_hz=fEst; r.frequency_error_hz=fEst-fTrue; r.amplitude_est_mm=fit1.A; r.rmse_V=fit1.rmse; r.elapsed_s=elapsed; r.model_order=ternary(chooseSingle,1,2); r.bic_single=bic1; r.bic_dual=bic2; r.success=chooseSingle&&abs(r.frequency_error_hz)<=1&&abs(r.g_error_mm)<=.05&&abs(r.amplitude_est_mm-.25)<=.06; rows(end+1,1)=r; %#ok<AGROW>
        end
    end
end
Audit=struct2table(rows); out=fullfile(root,'output',[char(runMode) '_single']); if ~exist(out,'dir'),mkdir(out);end; writetable(Audit,fullfile(out,'low_high_single_vp_bridge.csv')); save(fullfile(out,'low_high_single_vp_bridge.mat'),'Audit','-v7.3'); Result=struct('Audit',Audit,'wallTimeSeconds',toc(tAll),'outputDir',out);
end
function m=subset_map(m,N),if numel(m.t_v)<=N,return;end;n=numel(m.t_v);ii=round(linspace(1,n,N));fn=fieldnames(m);for k=1:numel(fn),v=m.(fn{k});if isnumeric(v)&&isvector(v)&&numel(v)==n,m.(fn{k})=v(ii);end;end;end
function fit=refit_single_model(m,lib,fit2)
t=m.t_v(:);x=m.x_v(:);V=m.V_a(:);f0=fit2.f_id(1);a0=fit2.A_id(1)*cos(fit2.phi_id(1));b0=fit2.A_id(1)*sin(fit2.phi_id(1));z0=[fit2.g_used fit2.dx_used a0 b0 f0];lo=max(5,f0-20);hi=min(1500,f0+20);obj=@(z)mean((V-eval_gap_template(lib,z(1),x-z(2)-z(3)*sin(2*pi*z(5)*t)-z(4)*cos(2*pi*z(5)*t))).^2)+1e3*(max(z(1)-1.5,0)^2+max(.05-z(1),0)^2+max(abs(z(2))-.5,0)^2+max(lo-z(5),0)^2+max(z(5)-hi,0)^2);opts=optimset('Display','off','MaxIter',160,'TolX',1e-7,'TolFun',1e-12);z=fminsearch(obj,z0,opts);fit=struct('g',z(1),'dx',z(2),'f',z(5),'A',hypot(z(3),z(4)),'phi',atan2(z(4),z(3)),'rmse',sqrt(obj(z)));
end
function r=row0(),r=struct('class',"",'eo',NaN,'snr_db',NaN,'g_low_true',NaN,'g_low_hat',NaN,'g_high_true',NaN,'g_high_hat',NaN,'g_error_mm',NaN,'frequency_true_hz',NaN,'frequency_est_hz',NaN,'frequency_error_hz',NaN,'amplitude_est_mm',NaN,'rmse_V',NaN,'elapsed_s',NaN,'model_order',NaN,'bic_single',NaN,'bic_dual',NaN,'success',false);end
function y=ternary(q,a,b),if q,y=a;else,y=b;end;end
