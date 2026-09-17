function Result = Run_LowHighVPBridge(runMode)
%RUN_LOWHIGHVPBRIDGE Low-speed gap calibration feeding high-speed VP fitting.
if nargin < 1, runMode = "smoke"; end
runMode = lower(string(runMode));
root = fileparts(mfilename('fullpath')); mainDir = fileparts(root);
addpath(mainDir, '-begin'); addpath(fullfile(mainDir, 'local_func'), '-begin');
ctx = load_inv_log_2_project_context(mainDir); trust = load_fixed_trust_domain(mainDir);
lib = make_fixed_trust_template_library(ctx.gapList, ctx.xCell, ctx.yCell, NaN, ...
    ctx.cfgAna.xGridN, get_inv_log_2_model_def(), trust);
base = make_inv_log_2_demo_case(ctx, .2);
cases = {[500 1300], [700 1200], [900 1400]};
freqStep = 25;
refineCount = 5;
candidateDiversityHz = 8;
seedList = 1;
adaptiveDense = false;
if runMode == "smoke"
    gapPairs = [.2 .2; .8 .8; .8 .5]; N = 600; snrList = 20;
elseif runMode == "non_grid"
    cases = {[531 1300], [733 1217], [917 1383]};
    gapPairs = [.2 .2; .8 .5]; N = 900; snrList = 20; freqStep = 5; refineCount = 5; candidateDiversityHz = 8;
elseif runMode == "non_grid_repeat"
    cases = {[531 1300], [733 1217], [917 1383]};
    gapPairs = [.2 .2; .8 .5]; N = 900; snrList = 20; freqStep = 5; seedList = 1:3;
elseif runMode == "low_snr"
    cases = {[531 1300], [733 1217], [917 1383]};
    gapPairs = [.2 .2; .8 .5]; N = 900; snrList = [15 10]; freqStep = 5;
elseif runMode == "very_low_snr"
    cases = {[531 1300], [733 1217], [917 1383]};
    gapPairs = [.2 .2; .8 .5]; N = 1500; snrList = 5; freqStep = 5; adaptiveDense = true;
elseif runMode == "very_low_snr_repeat"
    cases = {[531 1300], [733 1217], [917 1383]};
    gapPairs = [.2 .2; .8 .5]; N = 1500; snrList = 5; freqStep = 5; seedList = 1:3; adaptiveDense = true;
elseif runMode == "very_low_snr_seed2"
    cases = {[531 1300], [733 1217], [917 1383]};
    gapPairs = [.2 .2; .8 .5]; N = 1500; snrList = 5; freqStep = 5; seedList = 2; adaptiveDense = true;
elseif runMode == "very_low_snr_seed2_dense"
    cases = {[531 1300], [733 1217], [917 1383]};
    gapPairs = [.2 .2; .8 .5]; N = 2500; snrList = 5; freqStep = 5; seedList = 2;
elseif runMode == "random_frequency"
    cases = {[347 863], [428 1129], [612 1471], [689 1043], [842 1267], [991 1436]};
    gapPairs = [.2 .2; .8 .5]; N = 900; snrList = 20; freqStep = 5;
elseif runMode == "random_parameters"
    cases = {[383 941], [527 1189], [671 1327], [809 1453]};
    ampCases = [.10 .08; .35 .08; .12 .30; .30 .25];
    phaseCases = [.10 2.20; pi/2 -.20; -2.0 1.0; 2.8 -2.4];
    gapPairs = [.2 .2; .8 .5]; N = 900; snrList = 20; freqStep = 5;
elseif runMode == "random_parameters_dense"
    cases = {[383 941], [527 1189], [671 1327], [809 1453]};
    ampCases = [.10 .08; .35 .08; .12 .30; .30 .25];
    phaseCases = [.10 2.20; pi/2 -.20; -2.0 1.0; 2.8 -2.4];
    gapPairs = [.2 .2; .8 .5]; N = 1500; snrList = 20; freqStep = 5;
elseif runMode == "full_robust"
    gapPairs = [.2 .2; .5 .5; .8 .8; .5 .2; .8 .5]; N = 900; snrList = [Inf 20]; freqStep = 5;
else
    gapPairs = [.2 .2; .5 .5; .8 .8; .5 .2; .8 .5]; N = 900; snrList = [Inf 20];
end
if ~exist('ampCases','var'),ampCases=repmat([.25 .15],numel(cases),1);end
if ~exist('phaseCases','var'),phaseCases=repmat([pi/4 -pi/3],numel(cases),1);end
rows = repmat(row0(), 0, 1); tAll = tic;
for ig = 1:size(gapPairs,1)
    gLow = gapPairs(ig,1); gHigh = gapPairs(ig,2);
    lowCfg = base.cfgCase; lowCfg.RPM_low = min(lowCfg.RPM_high, 300); lowCfg.NumRevs_low = 8;
    for isnr = 1:numel(snrList)
        snrDb = snrList(isnr);
        for irep=1:numel(seedList)
            noiseSeed=seedList(irep);
            lowData = simulate_low_speed_template(@(z) eval_gap_template(lib, gLow, z), lowCfg, lib.domain, snrDb, 880000+1000*noiseSeed+100*ig+isnr);
            low = aggregate_low_speed_template(lowData, lowCfg.alpha_k, lowCfg.R_tip, lib.domain, lib.xGrid);
            lowState = estimate_highspeed_static_gap_raw(low.mapped, lib, lowCfg);
            for ic = 1:numel(cases)
        fTrue = cases{ic}; ATrue=ampCases(ic,:); phiTrue=phaseCases(ic,:); c = base.cfgCase;
        d = simulate_rotating_waveform_from_template(@(z) eval_gap_template(lib, gHigh, z), c.RPM_high, c.NumRevs_high, c.fs, c.R_tip, c.alpha_k, 0, lib.domain, ATrue, fTrue, phiTrue, 'noise_ratio');
        sig = rms(d.V_clean - min(d.V_clean));
        noiseStd = ternary(isinf(snrDb), 0, sig / 10^(snrDb/20));
        rng(990000+1000*noiseSeed+100*ig+isnr, 'twister'); d.V_cap = d.V_clean + noiseStd * randn(size(d.V_clean));
        mFull = map_highspeed_to_space(d, c.alpha_k, c.R_tip, lib.domain, .02); m = subset_map(mFull, N);
        % Use an unordered full-band coarse pair grid. The legacy split grid
        % excludes 700 and 900 Hz from the first component by construction.
        c.f1Grid = 300:freqStep:1500; c.f2Grid = 300:freqStep:1500; c.vpUniqueUnorderedPairs = true; c.vpUseStagedGrid = false;
        c.numVarproCandidates = 12;
        c.mainVpGapHalfWidth = .45; c.mainVpGapN = 61; c.mainVpCoarseCount = 31; c.mainVpMaxKeep = 9;
        c.mainVpGlobalFreqRefine = true; c.mainVpGlobalFreqRefineGapCount = 3;
        c.mainVpRefineCount = refineCount;
        c.vpCandidateDiversityHz = candidateDiversityHz;
        c.mainUseProjectedJointMuList = 0;
        st = lowState; st.gHat = lowState.gHat; st.dx0 = lowState.dx0;
        ticCase = tic; fit = run_main_vp_joint(m, lib, c, st); fallbackUsed=false;
        weakRatio=min(fit.A_id)/max(max(fit.A_id),eps);
        needDense=(adaptiveDense && fit.rhoJ<0.003) || (snrDb<=20 && weakRatio<0.30);
        if needDense && numel(mFull.t_v)>N
            denseN=ternary(snrDb<=5,2500,1500); mDense=subset_map(mFull,denseN); fit=run_main_vp_joint(mDense,lib,c,st); fallbackUsed=true;
        end
        elapsed = toc(ticCase);
        r = row0(); r.case_name = sprintf('%d+%d', fTrue); r.noise_seed=noiseSeed; r.g_low_true = gLow; r.g_low_hat = lowState.gHat; r.g_high_true = gHigh; r.g_high_hat = fit.g_used; r.g_error_mm = fit.g_used-gHigh; r.snr_db=snrDb; r.f1_true=fTrue(1); r.f2_true=fTrue(2); r.f1_est=fit.f_id(1); r.f2_est=fit.f_id(2); r.max_frequency_error=max(abs(fit.f_id(1:2)-fTrue)); r.a1_true=ATrue(1);r.a2_true=ATrue(2);r.a1_est=fit.A_id(1);r.a2_est=fit.A_id(2);r.max_amplitude_error=max(abs(fit.A_id(1:2)-ATrue));r.max_phase_error=max(abs(angle(exp(1i*(fit.phi_id(1:2)-phiTrue))))); r.rmse_V=fit.rmse; r.deltaJ=fit.deltaJ; r.rhoJ=fit.rhoJ; r.fallback_used=fallbackUsed; r.elapsed_s=elapsed; r.success=r.max_frequency_error<=10 && abs(r.g_error_mm)<=.05 && r.max_amplitude_error<=.06 && r.max_phase_error<=.5; rows(end+1,1)=r; %#ok<AGROW>
        end
        end
    end
end
Audit = struct2table(rows); out=fullfile(root,'output',char(runMode)); if ~exist(out,'dir'),mkdir(out);end; writetable(Audit,fullfile(out,'low_high_vp_bridge.csv')); save(fullfile(out,'low_high_vp_bridge.mat'),'Audit','-v7.3'); Result=struct('Audit',Audit,'wallTimeSeconds',toc(tAll),'outputDir',out);
end
function m=subset_map(m,N), if numel(m.t_v)<=N,return;end; n=numel(m.t_v); ii=round(linspace(1,n,N)); fn=fieldnames(m); for k=1:numel(fn),v=m.(fn{k});if isnumeric(v)&&isvector(v)&&numel(v)==n,m.(fn{k})=v(ii);end;end; end
function r=row0(),r=struct('case_name',"",'noise_seed',NaN,'g_low_true',NaN,'g_low_hat',NaN,'g_high_true',NaN,'g_high_hat',NaN,'g_error_mm',NaN,'snr_db',NaN,'f1_true',NaN,'f2_true',NaN,'f1_est',NaN,'f2_est',NaN,'max_frequency_error',NaN,'a1_true',NaN,'a2_true',NaN,'a1_est',NaN,'a2_est',NaN,'max_amplitude_error',NaN,'max_phase_error',NaN,'rmse_V',NaN,'deltaJ',NaN,'rhoJ',NaN,'fallback_used',false,'elapsed_s',NaN,'success',false);end
function y=ternary(q,a,b),if q,y=a;else,y=b;end;end
