function [Detail,Summary,Selection]=Run_AdaptiveLowSpeedTemplateBenchmark(snrList,seeds,numLowRevs,numHighRevs)
%RUN_ADAPTIVELOWSPEEDTEMPLATEBENCHMARK Compare fixed/adaptive templates.
% Only Gaussian noise is simulated. All template hyperparameters use the
% low-speed record alone; noiseless truth is touched only by pack_metrics.
if nargin<1||isempty(snrList),snrList=[25 15 5];end
if nargin<2||isempty(seeds),seeds=1:3;end
if nargin<3||isempty(numLowRevs),numLowRevs=20;end
if nargin<4||isempty(numHighRevs),numHighRevs=8;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','adaptive_low_speed_template_benchmark');if ~exist(outDir,'dir'),mkdir(outDir);end
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=numLowRevs;cfg0.NumRevs_high=numHighRevs;cfg0.route30ForwardModel="low_increment";cfg0.A_true=[.25 .15];
fRot=cfg0.RPM_high/60;eoTrue=[10 26];cfg0.f_true=eoTrue*fRot;cfg0.phi_true=[pi/4 -pi/3];cfg0.route30DualFrequencyRangeHz=eoTrue*fRot;cfg0.route30GapHalfWidthMm=.20;cfg0.dualSyncGapCount=11;cfg0.dualSyncKeepPerGap=80;cfg0.dualSyncIteratedReplayCount=200;cfg0.dualSyncRefineCount=20;cfg0.structuredComponentAmplitudeFloorMm=.05;cfg0.returnDualSyncDiagnostics=false;
g0=.50;delta=.10;nSensor=numel(cfg0.alpha_k);truthT=repmat(eval_gap_template(lib,g0,lib.xGrid),1,nSensor);dx=median(diff(lib.xGrid));truthD=gradmat(truthT,dx);
truthCal=struct('g0',g0,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:nSensor);truthModel=build_path_template_model(lib,truthT,truthCal);
methods={'A_fixed_SG9','B_groupCV_SG','C_groupCV_spline','D_split_forward_screen'};n=numel(snrList)*numel(seeds)*numel(methods);Detail=repmat(empty_row(),n,1);Selection=table();ir=0;
for isnr=1:numel(snrList)
  snrDb=snrList(isnr);
  for js=1:numel(seeds)
    seed=seeds(js);lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,g0,z),cfg0,lib.domain,snrDb,941000+seed);
    high=simulate_highspeed_from_low_increment(truthModel,delta,cfg0,snrDb,'snr_db',942000+seed);highMap=map_highspeed_to_space(high,cfg0.alpha_k,cfg0.R_tip,lib.domain,.02);
    % Use the acquisition support contract used by the current simulator.
    % Keeping this explicit makes the benchmark independent of legacy defaults.
    opt=struct('gridSpacingMm',.02,'minBinCount',5,'minCoverage',.90,'numFolds',5,'mapWindowRatio',.52);
    lowA=build_low_speed_templates_binned(lowData,cfg0.alpha_k,cfg0.R_tip,lib.domain,lib.xGrid,'mean',setfield(opt,'smoothSpan',9)); %#ok<SFLD>
    lowA.derivativeBySensor=gradmat(lowA.templateBySensor,dx);lowA.selected_window_mm=9*.02;
    lowB=build_low_speed_templates_adaptive(lowData,cfg0.alpha_k,cfg0.R_tip,lib.domain,lib.xGrid,'cv_sg',opt);
    lowC=build_low_speed_templates_adaptive(lowData,cfg0.alpha_k,cfg0.R_tip,lib.domain,lib.xGrid,'cv_spline',opt);
    % D: retain B's forward template, but allow a stronger derivative only
    % when its held-out voltage error lies within a two-SE low-speed guard.
    cv=lowB.cv_candidates;[m,imin]=min(cv.cv_rmse_V);ok=cv.cv_rmse_V<=m+2*cv.cv_se_V(imin);iScreen=find(ok,1,'last');
    lowScreen=build_low_speed_templates_binned(lowData,cfg0.alpha_k,cfg0.R_tip,lib.domain,lib.xGrid,'mean',setfield(opt,'smoothSpan',cv.span_bins(iScreen))); %#ok<SFLD>
    lowScreen.derivativeBySensor=gradmat(lowScreen.templateBySensor,dx);screenWindow=cv.window_mm(iScreen);
    lows={lowA,lowB,lowC,lowB};screen={lowA,lowB,lowC,lowScreen};
    sel=table(repmat(snrDb,3,1),repmat(seed,3,1),["B_groupCV_SG";"C_groupCV_spline";"D_EO_screen"],...
      [lowB.selected_window_mm;NaN;screenWindow],[lowB.selected_span;NaN;cv.span_bins(iScreen)],...
      [NaN;lowC.selected_tolerance_multiplier;NaN],'VariableNames',{'snr_db','seed','method','selected_window_mm','selected_span_bins','spline_tolerance_multiplier'});Selection=[Selection;sel]; %#ok<AGROW>
    for im=1:numel(methods)
      tStart=tic;low=lows{im};opts=struct('xDomain',lib.domain,'fitMask',[true false false false],'fitGainOffset',false,'g0Init',g0,'sigmaFloor',1);[pc,cal]=calibrate_low_speed_path(lib,low,opts);pc.sensorIds=low.sensorIds;
      pc.lowDerivativeBySensor=low.derivativeBySensor;model=build_path_template_model(lib,low.templateBySensor,pc);
      pcs=pc;pcs.lowDerivativeBySensor=screen{im}.derivativeBySensor;screenModel=build_path_template_model(lib,low.templateBySensor,pcs);
      state=estimate_highspeed_static_gap_raw(low.mapped,lib,cfg0);state.gHat=pc.g0;state.g_low_hat=pc.g0;state.dx0=0;
      cfgFixed=cfg0;cfgFixed.dualSyncCandidateEOPairs=eoTrue;tf=tic;if im==4,fitFixed=run_dual_sync_voltage_vp_split_derivative(highMap,model,screenModel,cfgFixed,state);else,fitFixed=run_dual_sync_voltage_vp(highMap,model,cfgFixed,state);end;fixedTime=toc(tf);
      cfgFull=cfg0;cfgFull.dualSyncCandidateEOPairs=[];te=tic;if im==4,fitFull=run_dual_sync_voltage_vp_split_derivative(highMap,model,screenModel,cfgFull,state);else,fitFull=run_dual_sync_voltage_vp(highMap,model,cfgFull,state);end;fullTime=toc(te);
      if im==4,sw=screenWindow;else,sw=NaN;end
      ir=ir+1;Detail(ir)=pack_metrics(snrDb,seed,methods{im},low,truthT,truthD,lib.xGrid,pc,cal,fitFixed,fitFull,g0,delta,cfg0.A_true,eoTrue,fixedTime,fullTime,toc(tStart),sw);
      writetable(struct2table(Detail(1:ir)),fullfile(outDir,'checkpoint_detail.csv'));save(fullfile(outDir,'checkpoint.mat'),'Detail','Selection','ir');
    end
  end
end
T=struct2table(Detail);Summary=summarize(T,snrList,methods);writetable(T,fullfile(outDir,'detail.csv'));writetable(Summary,fullfile(outDir,'summary.csv'));writetable(Selection,fullfile(outDir,'selection.csv'));save(fullfile(outDir,'benchmark.mat'),'Detail','Summary','Selection','truthT','truthD','snrList','seeds');
plot_results(Summary,snrList,methods,outDir);write_readme(Summary,Selection,outDir);disp(Summary);
end

function r=empty_row()
r=struct('snr_db',NaN,'seed',NaN,'method','','selected_window_mm',NaN,'screen_window_mm',NaN,'template_rmse_V',NaN,'weighted_derivative_rmse_V_per_mm',NaN,'derivative_correlation',NaN,'derivative_gain',NaN,'peak_bias_pct',NaN,'width_bias_pct',NaN,'area_bias_pct',NaN,'slope_bias_pct',NaN,'g0_error_mm',NaN,'fixed_delta_error_mm',NaN,'fixed_A1_error_mm',NaN,'fixed_A2_error_mm',NaN,'full_eo1',NaN,'full_eo2',NaN,'full_eo_correct',false,'full_delta_error_mm',NaN,'full_A1_error_mm',NaN,'full_A2_error_mm',NaN,'fixed_runtime_s',NaN,'full_runtime_s',NaN,'total_runtime_s',NaN,'calibration_rmse_V',NaN);
end
function r=pack_metrics(snr,seed,name,low,Tt,Dt,x,pc,cal,ff,fe,g0,delta,A,eo,tf,te,tt,screenW)
D=low.derivativeBySensor;w=Dt.^2;r=empty_row();r.snr_db=snr;r.seed=seed;r.method=name;if isfield(low,'selected_window_mm'),r.selected_window_mm=low.selected_window_mm;end;r.screen_window_mm=screenW;
r.template_rmse_V=sqrt(mean((low.templateBySensor-Tt).^2,'all'));r.weighted_derivative_rmse_V_per_mm=sqrt(sum(w.*(D-Dt).^2,'all')/sum(w,'all'));r.derivative_correlation=sum(D.*Dt,'all')/sqrt(sum(D.^2,'all')*sum(Dt.^2,'all'));r.derivative_gain=sum(D.*Dt,'all')/sum(Dt.^2,'all');
[p0,w0,a0,s0]=shape(Tt,Dt,x);[p1,w1,a1,s1]=shape(low.templateBySensor,D,x);r.peak_bias_pct=100*mean((p1-p0)./p0);r.width_bias_pct=100*mean((w1-w0)./w0);r.area_bias_pct=100*mean((a1-a0)./a0);r.slope_bias_pct=100*mean((s1-s0)./s0);
r.g0_error_mm=abs(pc.g0-g0);r.fixed_delta_error_mm=abs((ff.g_used-pc.g0)-delta);r.fixed_A1_error_mm=abs(ff.A_id(1)-A(1));r.fixed_A2_error_mm=abs(ff.A_id(2)-A(2));q=sort(fe.eo_id(:)).';r.full_eo1=q(1);r.full_eo2=q(2);r.full_eo_correct=isequal(q,eo);r.full_delta_error_mm=abs((fe.g_used-pc.g0)-delta);r.full_A1_error_mm=abs(fe.A_id(1)-A(1));r.full_A2_error_mm=abs(fe.A_id(2)-A(2));r.fixed_runtime_s=tf;r.full_runtime_s=te;r.total_runtime_s=tt;r.calibration_rmse_V=cal.rmse_V;
end
function [pk,wd,ar,sl]=shape(T,D,x)
n=size(T,1);edge=max(3,round(.1*n));pk=zeros(1,size(T,2));wd=pk;ar=pk;sl=pk;dx=median(diff(x));for i=1:size(T,2),b=median([T(1:edge,i);T(end-edge+1:end,i)]);q=abs(T(:,i)-b);pk(i)=max(q);wd(i)=nnz(q>=.5*pk(i))*dx;ar(i)=trapz(x,q);sl(i)=max(abs(D(:,i)));end
end
function S=summarize(T,snrList,methods)
vars={'template_rmse_V','weighted_derivative_rmse_V_per_mm','derivative_correlation','derivative_gain','peak_bias_pct','width_bias_pct','area_bias_pct','slope_bias_pct','g0_error_mm','fixed_delta_error_mm','fixed_A1_error_mm','fixed_A2_error_mm','full_delta_error_mm','full_A1_error_mm','full_A2_error_mm','full_runtime_s'};
S=table();for s=snrList,for m=1:numel(methods),idx=T.snr_db==s&strcmp(T.method,methods{m});row=table(s,string(methods{m}),nnz(idx),mean(T.full_eo_correct(idx)),'VariableNames',{'snr_db','method','n','full_eo_correct_rate'});for v=vars,q=T.(v{1})(idx);row.([v{1},'_median'])=median(q);row.([v{1},'_iqr'])=iqr(q);end;S=[S;row];end,end %#ok<AGROW>
end
function plot_results(S,snrList,methods,outDir)
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 14]);tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');cc=lines(numel(methods));fields={'template_rmse_V_median','weighted_derivative_rmse_V_per_mm_median','derivative_correlation_median','full_eo_correct_rate'};scale=[1000 1 1 100];yl={'Template RMSE (mV)','Weighted derivative RMSE (V/mm)','Derivative correlation','Correct EO pair (%)'};
for p=1:4,nexttile;hold on;box on;for m=1:numel(methods),q=zeros(size(snrList));for i=1:numel(snrList),idx=S.snr_db==snrList(i)&strcmp(S.method,methods{m});q(i)=S.(fields{p})(idx);end;plot(1:numel(snrList),scale(p)*q,'-o','LineWidth',1.1,'MarkerSize',4,'Color',cc(m,:),'DisplayName',strrep(methods{m},'_',' '));end;set(gca,'XTick',1:numel(snrList),'XTickLabel',string(snrList));xlabel('Low/high SNR (dB)');ylabel(yl{p});if p==1,legend('Location','best','FontSize',7);end;end
title(tl,'Low-speed-only regularization benchmark','FontWeight','normal');set(findall(fig,'-property','FontName'),'FontName','Times New Roman');set(findall(fig,'Type','axes'),'TickDir','in','FontSize',8,'Box','on');exportgraphics(fig,fullfile(outDir,'benchmark_overview.png'),'Resolution',300);exportgraphics(fig,fullfile(outDir,'benchmark_overview.pdf'),'ContentType','vector');
fig2=figure('Color','w','Units','centimeters','Position',[2 2 17 7]);tiledlayout(1,4,'TileSpacing','compact','Padding','compact');fs={'peak_bias_pct_median','width_bias_pct_median','area_bias_pct_median','slope_bias_pct_median'};labs={'Peak bias (%)','Half-height width bias (%)','Area bias (%)','Maximum-slope bias (%)'};for p=1:4,nexttile;hold on;box on;for m=1:numel(methods),q=zeros(size(snrList));for i=1:numel(snrList),idx=S.snr_db==snrList(i)&strcmp(S.method,methods{m});q(i)=S.(fs{p})(idx);end;plot(1:numel(snrList),q,'-o','LineWidth',1.1,'MarkerSize',4,'Color',cc(m,:));end;yline(0,'k:');set(gca,'XTick',1:numel(snrList),'XTickLabel',string(snrList));xlabel('SNR (dB)');ylabel(labs{p});end;set(findall(fig2,'-property','FontName'),'FontName','Times New Roman');set(findall(fig2,'Type','axes'),'TickDir','in','FontSize',8,'Box','on');exportgraphics(fig2,fullfile(outDir,'shape_bias.png'),'Resolution',300);exportgraphics(fig2,fullfile(outDir,'shape_bias.pdf'),'ContentType','vector');
end
function write_readme(~,Sel,outDir)
fid=fopen(fullfile(outDir,'RESULTS.md'),'w');fprintf(fid,'# Adaptive low-speed template benchmark\n\n');fprintf(fid,'All hyperparameters were selected from grouped low-speed revolutions. No noiseless template or high-speed result entered selection.\n\n');fprintf(fid,'- A: 0.02 mm bin mean + fixed SG9 (0.18 mm).\n- B: grouped-revolution CV over physical SG windows; minimum held-out prediction error.\n- C: cubic smoothing spline; grouped-CV noise-tolerance multiplier and analytic derivative.\n- D: B forward template plus a separately regularized EO-screen derivative, admitted by a two-SE voltage-prediction guard.\n\n');fprintf(fid,'Selected SG windows (mm), median [range]:\n');for s=unique(Sel.snr_db).',q=Sel.selected_window_mm(Sel.snr_db==s&Sel.method=="B_groupCV_SG");fprintf(fid,'- %g dB: %.2f [%.2f, %.2f]\n',s,median(q),min(q),max(q));end;fprintf(fid,'\nSee summary.csv and detail.csv for all metrics.\n');fclose(fid);
end
function D=gradmat(T,dx)
D=zeros(size(T));for i=1:size(T,2),D(:,i)=gradient(T(:,i),dx);end
end
