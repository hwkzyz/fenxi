function [Detail,BinBias]=Run_X01SpatialDiscretizationAudit()
%RUN_X01SPATIALDISCRETIZATIONAUDIT Test hard bins against continuous estimators.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','x01_spatial_discretization_audit');if ~exist(outDir,'dir'),mkdir(outDir);end
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
cfg=make_cfg(make_inv_log_2_demo_case(ctx,.2).cfgCase);c=make_case();cfg.f_true=c.eo*cfg.RPM_high/60;cfg.A_true=c.A;cfg.phi_true=c.phi;
reference=identity_model(lib,c.gLow,cfg);geom=simulate_highspeed_from_low_increment(reference,.1,cfg,Inf,'fixed_std',1);daq=make_daq(geom,lib,c,cfg);map=map_highspeed_to_space(daq,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
lowRaw=simulate_low_speed_template(@(z)eval_gap_template(lib,c.gLow,z),cfg,lib.domain,Inf,921001,'fixed_std');mapped=map_highspeed_to_space(lowRaw,cfg.alpha_k,cfg.R_tip,lib.domain,0);
mappedWide=map_with_window(lowRaw,cfg,lib.domain,.52);
% The present simulator passes T_opr_truth directly to map_highspeed_to_space.
% Thus its OPR coordinate is already exact; this audit isolates aggregation,
% while a future timing-jitter model can populate the mapping-error cell.
mapErr=coordinate_replay_error(lowRaw,mapped,cfg,lib.domain);
models={};labels=string.empty;meta=table();
widths=[.005 .01 .02 .04 .08];
for w=widths
    models{end+1}=hard_bin_model(mapped,lib,cfg,w); %#ok<AGROW>
    labels(end+1)="hard-bin "+string(w)+" mm"; %#ok<AGROW>
end
models{end+1}=hard_bin_model(mappedWide,lib,cfg,.02);labels(end+1)="hard-bin 0.02 mm, low map window 0.52";
models{end+1}=turnwise_pchip_model(lowRaw,lib,cfg);labels(end+1)="OPR mapped, no bin (turnwise PCHIP)";
models{end+1}=local_linear_model(mapped,lib,cfg,.02);labels(end+1)="OPR mapped, local-linear h=0.02 mm";
models{end+1}=local_linear_model(mapped,lib,cfg,.04);labels(end+1)="OPR mapped, local-linear h=0.04 mm";
models{end+1}=local_linear_model(mapped,lib,cfg,.06);labels(end+1)="OPR mapped, local-linear h=0.06 mm";
models{end+1}=local_linear_model(mapped,lib,cfg,.08);labels(end+1)="OPR mapped, local-linear h=0.08 mm";
rows=repmat(row0(),numel(models),1);exact=repmat(eval_gap_template(lib,c.gLow,lib.xGrid),1,numel(cfg.alpha_k));
for i=1:numel(models)
    m=models{i};tf=fixed_pair(map,m,cfg,c.eo,c.gHigh);wf=fixed_pair(map,m,cfg,[10 18],c.gHigh);
    rows(i).estimator=labels(i);rows(i).template_rmse_V=sqrt(mean((m.templateLow-exact).^2,'all'));rows(i).true_rmse_V=tf.rmse;rows(i).wrong_rmse_V=wf.rmse;rows(i).delta_rmse_V=wf.rmse-tf.rmse;rows(i).winner=ternary(rows(i).delta_rmse_V>=0,"(10,12)","(10,18)");rows(i).g0_est_mm=m.pathCal.g0;
end
Detail=struct2table(rows);BinBias=bin_bias_metrics(mapped,lib,c.gLow,.02);
% Support-aware control: the benchmark-wide physical envelope is declared
% before EO selection.  It is intentionally separate from loose solver caps.
physical=struct('dx_bounds_mm',[-.05 .05],'amplitude_max_mm',[.25 .25],...
    'gap_bounds_mm',[.40 .70],...
    'source',"Predeclared audit envelope; independent of EO winner");
Clegacy=derive_support_aware_domain(map,mapped,lib,physical);
Cwide=derive_support_aware_domain(map,mappedWide,lib,physical);
wideModel=hard_bin_model(mappedWide,lib,cfg,.02);certMap=restrict_map(map,Cwide.high_certified_mm);
tf=fixed_pair(certMap,wideModel,cfg,c.eo,c.gHigh);wf=fixed_pair(certMap,wideModel,cfg,[10 18],c.gHigh);
supportRow=struct2table(struct('estimator',"support-aware hard-bin 0.02 mm",...
    'template_rmse_V',sqrt(mean((wideModel.templateLow-exact).^2,'all')),...
    'true_rmse_V',tf.rmse,'wrong_rmse_V',wf.rmse,'delta_rmse_V',wf.rmse-tf.rmse,...
    'winner',ternary(wf.rmse>=tf.rmse,"(10,12)","(10,18)"),'g0_est_mm',wideModel.pathCal.g0));
Detail=[Detail;supportRow];
Contract=table(string(["legacy_low_window";"support_aware_low_window"]),...
    [Clegacy.low_observed_mm(1);Cwide.low_observed_mm(1)],...
    [Clegacy.low_observed_mm(2);Cwide.low_observed_mm(2)],...
    [Clegacy.high_certified_mm(1);Cwide.high_certified_mm(1)],...
    [Clegacy.high_certified_mm(2);Cwide.high_certified_mm(2)],...
    [Clegacy.full_high_domain_supported;Cwide.full_high_domain_supported],...
    'VariableNames',{'calibration','low_support_min_mm','low_support_max_mm','high_certified_min_mm','high_certified_max_mm','observed_high_fully_supported'});
writetable(Detail,fullfile(outDir,'summary.csv'));writetable(BinBias,fullfile(outDir,'bin_bias.csv'));save(fullfile(outDir,'audit.mat'),'Detail','BinBias','mapErr','-v7.3');
 writetable(Contract,fullfile(outDir,'support_contract.csv'));
write_report(outDir,mapErr);disp(Detail);disp(mapErr);
end

function cfg=make_cfg(cfg)
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=20;cfg.NumRevs_high=8;cfg.route30ForwardModel="low_increment";cfg.route30GapHalfWidthMm=.20;cfg.route30DualFrequencyRangeHz=[500 1300];cfg.dualSyncGapCount=11;cfg.funnelAllPairsFirstJointStep=true;cfg.funnelSecondJointPairCount=24;cfg.funnelFullRefineCount=1;cfg.dualSyncMaxIter=600;cfg.structuredComponentAmplitudeFloorMm=.05;cfg.returnDualSyncDiagnostics=true;
end
function c=make_case(),c=struct('eo',[10 12],'A',[.25 .10],'phi',[0 pi/2],'gLow',.5,'gHigh',.7);end
function m=identity_model(lib,g,cfg),T=repmat(eval_gap_template(lib,g,lib.xGrid),1,numel(cfg.alpha_k));pc=struct('g0',g,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg.alpha_k));m=build_path_template_model(lib,T,pc);end

function m=hard_bin_model(M,lib,cfg,w)
[xBin,T,C,S,sids]=raw_bins(M,lib.domain,w);Tout=interp1(xBin,T,lib.xGrid,'pchip','extrap');m=make_model(Tout,lib,C,S,xBin,sids,.5);
end
function m=turnwise_pchip_model(lowRaw,lib,cfg)
low=aggregate_low_speed_template(lowRaw,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid);m=make_model(low.templateLow,lib,low.countLow,low.stdLow,low.xGrid,1,.5);
end
function m=local_linear_model(M,lib,~,h)
xout=lib.xGrid(:);sids=unique(M.S_v(:),'stable');T=nan(numel(xout),numel(sids));C=zeros(size(T));S=nan(size(T));
for is=1:numel(sids)
    q=M.S_v==sids(is);x=M.x_v(q);v=M.V_a(q);
    for j=1:numel(xout)
        d=x-xout(j);keep=abs(d)<=h;z=d(keep);y=v(keep);C(j,is)=numel(y);S(j,is)=std(y);
        if numel(y)>=3
            w=(1-(z/h).^3).^3;X=[ones(numel(z),1) z];b=(X.*w)\(y.*w);T(j,is)=b(1);
        end
    end
    T(:,is)=fillmissing(T(:,is),'linear','EndValues','nearest');S(:,is)=fillmissing(S(:,is),'nearest');
end
m=make_model(T,lib,C,S,xout,sids,.5);
end
function m=make_model(T,lib,C,S,x,sids,g0)
if isvector(T),T=T(:);end
C=match_size(C,size(T));S=match_size(S,size(T));low=struct('xGrid',lib.xGrid,'templateLow',mean(T,2),'templateBySensor',T,'countLow',sum(C,2),'stdLow',mean(S,2));opts=struct('xDomain',lib.domain,'fitMask',[true false false false],'fitGainOffset',false,'g0Init',g0,'sigmaFloor',1);[pc,~]=calibrate_low_speed_path(lib,low,opts);pc.sensorIds=sids;m=build_path_template_model(lib,T,pc);
end
function A=match_size(A,sz)
if isvector(A),A=repmat(A(:),1,sz(2));elseif size(A,1)~=sz(1),A=interp1(linspace(0,1,size(A,1)),A,linspace(0,1,sz(1))','nearest','extrap');end
end
function [xBin,T,C,S,sids]=raw_bins(M,domain,w)
xBin=(domain(1):w:domain(2)).';if xBin(end)<domain(2),xBin(end+1)=domain(2);end;sids=unique(M.S_v(:),'stable');T=nan(numel(xBin),numel(sids));C=zeros(size(T));S=nan(size(T));
for is=1:numel(sids),q=M.S_v==sids(is);ib=round((M.x_v(q)-xBin(1))/w)+1;v=M.V_a(q);ok=ib>=1&ib<=numel(xBin);ib=ib(ok);v=v(ok);for j=1:numel(xBin),z=v(ib==j);C(j,is)=numel(z);if numel(z)>=3,T(j,is)=mean(z);S(j,is)=std(z);end,end;T(:,is)=fillmissing(T(:,is),'linear','EndValues','nearest');S(:,is)=fillmissing(S(:,is),'nearest');end
end
function D=make_daq(geom,lib,c,cfg)
t=geom.t(:);om=cfg.RPM_high*2*pi/60;h=.52*diff(lib.domain);u=c.A(1)*sin(2*pi*cfg.f_true(1)*t+c.phi(1))+c.A(2)*sin(2*pi*cfg.f_true(2)*t+c.phi(2));V=min(eval_gap_template(lib,c.gLow,lib.xGrid))*ones(size(t));Vs=V;for r=1:cfg.NumRevs_high,for s=1:numel(cfg.alpha_k),te=geom.T_opr_truth(r)+cfg.alpha_k(s)/om;idx=find(abs(t-te)<=h/(om*cfg.R_tip));x=om*cfg.R_tip*(t(idx)-te);q=x-u(idx);ok=q>=lib.domain(1)&q<=lib.domain(2);if any(ok),V(idx(ok))=eval_gap_template(lib,c.gHigh,q(ok),s);end;ok=x>=lib.domain(1)&x<=lib.domain(2);if any(ok),Vs(idx(ok))=eval_gap_template(lib,c.gHigh,x(ok),s);end,end,end;D=geom;D.V_clean=V;D.V_static=Vs;D.V_cap=V;
end
function fit=fixed_pair(map,m,cfg,pair,g),q=cfg;q.dualSyncCandidateEOPairs=pair;state=struct('gHat',g,'g_low_hat',g,'dx0',0);fit=run_dual_sync_voltage_vp_funnel(map,m,q,state).fit;end
function B=bin_bias_metrics(M,lib,g,w)
[x,T,C,~,~]=raw_bins(M,lib.domain,w);truth=eval_gap_template(lib,g,x);bias=mean(T,2)-truth;mu=zeros(numel(x),1);v=mu;for j=1:numel(x),q=round((M.x_v-x(1))/w)+1==j;d=M.x_v(q)-x(j);mu(j)=mean(d,'omitnan');v(j)=mean(d.^2,'omitnan');end;fx=eval_gap_derivative(lib,g,x);B=table(x,mean(C,2),mu,v,bias,fx.*mu,'VariableNames',{'x_mm','mean_count','mean_offset_mm','mean_square_offset_mm2','template_bias_V','first_order_prediction_V'});
end
function E=coordinate_replay_error(D,M,cfg,domain)
t=D.t(:);T=D.T_opr_truth(:);err=[];for r=1:numel(T)-1,om=2*pi/(T(r+1)-T(r));for s=1:numel(cfg.alpha_k),te=T(r)+cfg.alpha_k(s)/om;idx=find(abs(t-te)<=.48*diff(domain)/(om*cfg.R_tip));x=om*cfg.R_tip*(t(idx)-te);err=[err;x];end,end;E=struct('mapping_uses_truth_opr',true,'coordinate_replay_max_abs_error_mm',0,'note',"The current simulator maps with T_opr_truth; a separate OPR perturbation model is required to study mapping error.");
end
function M=restrict_map(M,domain)
keep=M.x_v>=domain(1)&M.x_v<=domain(2);names={'t_v','V_a','x_v','rev_v','S_v','theta_v'};for i=1:numel(names),M.(names{i})=M.(names{i})(keep);end
end
function M=map_with_window(D,cfg,domain,ratio)
T=D.T_opr_truth(:);t=D.t(:);V=D.V_cap(:);h=ratio*diff(domain);tv=[];vv=[];xv=[];rv=[];sv=[];th=[];
for r=1:numel(T)-1
    om=2*pi/(T(r+1)-T(r));
    for s=1:numel(cfg.alpha_k)
        te=T(r)+cfg.alpha_k(s)/om;idx=abs(t-te)<=h/(om*cfg.R_tip);x=om*cfg.R_tip*(t(idx)-te);keep=x>=domain(1)&x<=domain(2);idx=find(idx);idx=idx(keep);x=x(keep);
        tv=[tv;t(idx)];vv=[vv;V(idx)];xv=[xv;x];rv=[rv;r*ones(numel(x),1)];sv=[sv;s*ones(numel(x),1)];th=[th;2*pi*(r-1)+om*(t(idx)-T(r))];
    end
end
M=struct('t_v',tv,'V_a',vv,'x_v',xv,'rev_v',rv,'S_v',sv,'theta_v',th);
end
function r=row0(),r=struct('estimator',"",'template_rmse_V',NaN,'true_rmse_V',NaN,'wrong_rmse_V',NaN,'delta_rmse_V',NaN,'winner',"",'g0_est_mm',NaN);end
function y=ternary(tf,a,b),if tf,y=a;else,y=b;end,end
function write_report(outDir,E),fid=fopen(fullfile(outDir,'README.md'),'w');fprintf(fid,'# X01 Spatial Discretization Audit\n\nThe synthetic OPR mapping currently consumes `T_opr_truth`, so the OPR-coordinate error is exactly zero by construction. This audit therefore identifies hard-bin versus continuous-template effects; a timing-jitter generator is needed before a mapping-error claim can be tested.\n');fclose(fid);end
