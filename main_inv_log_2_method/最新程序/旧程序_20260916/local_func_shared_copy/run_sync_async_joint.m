function R=run_sync_async_joint(highMap,lib,cfg,staticState,syncEO)
%RUN_SYNC_ASYNC_JOINT  Fixed synchronous EO plus inverse-scanned async f.
fsync=syncEO*(cfg.RPM_high/60); gapState=estimate_projected_gap_only(highMap,lib,cfg,staticState);
inv=inverse_map_local_displacement(highMap,lib,gapState.g_used,0,cfg);ok=inv.valid&inv.weight>0;
t=highMap.t_v(ok);u=inv.u_inv(ok);sw=sqrt(inv.weight(ok));
fGrid=cfg.singleFreqGrid(:); coarse=inf(size(fGrid));
for k=1:numel(fGrid)
 if abs(fGrid(k)-fsync)<20,continue;end
 X=[ones(nnz(ok),1),sin(2*pi*fsync*t),cos(2*pi*fsync*t),sin(2*pi*fGrid(k)*t),cos(2*pi*fGrid(k)*t)];b=(X.*sw)\(u.*sw);r=(u-X*b).*sw;coarse(k)=dot(r,r);
end
[~,ord]=sort(coarse);nCand=3;nFine=get_cfg(cfg,'singleSAFineKeep',2);snr=NaN;
if isfield(highMap,'snr_db_equiv'),snr=highMap.snr_db_equiv;end
if ~isfinite(snr)&&isfield(cfg,'snrDb')&&isscalar(cfg.snrDb),snr=cfg.snrDb;end
if isfinite(snr)&&snr<get_cfg(cfg,'singleSALowSNRThreshold',20)
 nCand=get_cfg(cfg,'singleSALowSNRTopK',5);nFine=get_cfg(cfg,'singleSALowSNRFineKeep',3);
end
if isfinite(snr)&&snr<get_cfg(cfg,'singleSAVeryLowSNRThreshold',10)
 nCand=get_cfg(cfg,'singleSAVeryLowSNRTopK',8);nFine=get_cfg(cfg,'singleSAVeryLowSNRFineKeep',4);
end
cand=zeros(1,min(nCand,nnz(isfinite(coarse))));n=0;
for i=1:numel(ord)
 if ~isfinite(coarse(ord(i))),continue;end
 if n==0||all(abs(fGrid(ord(i))-cand(1:n))>=20),n=n+1;cand(n)=fGrid(ord(i));if n==numel(cand),break;end,end
end
fine=[];for k=1:numel(cand),fine=[fine,(cand(k)-10:0.1:cand(k)+10)];end %#ok<AGROW>
fine=unique(fine(fine>=min(fGrid)&fine<=max(fGrid)&abs(fine-fsync)>=20));
[~,q]=sort(arrayfun(@(f) inverse_pair_score(f,t,u,sw,fsync),fine));
fine=fine(q(1:min(nFine,numel(q))));best=[];bestRmse=inf;
for k=1:numel(fine)
 fAsync=fine(k);X=[ones(nnz(ok),1),sin(2*pi*fsync*t),cos(2*pi*fsync*t),sin(2*pi*fAsync*t),cos(2*pi*fAsync*t)];b=(X.*sw)\(u.*sw);
 seed=struct('dx',b(1),'A',[hypot(b(2),b(3)),hypot(b(4),b(5))],'phi',[atan2(b(3),b(2)),atan2(b(5),b(4))]);
 fit=fit_dual_fixed_frequency_ab_multistart(highMap,lib,cfg,gapState,[fsync,fAsync],seed);
 if fit.rmse<bestRmse,bestRmse=fit.rmse;best=fit;end
end
% Release only the asynchronous frequency after the fixed-grid voltage fit.
% The synchronous EO remains an exact physical constraint.
if ~isempty(best)
 tAll=highMap.t_v(:);xAll=highMap.x_v(:);VAll=highMap.V_a(:);au=get_cfg(cfg,'finalAmplitudeUpperBoundMm',1.5);
 gh=get_cfg(cfg,'finalProjectedJointGapHalfWidth',0.12);dh=get_cfg(cfg,'finalProjectedJointDxHalfWidth',0.20);fh=get_cfg(cfg,'singleSAContinuousHalfWidthHz',0.5);
 z0=[best.theta,best.f(2)];lb=[max(0.05,gapState.g_used-gh),-dh,-au*ones(1,4),max(min(fGrid),best.f(2)-fh)];ub=[gapState.g_used+gh,dh,au*ones(1,4),min(max(fGrid),best.f(2)+fh)];
 lb(1)=max(lb(1),min(lib.gapTrain)-cfg.rawGapSearchMargin);ub(1)=min(ub(1),max(lib.gapTrain)+cfg.rawGapSearchMargin);
 [z,si]=solve_lsq_bounded(@(q)sa_cont_res(q,tAll,VAll,xAll,lib,fsync),z0,lb,ub,400,1e-12,1e-12);
 A1=hypot(z(3),z(4));A2=hypot(z(5),z(6));p=[A1,atan2(z(4),z(3)),fsync,A2,atan2(z(6),z(5)),z(7)];Vf=eval_gap_template(lib,z(1),xAll-z(2)-fit_u(p,tAll));rr=VAll-Vf;rm=sqrt(mean(rr(isfinite(rr)).^2));
 if rm<bestRmse,best.g=z(1);best.dx=z(2);best.f=[fsync,z(7)];best.A=[A1,A2];best.phi=p([2,5]);best.p=p;best.theta=z;best.VFit=Vf;best.rmse=rm;best.solve_info=si;end
end
ident=compute_dx_vibration_identifiability(highMap,lib,best.g,best.dx,best.A,best.phi,best.f);
R=struct('method',"inverse_SA_direct_ab",'mode',"SA",'g_used',best.g,'dx_used',best.dx,'f_id',best.f, ...
 'A_id',best.A,'phi_id',best.phi,'rmse',best.rmse,'fit',best,'sync_eo',syncEO,'inverseMap',inv, ...
 'async_candidates',fine,'async_coarse_candidates',cand, ...
 'identifiability',ident,'projectedGapState',gapState,'staticState',staticState,'vpUsed',false);
end
function r=sa_cont_res(z,t,V,x,lib,fs)
A1=hypot(z(3),z(4));A2=hypot(z(5),z(6));u=A1*sin(2*pi*fs*t+atan2(z(4),z(3)))+A2*sin(2*pi*z(7)*t+atan2(z(6),z(5)));
r=V-eval_gap_template(lib,z(1),x-z(2)-u);r(~isfinite(r))=10*max(std(V),1e-3);
end
function s=inverse_pair_score(f,t,u,sw,fs)
X=[ones(numel(t),1),sin(2*pi*fs*t),cos(2*pi*fs*t),sin(2*pi*f*t),cos(2*pi*f*t)];b=(X.*sw)\(u.*sw);r=(u-X*b).*sw;s=dot(r,r);
end
function v=get_cfg(s,n,d),if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end,end
