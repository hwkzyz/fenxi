function R=run_dual_sync_joint(highMap,lib,cfg,staticState,eoPair)
%RUN_DUAL_SYNC_JOINT  Two synchronous components with fixed integer EOs.
f=sort(eoPair(:).')*(cfg.RPM_high/60);gapState=estimate_projected_gap_only(highMap,lib,cfg,staticState);
inv=inverse_map_local_displacement(highMap,lib,gapState.g_used,0,cfg);ok=inv.valid&inv.weight>0;
t=highMap.t_v(ok);sw=sqrt(inv.weight(ok));X=[ones(nnz(ok),1),sin(2*pi*f(1)*t),cos(2*pi*f(1)*t),sin(2*pi*f(2)*t),cos(2*pi*f(2)*t)];
b=(X.*sw)\(inv.u_inv(ok).*sw);seed=struct('dx',b(1),'A',[hypot(b(2),b(3)),hypot(b(4),b(5))], ...
 'phi',[atan2(b(3),b(2)),atan2(b(5),b(4))]);fit=fit_dual_fixed_frequency_ab_multistart(highMap,lib,cfg,gapState,f,seed);
ident=compute_dx_vibration_identifiability(highMap,lib,fit.g,fit.dx,fit.A,fit.phi,fit.f);
R=struct('method',"inverse_SS_direct_ab",'mode',"SS",'g_used',fit.g,'dx_used',fit.dx,'f_id',fit.f, ...
 'A_id',fit.A,'phi_id',fit.phi,'rmse',fit.rmse,'fit',fit,'eo_id',sort(eoPair(:).'), ...
 'inverseMap',inv,'inverseSeed',seed,'identifiability',ident,'projectedGapState',gapState,'staticState',staticState,'vpUsed',false);
end
