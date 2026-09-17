function fit = fit_vp_main_staged(t,V,x,templateLib,gHat,cfg,precomp,gapIdx)
%FIT_VP_MAIN_STAGED Coarse frequency VP followed by local 5 Hz refinement.
cfg0 = cfg; cfg0.vpUseStagedGrid = false;
coarse = fit_vp_main_base(t,V,x,templateLib,gHat,cfg0,precomp,gapIdx);
keep = min(6,numel(coarse.candidates));
f1 = []; f2 = [];
for k=1:keep
    f1 = [f1, coarse.candidates(k).f(1) + (-20:5:20)]; %#ok<AGROW>
    f2 = [f2, coarse.candidates(k).f(2) + (-20:5:20)]; %#ok<AGROW>
end
f1 = unique(f1(f1>=min(cfg.f1Grid)&f1<=max(cfg.f1Grid)));
f2 = unique(f2(f2>=min(cfg.f2Grid)&f2<=max(cfg.f2Grid)));
if isempty(f1) || isempty(f2)
    fit = coarse; return;
end
cfg1 = cfg0; cfg1.f1Grid=f1; cfg1.f2Grid=f2;
fine = fit_vp_main_base(t,V,x,templateLib,gHat,cfg1,[],[]);
allCand = [coarse.candidates(:); fine.candidates(:)];
[~,ord] = sort([allCand.linear_sse],'ascend');
allCand = allCand(ord);
fit = fine;
fit.candidates = allCand(1:min(numel(allCand),max(cfg.numVarproCandidates,10)));
fit = fit.candidates(1);
fit.candidates = allCand(1:min(numel(allCand),max(cfg.numVarproCandidates,10)));
fit.J1 = fit.candidates(1).linear_sse;
fit.J2 = inf;
if numel(fit.candidates)>=2, fit.J2=fit.candidates(2).linear_sse; end
fit.deltaJ=fit.J2-fit.J1; fit.rhoJ=fit.deltaJ/max(fit.J1,eps);
end
