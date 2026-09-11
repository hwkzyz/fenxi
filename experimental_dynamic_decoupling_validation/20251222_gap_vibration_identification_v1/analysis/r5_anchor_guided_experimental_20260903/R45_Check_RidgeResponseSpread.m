function T=R45_Check_RidgeResponseSpread()
% Evaluate every candidate at the SAME control trajectory and gap increments.
% This distinguishes nontrivial observation-model changes from fitted
% parameter compensation. Formula and clipping follow unchanged Main10.
root=fileparts(mfilename('fullpath')); pkg=fileparts(fileparts(root));
out=fullfile(root,'results','unified_closure_20260905'); rows={}; n=0;
for blade=[1 5]
 if blade==1, suffix='r5_control_b1_contract_20260904'; else, suffix='r5_control_b5_contract_20260904'; end
 S=load(fullfile(pkg,'results','gap_aware',sprintf('Main_GapAware_VPTopK_FullWave_20251222_B%d_S123_%s.mat',blade,suffix)),'Result');
 W=S.Result.WindowResult;
 libraries=cell(1,5);
 for k=1:5
  L=load(fullfile(out,'ridge_sidecars',sprintf('R5_RidgeSidecar_B%d_P%d.mat',blade,k)),'CorrectedGapLibrary'); libraries{k}=L.CorrectedGapLibrary;
 end
 for w=1:numel(W)
  fit=W(w).modelFits.gap_only; bundle=W(w).bundle;
  delta=zeros(numel(bundle.X),5);
  for k=1:5
   L=libraries{k}; R=L.responseSurface;
   for is=1:numel(bundle.sensorIds)
    sid=bundle.sensorIds(is); c=L.sensor([L.sensor.sensorId]==sid); m=bundle.sensorIndex==is;
    x=bundle.X(m)-fit.dxMm-fit.uMm(m)-fit.sensorEtaMm(is);
    delta(m,k)=c.voltageGain*(raw(R,c,x,fit.deltaGapMm(is))-raw(R,c,x,0));
   end
  end
  assert(all(isfinite(delta(:))));
  rmsByCandidate=sqrt(mean(delta.^2,1)); maxPair=0;
  for i=1:5
   for j=i+1:5, maxPair=max(maxPair,sqrt(mean((delta(:,i)-delta(:,j)).^2))); end
  end
  n=n+1; rows{n}=table(blade,w,min(rmsByCandidate),max(rmsByCandidate),maxPair, ...
   'VariableNames',{'blade','window','min_increment_rms_mv','max_increment_rms_mv','max_pair_increment_difference_rms_mv'});
 end
end
T=vertcat(rows{:}); writetable(T,fullfile(out,'Ridge_FixedTrajectoryResponseSpread.csv')); disp(groupsummary(T,'blade','max','max_pair_increment_difference_rms_mv'));
end
function v=raw(R,c,x,dg)
xq=c.xScale*(x-c.tauMm); g=c.g0Mm+dg+c.muGapPerXMm*(x-c.tauMm);
xq=min(max(xq,min(R.xGrid)),max(R.xGrid)); g=min(max(g,min(R.gTrainMm)),max(R.gTrainMm));
B=interp1(R.xGrid,R.coeff,xq,'linear'); v=B(:,1)+B(:,2)./g+B(:,3).*log(g/R.g0Mm);
end
