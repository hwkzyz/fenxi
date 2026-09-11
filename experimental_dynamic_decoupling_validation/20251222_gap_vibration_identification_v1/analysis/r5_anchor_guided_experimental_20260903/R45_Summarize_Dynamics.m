function [T,Ridge]=R45_Summarize_Dynamics()
% Require all 16 complete runs and exact per-window measured bundle alignment.
root=fileparts(mfilename('fullpath')); pkg=fileparts(fileparts(root));
out=fullfile(root,'results','unified_closure_20260905'); rows={}; nr=0; candidateRows={};
labels={'p1','p2','p3','p4','p5','oracle','commonb2','r5lobo'};
for blade=[1 5]
 if blade==1, suffix='r5_control_b1_contract_20260904'; else, suffix='r5_control_b5_contract_20260904'; end
 base=load_result(pkg,blade,suffix); assert(numel(base.WindowResult)==18);
 libBase=load(fullfile(pkg,'inputs','calibration',sprintf('Step06I_OffsetTiltShared_GapLibrary_20251222_B%d_S123.mat',blade)),'CorrectedGapLibrary');
 for k=1:numel(labels)
  label=labels{k}; R=load_result(pkg,blade,sprintf('r45_final_b%d_%s_20260905',blade,label));
  assert(R.cfg.targetBlade==blade && numel(R.WindowResult)==18);
  assert(isequal(R.cfg.analysisSensors,[1 2 3]));
  assert(R.ZeroGapEquivalenceGlobalMaxAbsMv<1e-9);
  L=load(R.correctedLibFile,'CorrectedGapLibrary');
  assert(isequaln(L.CorrectedGapLibrary.sensor,libBase.CorrectedGapLibrary.sensor));
  for w=1:18
   W=R.WindowResult(w); B=base.WindowResult(w);
   assert(isequaln(W.timeWindow,B.timeWindow));
   for field={'X','T','V','sensorIndex','sensorIds','lapRange'}
    f=field{1}; assert(isequaln(W.bundle.(f),B.bundle.(f)),'Measured bundle differs: %s',f);
   end
   f=W.modelFits.gap_only; b=B.modelFits.gap_only;
   c=f.CandidateTable;
   c=addvars(c,repmat(blade,height(c),1),repmat(string(label),height(c),1),repmat(w,height(c),1), ...
    'Before',1,'NewVariableNames',{'blade','surface_candidate','window'});
   candidateRows{end+1}=c; %#ok<AGROW>
   dg=double(f.deltaGapMm(:)); g0=double(W.bundle.g0BySensor(:));
   assert(numel(dg)==3 && numel(g0)==3);
   nr=nr+1; rows{nr}=table(blade,string(label),w,mean(W.timeWindow),f.EO,f.freqHz, ...
    f.amplitudeMm,f.dxMm,f.plainRmseMv,dg(1),dg(2),dg(3), ...
    g0(1),g0(2),g0(3),g0(1)+dg(1),g0(2)+dg(2),g0(3)+dg(3), ...
    b.EO,b.freqHz,f.EO==b.EO,f.freqHz-b.freqHz,f.plainRmseMv-b.plainRmseMv, ...
    'VariableNames',{'blade','candidate','window','center_sec','eo','frequency_hz','amplitude_mm', ...
    'dx_mm','rmse_mv','dg1_mm','dg2_mm','dg3_mm','g01_mm','g02_mm','g03_mm', ...
    'g_total1_mm','g_total2_mm','g_total3_mm','control_eo','control_frequency_hz', ...
    'eo_matches_control','delta_frequency_hz','delta_rmse_mv'});
  end
 end
end
T=vertcat(rows{:}); writetable(T,fullfile(out,'Dynamics_AllWindows.csv'));
CandidateCosts=vertcat(candidateRows{:}); writetable(CandidateCosts,fullfile(out,'Dynamics_AllEOCandidateCosts.csv'));
Summary=groupsummary(T,{'blade','candidate'},{'mean','min','max'}, ...
 {'frequency_hz','amplitude_mm','rmse_mv','delta_frequency_hz','eo_matches_control'});
writetable(Summary,fullfile(out,'Dynamics_Summary.csv'));
rr={}; n=0;
for blade=[1 5]
 for w=1:18
  A=T(T.blade==blade & T.window==w & startsWith(T.candidate,'p'),:);
  assert(height(A)==5); n=n+1;
  rr{n}=table(blade,w,numel(unique(A.eo)),range(A.frequency_hz),range(A.amplitude_mm), ...
   range(A.dx_mm),range(A.rmse_mv),range(A.dg1_mm),range(A.dg2_mm),range(A.dg3_mm), ...
   'VariableNames',{'blade','window','eo_count','frequency_span_hz','amplitude_span_mm', ...
   'dx_span_mm','rmse_span_mv','dg1_span_mm','dg2_span_mm','dg3_span_mm'});
 end
end
Ridge=vertcat(rr{:}); writetable(Ridge,fullfile(out,'Ridge_PerWindowSpans.csv'));
save(fullfile(out,'Dynamics_Closure.mat'),'T','Summary','Ridge','-v7.3');
disp(groupsummary(Ridge,'blade','max',{'frequency_span_hz','amplitude_span_mm','eo_count'}));
end
function R=load_result(pkg,blade,suffix)
f=fullfile(pkg,'results','gap_aware',sprintf('Main_GapAware_VPTopK_FullWave_20251222_B%d_S123_%s.mat',blade,suffix));
S=load(f,'Result'); R=S.Result;
end
