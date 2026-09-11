function T = R5_V4_LeaveBlade(S, outputDir)
%R5_V4_LEAVEBLADE Cross-blade transfer diagnostic; never validates local gap.
arguments, S (1,1) struct; outputDir (1,:) char = ''; end
nB=numel(S.blade_id); nG=size(S.waveforms_mv,3); L=R5_Build_TransitionLibrary(S,'effective_state',1:nG);
rows=cell(nB*(nG-1),1); ir=0;
savedWaveform=false;
for ib=1:nB
  for it=2:nG
    a=struct('x_mm',S.x_mm,'y_mv',S.waveforms_mv(:,ib,1),'supportMask',S.support_mask);
    keep=[L.sourceBladeId].' ~= S.blade_id(ib);
    try
      C=R5_Generate_CandidateStates(a,L(keep),it,struct('topKDonor',3)); y=S.waveforms_mv(:,ib,it); m=C.supportMask & isfinite(y);
      e=sqrt(mean((C.waveform_mv(m)-y(m)).^2)); rel=e/max(range(y(m)),eps); st=C.status;
      if ~isempty(outputDir) && ~savedWaveform && nnz(m)>=3
        figDir=fullfile(outputDir,'figures'); if ~exist(figDir,'dir'),mkdir(figDir);end
        R5_Plot_StaticWaveformComparison(S.x_mm,a.y_mv,y,C.waveform_mv, ...
          sprintf('R5 leave blade: held B%d, state %d',S.blade_id(ib),it), ...
          fullfile(figDir,'R5_leave_blade_waveform_comparison.png'));
        savedWaveform=true;
      end
    catch ME, e=NaN; rel=NaN; m=false(size(S.x_mm)); st=['error:',ME.identifier]; end
    ir=ir+1; rows{ir}=table(S.blade_id(ib),it,nnz(m),e,rel,string(st), ...
      'VariableNames',{'heldout_blade_id','target_state','n_support','rmse_mv','relative_rmse','status'});
  end
end
T=vertcat(rows{1:ir});
if ~isempty(outputDir), if ~exist(outputDir,'dir'),mkdir(outputDir);end; writetable(T,fullfile(outputDir,'R5_V4_leave_blade_metrics.csv')); end
end
