function T = R5_Evaluate_StaticMigration(S, mode, outputDir)
%R5_EVALUATE_STATICMIGRATION Static evidence for R5 migration operators.
% mode=same_blade excludes the direct anchor->target operator.
% mode=leave_gap removes the target state from every donor and interpolates
% its increment using only the remaining state grid.
arguments
    S (1,1) struct
    mode (1,:) char {mustBeMember(mode,{'same_blade','leave_gap'})}
    outputDir (1,:) char = ''
end
nB = numel(S.blade_id); nG = size(S.waveforms_mv,3);
stateValues = 1:nG;
L = R5_Build_TransitionLibrary(S,'effective_state',stateValues);
rows = cell(nB*(nG-1),1); ir = 0;
savedWaveform = false;
for ib = 1:nB
    for it = 2:nG
        anchor = struct('x_mm',S.x_mm,'y_mv',S.waveforms_mv(:,ib,1), ...
            'supportMask',S.support_mask);
        keep = true(numel(L),1);
        if strcmp(mode,'same_blade')
            keep = [L.sourceBladeId].' == S.blade_id(ib) & ...
                ~([L.anchorState].' == 1 & [L.targetState].' == it);
        end
        try
            if strcmp(mode,'leave_gap')
                [yhat, support, status] = migrate_holdout_gap(anchor,S,it);
            else
                C = R5_Generate_CandidateStates(anchor,L(keep),it,struct('topKDonor',3));
                yhat = C.waveform_mv; support = C.supportMask; status = C.status;
            end
            ytrue = double(S.waveforms_mv(:,ib,it));
            mask = support & isfinite(yhat) & isfinite(ytrue);
            if nnz(mask) < 3
                rmse=NaN; rel=NaN; xp=NaN; xt=NaN; status='insufficient_holdout_support';
            else
                rmse = sqrt(mean((yhat(mask)-ytrue(mask)).^2));
                rel = rmse / max(range(ytrue(mask)),eps);
                [~,ip] = max(yhat(mask)); xx = S.x_mm(mask); xp = xx(ip);
                [~,itp] = max(ytrue(mask)); xt = xx(itp);
                if ~isempty(outputDir) && ~savedWaveform
                    figDir=fullfile(outputDir,'figures'); if ~exist(figDir,'dir'),mkdir(figDir);end
                    R5_Plot_StaticWaveformComparison(S.x_mm,anchor.y_mv,ytrue,yhat, ...
                        sprintf('R5 %s: B%d, held state %d',strrep(mode,'_',' '),S.blade_id(ib),it), ...
                        fullfile(figDir,['R5_',mode,'_waveform_comparison.png']));
                    savedWaveform=true;
                end
            end
        catch ME
            rmse=NaN; rel=NaN; xp=NaN; xt=NaN; status=['error:',ME.identifier]; mask=false(size(S.x_mm));
        end
        ir=ir+1;
        rows{ir}=table(S.blade_id(ib),it,S.physical_gap_mm(it),nnz(mask),rmse,rel,xp-xt,string(status), ...
            'VariableNames',{'blade_id','target_state','nominal_gap_mm','n_support','rmse_mv','relative_rmse','peak_x_error_mm','status'});
    end
end

function [yhat, support, status] = migrate_holdout_gap(anchor, S, targetState)
x = S.x_mm(:); states = S.physical_gap_mm(:).'; keep = true(size(states)); keep(targetState)=false;
deltas = nan(numel(x),numel(S.blade_id));
for ib=1:numel(S.blade_id)
    Y = double(squeeze(S.waveforms_mv(:,ib,:)));
    target = interp1(states(keep),Y(:,keep).',states(targetState),'pchip',NaN).';
    target = target(:);
    deltas(:,ib) = target - Y(:,1);
end
delta=mean(deltas,2,'omitnan'); yhat=anchor.y_mv(:)+delta;
support=anchor.supportMask(:) & isfinite(delta); status='holdout_gap_interpolation';
end
T=vertcat(rows{1:ir});
if ~isempty(outputDir)
    if ~exist(outputDir,'dir'), mkdir(outputDir); end
    writetable(T,fullfile(outputDir,['R5_',mode,'_migration_metrics.csv']));
    R5_Plot_StaticMigrationSummary(T,mode,fullfile(outputDir,['R5_',mode,'_migration_summary.png']));
end
end
