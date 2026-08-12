function [Detail,Competition]=Run_DualSyncFunnelV2ConfidenceAnalysis(seeds,snrLevels)
%RUN_DUALSYNCFUNNELV2CONFIDENCEANALYSIS Audit final-basin competition.
% This intentionally targets the three pressure cases where Funnel V2 has
% already retained the true EO pair in the final three candidates.
if nargin<1||isempty(seeds),seeds=1:20;end
if nargin<2||isempty(snrLevels),snrLevels=[25 15 5];end

root=fileparts(mfilename('fullpath'));
outDir=fullfile(root,'output','dual_sync_funnel_v2_generalization_daq_snr');
ids=["C02_eo_10_12";"X01_corner_near";"X02_corner_lowamp"];
Detail=table();
for snr=snrLevels(:).'
    [D,~]=Run_DualSyncFunnelV1Generalization(snr,seeds,ids,false,40,...
        "v2_all_joint1",24,"confidence_stableseed");
    Detail=[Detail;D]; %#ok<AGROW>
end

Competition=summarize_competition(Detail,ids,snrLevels);
writetable(Detail,fullfile(outDir,'detail_final_competition_v2.csv'));
writetable(Competition,fullfile(outDir,'summary_final_competition_v2.csv'));
plot_competition(Detail,Competition,outDir);
save(fullfile(outDir,'final_competition_v2.mat'),'Detail','Competition','ids','snrLevels','seeds');
disp(Competition);
end

function S=summarize_competition(D,ids,snrLevels)
n=numel(ids)*numel(snrLevels);
row=struct('nominal_snr_db',NaN,'condition_id',"",'n',0,...
    'eo_correct_rate',NaN,'recall_at_3_rate',NaN,'confidence_rate',NaN,...
    'ambiguous_rate',NaN,'correct_and_confident_rate',NaN,...
    'wrong_and_confident_rate',NaN,'median_z_all',NaN,'median_z_correct',NaN,...
    'median_z_wrong',NaN,'median_frequency_margin_V',NaN,...
    'most_common_competitor_eo1',NaN,'most_common_competitor_eo2',NaN);
rows=repmat(row,n,1);k=0;
for snr=snrLevels(:).'
    for i=1:numel(ids)
        k=k+1;q=D(D.nominal_snr_db==snr&D.condition_id==ids(i),:);
        rows(k).nominal_snr_db=snr;rows(k).condition_id=ids(i);rows(k).n=height(q);
        if isempty(q),continue;end
        correct=logical(q.eo_correct);
        rows(k).eo_correct_rate=mean(correct);
        rows(k).recall_at_3_rate=mean(q.recall_at_3);
        confident=logical(q.identification_confident);
        rows(k).confidence_rate=mean(confident);
        rows(k).ambiguous_rate=mean(q.identification_status=="ambiguous");
        rows(k).correct_and_confident_rate=mean(correct&confident);
        rows(k).wrong_and_confident_rate=mean(~correct&confident);
        rows(k).median_z_all=median(q.noise_normalized_margin,'omitnan');
        rows(k).median_z_correct=median(q.noise_normalized_margin(correct),'omitnan');
        rows(k).median_z_wrong=median(q.noise_normalized_margin(~correct),'omitnan');
        rows(k).median_frequency_margin_V=median(q.frequency_margin_V,'omitnan');
        [rows(k).most_common_competitor_eo1,rows(k).most_common_competitor_eo2]=mode_pair(q);
    end
end
S=struct2table(rows);
end

function [a,b]=mode_pair(q)
a=NaN;b=NaN;p=[q.competitor_eo1 q.competitor_eo2];
p=round(p,6);p=p(all(isfinite(p),2),:);if isempty(p),return;end
[u,~,ix]=unique(p,'rows');[~,j]=max(accumarray(ix,1));a=u(j,1);b=u(j,2);
end

function plot_competition(D,S,outDir)
labels=string(S.nominal_snr_db)+" dB | "+S.condition_id;
f=figure('Color','w','Units','centimeters','Position',[1 1 18 12.5]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');x=1:height(S);
nexttile;bar(x,[S.eo_correct_rate S.confidence_rate S.recall_at_3_rate],'grouped');
ylim([0 1.05]);ylabel('Rate');legend({'Final EO correct','Confidence flag','True EO in Top-3'},'Location','southwest');title('(a) Final selection and coverage');
nexttile;bar(x,[S.correct_and_confident_rate S.wrong_and_confident_rate S.ambiguous_rate],'stacked');
ylim([0 1.05]);ylabel('Rate');legend({'Correct and confident','Wrong but confident','Ambiguous'},'Location','best');title('(b) Confidence outcome');
nexttile;hold on;
plot(x,S.median_z_all,'ko-','LineWidth',1.0,'MarkerFaceColor','w');yline(2,'--r','2 sigma');
ylabel('Noise-normalized margin');title('(c) Median basin separation');
nexttile;bar(x,S.median_frequency_margin_V);ylabel('RMSE margin (V)');title('(d) Final-basin RMSE separation');
ax=findall(f,'Type','axes');for k=1:numel(ax),set(ax(k),'FontName','Times New Roman','FontSize',8,...
        'TickDir','in','Box','on','XTick',x,'XTickLabel',labels,'XTickLabelRotation',45);end
xlabel(tl,'Noise level and pressure case','FontName','Times New Roman','FontSize',9);
exportgraphics(f,fullfile(outDir,'final_competition_v2.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'final_competition_v2.pdf'),'ContentType','vector');
close(f);
end
