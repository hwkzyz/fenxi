function [S,T]=V4_HighSpeedClosedLoop_Framework_20251222(A,B,methodId)
% Model-independent V4 closure metrics for two independently generated
% high-speed result tables with identical window ordering.
required={'window_id','gap_EO','gap_frequency_hz','gap_amplitude_mm','gap_rmse_mV','gap_mean_delta_gap_mm'};
if isempty(A)||isempty(B)
    T=table; S=table(0,0,NaN,NaN,NaN,NaN,string(methodId),'VariableNames',{'nWindows','nEOMatches','meanAbsFrequencyDiffHz','meanAbsAmplitudeDiffMm','meanAbsRMSEDiffmV','meanAbsDeltaGapDiffMm','methodId'}); return;
end
assert(all(ismember(required,A.Properties.VariableNames))&&all(ismember(required,B.Properties.VariableNames)), ...
    'V4:MissingColumns','Both V4 tables must contain the formal Main10 metric columns.');
assert(height(A)==height(B),'V4:WindowCountMismatch','V4 tables have different window counts (%d versus %d).',height(A),height(B));
assert(isequal(A.window_id,B.window_id),'V4:WindowIdMismatch','V4 tables do not contain identical ordered window IDs.');
n=height(A); rows=zeros(n,9);
for i=1:n
 rows(i,:)=[A.window_id(i),A.gap_EO(i)==B.gap_EO(i),A.gap_frequency_hz(i)-B.gap_frequency_hz(i),A.gap_amplitude_mm(i)-B.gap_amplitude_mm(i),A.gap_rmse_mV(i)-B.gap_rmse_mV(i),A.gap_mean_delta_gap_mm(i)-B.gap_mean_delta_gap_mm(i),A.gap_EO(i),B.gap_EO(i),A.gap_frequency_hz(i)];
end
T=array2table(rows,'VariableNames',{'windowId','EO_match','frequencyDiffHz_AminusB','amplitudeDiffMm_AminusB','RMSEDiffmV_AminusB','deltaGapDiffMm_AminusB','EO_A','EO_B','frequencyAHz'});
S=table(height(rows),sum(rows(:,2)),mean(abs(rows(:,3)),'omitnan'),mean(abs(rows(:,4)),'omitnan'),mean(abs(rows(:,5)),'omitnan'),mean(abs(rows(:,6)),'omitnan'),'VariableNames',{'nWindows','nEOMatches','meanAbsFrequencyDiffHz','meanAbsAmplitudeDiffMm','meanAbsRMSEDiffmV','meanAbsDeltaGapDiffMm'});
S.methodId=string(methodId);
end
