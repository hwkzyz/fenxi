function C=frequency_solution_confidence(bestFrequency,bestRmse,alternativeFrequency,alternativeRmse,nSamples,cfg)
%FREQUENCY_SOLUTION_CONFIDENCE Compare the selected voltage minimum with a distinct basin.
bestFrequency=sort(bestFrequency(:).');
alternativeRmse=alternativeRmse(:);
if isempty(alternativeFrequency)
    alternativeFrequency=zeros(0,numel(bestFrequency));
elseif isvector(alternativeFrequency)&&numel(bestFrequency)==1
    alternativeFrequency=alternativeFrequency(:);
end
tol=get_field(cfg,'identificationDistinctFrequencyHz',2);
valid=isfinite(alternativeRmse)&all(isfinite(alternativeFrequency),2);
distinct=false(size(valid));
for i=find(valid).'
    f=sort(alternativeFrequency(i,:));
    distinct(i)=numel(f)~=numel(bestFrequency)||max(abs(f-bestFrequency))>tol;
end
idx=find(valid&distinct);
if isempty(idx)
    competitor=nan(1,numel(bestFrequency));compRmse=NaN;
    rmseMargin=NaN;mseMargin=NaN;zMargin=NaN;confident=false;status="unvalidated";
else
    [compRmse,j]=min(alternativeRmse(idx));k=idx(j);
    competitor=sort(alternativeFrequency(k,:));
    rmseMargin=compRmse-bestRmse;
    mseMargin=compRmse^2-bestRmse^2;
    p=get_field(cfg,'identificationEffectiveParameterCount',numel(bestFrequency)*3+2);
    noiseMseStd=max(bestRmse^2,eps)*sqrt(2/max(nSamples-p,1));
    zMargin=mseMargin/noiseMseStd;
    confident=zMargin>=get_field(cfg,'identificationConfidenceSigma',2);
    if confident,status="identified";else,status="ambiguous";end
end
C=struct('frequency_margin',rmseMargin,'mse_margin',mseMargin,...
    'noise_normalized_margin',zMargin,'identification_confident',confident,...
    'identification_status',status,...
    'competing_frequency_or_order',competitor,'competing_rmse',compRmse,...
    'distinct_frequency_tolerance_hz',tol);
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
