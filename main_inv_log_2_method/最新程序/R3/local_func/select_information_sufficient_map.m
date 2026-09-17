function [selectedMap, info] = select_information_sufficient_map(map, lib, gRef, cfg)
%SELECT_INFORMATION_SUFFICIENT_MAP Deterministic information-preserving subset.
% A reference Jacobian dictionary covers gap, translation, and a coarse
% Fourier band. The smallest stratified subset whose normalized information
% matrix spectrally approximates the complete trusted set is returned.
n=numel(map.t_v);
if n<2,selectedMap=map;info=struct('mode',"complete",'used',n);return;end
t=map.t_v(:);x=map.x_v(:);
h=get_field(cfg,'informationGapStepMm',1e-3);
Fx=eval_gap_derivative(lib,gRef,x);
Fg=(eval_gap_template(lib,gRef+h,x)-eval_gap_template(lib,gRef-h,x))/(2*h);
anchors=get_field(cfg,'informationFrequencyAnchorsHz',300:300:1500);
D=[Fg,-Fx];
for f=anchors(:).'
    D=[D,-Fx.*sin(2*pi*f*t),-Fx.*cos(2*pi*f*t)]; %#ok<AGROW>
end
valid=all(isfinite(D),2);
if ~all(valid)
    map=filter_map(map,valid);D=D(valid,:);n=size(D,1);
end
scale=sqrt(sum(D.^2,1));keep=scale>max(scale)*1e-10;
D=D(:,keep)./max(scale(keep),eps)*sqrt(n);
[~,S,V]=svd(D,'econ');sv=diag(S);
r=nnz(sv>max(sv)*get_field(cfg,'informationRankTolerance',1e-8));
if r==0,error('invlog2:EmptyInformationDictionary','No informative trusted samples.');end
Vr=V(:,1:r);lambda=(sv(1:r).^2)/n;
W=Vr*diag(1./sqrt(lambda));
groupCount=count_groups(map);
oversample=get_field(cfg,'informationOversampleFactor',8);
minPerGroup=get_field(cfg,'informationMinPerGroup',4);
nMin=max([ceil(oversample*r),minPerGroup*max(groupCount,1),r+1]);
nMin=min(n,nMin);
minInfo=get_field(cfg,'informationMinEigenvalue',0.35);
maxInfo=get_field(cfg,'informationMaxEigenvalue',2.5);
growth=get_field(cfg,'informationGrowthFactor',1.5);
sizes=nMin;
while sizes(end)<n
    sizes(end+1)=min(n,max(sizes(end)+1,ceil(growth*sizes(end)))); %#ok<AGROW>
end
bestIdx=(1:n).';bestErr=0;chosen=n;
for k=1:numel(sizes)
    idx=stratified_indices(map,sizes(k),minPerGroup);
    Gs=(D(idx,:).'*D(idx,:))/numel(idx);
    ev=eig((W.'*Gs*W+W.'*Gs.'*W)/2);
    err=max(abs(real(ev)-1));
    if (min(real(ev))>=minInfo&&max(real(ev))<=maxInfo)||sizes(k)==n
        bestIdx=idx;bestErr=err;chosen=numel(idx);break;
    end
end
selectedMap=filter_map(map,ismember((1:n).',bestIdx));
info=struct('mode',"jacobian_spectral_stratified",'full_count',n,...
    'used',chosen,'dictionary_rank',r,'group_count',groupCount,...
    'spectral_error',bestErr,'information_min_eigenvalue',min(real(ev)),...
    'information_max_eigenvalue',max(real(ev)),...
    'information_min_threshold',minInfo,'information_max_threshold',maxInfo,...
    'frequency_anchors_hz',anchors);
selectedMap.sampling_info=info;
end

function idx=stratified_indices(map,targetN,minPerGroup)
n=numel(map.t_v);
if targetN>=n,idx=(1:n).';return;end
if ~isfield(map,'rev_v')||~isfield(map,'S_v')
    idx=unique(round(linspace(1,n,targetN))).';return;
end
[~,~,gid]=unique([map.rev_v(:),map.S_v(:)],'rows','stable');G=max(gid);
counts=accumarray(gid,1,[G,1]);raw=targetN*counts/n;
take=max(minPerGroup,floor(raw));take=min(take,counts);
while sum(take)<targetN
    room=counts-take;gain=raw-take;gain(room<=0)=-Inf;[~,j]=max(gain);
    if room(j)<=0,break;end,take(j)=take(j)+1;
end
while sum(take)>targetN
    eligible=find(take>min(minPerGroup,counts));if isempty(eligible),break;end
    [~,j]=max(take(eligible)-raw(eligible));take(eligible(j))=take(eligible(j))-1;
end
idx=zeros(sum(take),1);at=0;
for g=1:G
    members=find(gid==g);[~,ord]=sortrows([map.x_v(members),map.t_v(members)],[1 2]);
    members=members(ord);q=unique(round(linspace(1,numel(members),take(g))));
    sel=members(q);idx(at+(1:numel(sel)))=sel;at=at+numel(sel);
end
idx=sort(idx(1:at));
end

function out=filter_map(map,mask)
n=numel(map.t_v);out=map;names=fieldnames(map);
for i=1:numel(names)
    v=map.(names{i});
    if (isnumeric(v)||islogical(v)||isstring(v)||iscell(v))&&numel(v)==n&&~isscalar(v)
        out.(names{i})=v(mask);
    end
end
end

function n=count_groups(map)
if isfield(map,'rev_v')&&isfield(map,'S_v')
    n=size(unique([map.rev_v(:),map.S_v(:)],'rows'),1);
else,n=1;end
end

function v=get_field(s,name,d)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)),v=s.(name);else,v=d;end
end
