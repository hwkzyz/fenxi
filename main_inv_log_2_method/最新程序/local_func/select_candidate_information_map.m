function [out,info]=select_candidate_information_map(map,fit,lib,targetN,cfg)
%SELECT_CANDIDATE_INFORMATION_MAP Amplitude-aware c-optimal refinement set.
% The pilot fit supplies amplitudes and frequencies. Rows are scored by
% their contribution to the covariance of the frequency parameters, while
% sensor-by-turn coverage remains mandatory.
n=numel(map.t_v);targetN=min(n,max(1,round(targetN)));
t=map.t_v(:);x=map.x_v(:);t0=mean(t);
if isfield(fit,'fit'),q=fit.fit;else,q=fit;end
if isfield(fit,'g_used'),g=fit.g_used;elseif isfield(q,'g'),g=q.g;else,g=NaN;end
if isfield(fit,'dx_used'),dx=fit.dx_used;elseif isfield(q,'dx'),dx=q.dx;else,dx=0;end
p=q.p;u=fit_u(p,t);z=x-dx-u;Fx=eval_gap_derivative(lib,g,z);h=1e-3;
Fg=(eval_gap_template(lib,g+h,z)-eval_gap_template(lib,g-h,z))/(2*h);
if numel(p)>=6
    s1=sin(2*pi*p(3)*t+p(2));c1=cos(2*pi*p(3)*t+p(2));
    s2=sin(2*pi*p(6)*t+p(5));c2=cos(2*pi*p(6)*t+p(5));
    J=[Fg,-Fx,-Fx.*s1,-Fx.*p(1).*c1,-Fx.*p(1).*c1.*(2*pi*(t-t0)),...
        -Fx.*s2,-Fx.*p(4).*c2,-Fx.*p(4).*c2.*(2*pi*(t-t0))];
    fcols=[5 8];amplitudes=[p(1) p(4)];
else
    s=sin(2*pi*p(3)*t+p(2));c=cos(2*pi*p(3)*t+p(2));
    J=[Fg,-Fx,-Fx.*s,-Fx.*p(1).*c,-Fx.*p(1).*c.*(2*pi*(t-t0))];
    fcols=5;amplitudes=p(1);
end
valid=all(isfinite(J),2);J(~valid,:)=0;
C=pinv(J.'*J);targetDirection=J*C(:,fcols);
frequencyScore=sum(targetDirection.^2,2);
[U,~,~]=svd(J,'econ');generalScore=sum(U.^2,2);
score=frequencyScore/max(max(frequencyScore),eps)+...
    get_field(cfg,'candidateGeneralLeverageWeight',0.2)*...
    generalScore/max(max(generalScore),eps);
idx=select_stratified_score(map,score,targetN,...
    get_field(cfg,'informationMinPerGroup',4));
mask=false(n,1);mask(idx)=true;out=filter_map(map,mask);
info=struct('mode',"amplitude_aware_frequency_c_optimal",'full_count',n,...
    'used',numel(idx),'amplitudes_mm',amplitudes,...
    'minimum_amplitude_mm',min(abs(amplitudes)),'frequency_columns',fcols);
out.sampling_info=info;
end

function idx=select_stratified_score(map,score,targetN,minPerGroup)
n=numel(score);if targetN>=n,idx=(1:n).';return;end
if ~isfield(map,'rev_v')||~isfield(map,'S_v')
    [~,ord]=sort(score,'descend');idx=sort(ord(1:targetN));return;
end
[~,~,gid]=unique([map.rev_v(:),map.S_v(:)],'rows','stable');G=max(gid);
counts=accumarray(gid,1,[G,1]);raw=targetN*counts/n;
take=max(minPerGroup,floor(raw));take=min(take,counts);
while sum(take)<targetN
    room=counts-take;gain=raw-take;gain(room<=0)=-Inf;[~,j]=max(gain);
    if room(j)<=0,break;end,take(j)=take(j)+1;
end
while sum(take)>targetN
    e=find(take>min(minPerGroup,counts));if isempty(e),break;end
    [~,j]=max(take(e)-raw(e));take(e(j))=take(e(j))-1;
end
idx=zeros(sum(take),1);at=0;
for g=1:G
    m=find(gid==g);[~,ord]=sort(score(m),'descend');sel=m(ord(1:take(g)));
    idx(at+(1:numel(sel)))=sel;at=at+numel(sel);
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

function v=get_field(s,name,d)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)),v=s.(name);else,v=d;end
end
