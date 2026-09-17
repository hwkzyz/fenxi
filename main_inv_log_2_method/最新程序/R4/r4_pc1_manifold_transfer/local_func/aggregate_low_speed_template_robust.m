function low=aggregate_low_speed_template_robust(dataLow,alpha_k,R_tip,domain,xGrid,opts)
%AGGREGATE_LOW_SPEED_TEMPLATE_ROBUST Step04-like low-speed aggregation.
% Uses per-turn spatial interpolation, iterative center alignment, robust
% turn rejection, bin medians, and a smooth final reference template.
if nargin<6||isempty(opts),opts=struct();end
xGrid=xGrid(:);mapped=map_highspeed_to_space(dataLow,alpha_k,R_tip,domain,0);
keys=[mapped.rev_v(:),mapped.S_v(:)];
if isempty(keys),error('No low-speed samples available.');end
groupKeys=unique(keys,'rows','stable');Y=nan(size(groupKeys,1),numel(xGrid));coverage=zeros(size(Y,1),1);
for ig=1:size(groupKeys,1)
    idx=mapped.rev_v==groupKeys(ig,1)&mapped.S_v==groupKeys(ig,2);
    x=mapped.x_v(idx);v=mapped.V_a(idx);[x,ord]=sort(x(:));v=v(ord);[x,keep]=unique(x,'stable');v=v(keep);
    if numel(x)<2,continue;end
    inside=xGrid>=x(1)&xGrid<=x(end);coverage(ig)=mean(inside);
    Y(ig,inside)=interp1(x,v,xGrid(inside),'pchip');
end
validRows=coverage>=get_opt(opts,'minCoverage',.90);
if nnz(validRows)<3,validRows=coverage>0;end
% Iterative spatial registration against the provisional pointwise median.
ref=median(Y(validRows,:),1,'omitnan').';
shifts=zeros(size(Y,1),1);alignRmse=inf(size(Y,1),1);
shiftGrid=get_opt(opts,'shiftGrid',(-.08:.01:.08));
for ig=find(validRows).'
    y=Y(ig,:).';ok=isfinite(y)&isfinite(ref);if nnz(ok)<20,validRows(ig)=false;continue;end
    best=inf;bestShift=0;
    for ds=shiftGrid(:).'
        ys=interp1(xGrid(ok),y(ok),xGrid-ds,'pchip',NaN);q=isfinite(ref)&isfinite(ys);
        if nnz(q)<20,continue;end
        e=sqrt(mean((ys(q)-ref(q)).^2));
        if e<best,best=e;bestShift=ds;end
    end
    shifts(ig)=bestShift;alignRmse(ig)=best;
end
% Reject outlying turns by robust residual and keep adequate coverage.
q=alignRmse(validRows);med=median(q,'omitnan');sc=1.4826*median(abs(q-med),'omitnan');
if ~isfinite(sc)||sc<eps,sc=max(med*0.05,1e-9);end
keepRows=validRows & alignRmse<=med+4*sc;
if nnz(keepRows)<3,keepRows=validRows;end
Yal=nan(size(Y));
for ig=find(keepRows).'
    row=Y(ig,:).';ok=isfinite(row);
    if nnz(ok)>=2
        Yal(ig,:)=interp1(xGrid(ok),row(ok),xGrid-shifts(ig),'pchip',NaN).';
    end
end
template=median(Yal(keepRows,:),1,'omitnan').';
template=fillmissing(template,'linear','EndValues','nearest');
span=get_opt(opts,'smoothSpan',9);if mod(span,2)==0,span=span+1;end
if span>=3
    template=smoothdata(template,'movmedian',span,'omitnan');
    template=smoothdata(template,'sgolay',span,'omitnan');
end
template=fillmissing(template,'linear','EndValues','nearest');
low=struct('xGrid',xGrid,'templateLow',template,'stdLow',std(Yal(keepRows,:),0,1,'omitnan').',...
    'countLow',sum(isfinite(Yal(keepRows,:)),1).','perTurn',Yal,'mapped',mapped,...
    'num_turns',nnz(keepRows),'selected_rows',keepRows,'alignment_shift',shifts,...
    'alignment_rmse',alignRmse,'coverage',coverage,'raw_perTurn',Y);
end

function v=get_opt(s,n,d)
if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
