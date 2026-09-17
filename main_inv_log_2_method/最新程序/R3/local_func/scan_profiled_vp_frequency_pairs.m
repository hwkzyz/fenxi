function scan = scan_profiled_vp_frequency_pairs(t,V,x,templateLib,cfg,gGrid,dxGrid)
%SCAN_PROFILED_VP_FREQUENCY_PAIRS  Frequency scan profiled over g and dx.
% The displacement basis includes an explicit Fx*dx column, reducing the
% tendency of the fixed-point VP score to absorb offset error into vibration.
f1 = cfg.f1Grid(:).'; f2 = cfg.f2Grid(:).';
keepN = max(get_opt(cfg,'numVarproCandidates',10)*8,40);
rows = repmat(struct('f1',NaN,'f2',NaN,'sse',Inf,'g_seed',NaN,'dx_seed',NaN), ...
    numel(f1)*numel(f2),1); n=0;
for ig=1:numel(gGrid)
 for id=1:numel(dxGrid)
  F0=eval_gap_template(templateLib,gGrid(ig),x-dxGrid(id));
  Fx=eval_gap_derivative(templateLib,gGrid(ig),x-dxGrid(id));
  ok=isfinite(F0)&isfinite(Fx)&isfinite(V); if nnz(ok)<20, continue; end
  tw=t(ok); yw=V(ok)-F0(ok); fw=Fx(ok);
  for i1=1:numel(f1)
   s1=sin(2*pi*f1(i1)*tw); c1=cos(2*pi*f1(i1)*tw);
   for i2=1:numel(f2)
    if abs(f2(i2)-f1(i1))<20, continue; end
    B=[ones(size(tw)),fw,fw.*s1,fw.*c1,fw.*sin(2*pi*f2(i2)*tw),fw.*cos(2*pi*f2(i2)*tw)];
    beta=B\yw; rr=yw-B*beta; sse=dot(rr,rr);
    n=n+1; rows(n)=struct('f1',f1(i1),'f2',f2(i2),'sse',sse,'g_seed',gGrid(ig),'dx_seed',dxGrid(id));
   end
  end
 end
end
if n==0, scan=struct('candidates',rows([]),'all_rows',rows([])); return; end
rows=rows(1:n); [~,ord]=sort([rows.sse]); rows=rows(ord);
% Keep separated frequency pairs, as in the production VP scanner.
out=repmat(rows(1),0,1);
for i=1:numel(rows)
 if any(abs([out.f1]-rows(i).f1)<8 & abs([out.f2]-rows(i).f2)<8), continue; end
 out(end+1)=rows(i); %#ok<AGROW>
 if numel(out)>=keepN, break; end
end
scan=struct('candidates',out,'all_rows',rows);
end

function v=get_opt(S,n,d)
if isstruct(S)&&isfield(S,n)&&~isempty(S.(n)),v=S.(n);else,v=d;end
end
