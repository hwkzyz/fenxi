function [map,C]=restrict_map_to_support_contract(map,C)
%RESTRICT_MAP_TO_SUPPORT_CONTRACT Keep only the certified high-speed domain.
if ~C.certified_domain_nonempty,error('support:EmptyCertifiedDomain','Certified high-speed domain is empty.');end
domain=C.high_certified_mm;
keep=map.x_v(:)>=domain(1)&map.x_v(:)<=domain(2);
if ~any(keep),error('support:NoCertifiedSamples','No high-speed sample lies in the certified domain.');end
n=numel(map.x_v);names=fieldnames(map);
for i=1:numel(names)
    value=map.(names{i});
    if ~isscalar(value)&&numel(value)==n&&(isnumeric(value)||islogical(value)||isstring(value)||iscell(value))
        map.(names{i})=value(keep);
    end
end
used=[min(map.x_v) max(map.x_v)];
C.high_used_mm=used;C.left_margin_mm=used(1)-domain(1);
C.right_margin_mm=domain(2)-used(2);C.used_sample_count=nnz(keep);
C.removed_sample_count=n-nnz(keep);C.runtime_domain_valid=...
    C.left_margin_mm>=-10*eps(max(abs(domain)))&&C.right_margin_mm>=-10*eps(max(abs(domain)));
map.support_contract=C;
end
