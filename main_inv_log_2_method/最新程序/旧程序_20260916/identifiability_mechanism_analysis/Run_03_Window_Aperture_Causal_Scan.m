function Result=Run_03_Window_Aperture_Causal_Scan()
%RUN_03_WINDOW_APERTURE_CAUSAL_SCAN Test aperture-rank-information linkage.
analysisDir=fileparts(mfilename('fullpath'));
src=fullfile(analysisDir,'output','02_exact_frozen_conditional_geometry', ...
    'exact_frozen_conditional_geometry.mat');
if ~isfile(src)
    error('identifiability:MissingGeometry','Run Run_02 first.');
end
S=load(src,'map','model','cfg','cases');
map=S.map;model=S.model;cfg=S.cfg;cases=S.cases;
halfMax=min(max(abs(map.x_v)),3.5);
halfWidths=linspace(.12,halfMax,18);
rows=repmat(row0(),0,1);
for ic=1:numel(cases)
    c=cases(ic);
    for ih=1:numel(halfWidths)
        use=abs(map.x_v)<=halfWidths(ih);
        for im=1:2
            if im==1,mode="exact";phase=map.theta_v(:);
            else,mode="frozen";phase=map.theta_v(:)-map.x_v(:)/cfg.R_tip;end
            x=map.x_v(:);sid=map.S_v(:);
            u=c.A(1)*sin(c.eo(1)*phase+c.phi(1))+c.A(2)*sin(c.eo(2)*phase+c.phi(2));
            z=x-u;fx=eval_gap_derivative(model,.6,z,sid);h=.002;
            fg=(eval_gap_template(model,.6+h,z,sid)-eval_gap_template(model,.6-h,z,sid))/(2*h);
            q=-fx;
            Jp=[q.*sin(c.eo(1)*phase),q.*cos(c.eo(1)*phase), ...
                q.*sin(c.eo(2)*phase),q.*cos(c.eo(2)*phase)];
            Jn=[-fg,q];valid=use&all(isfinite([Jp Jn]),2);
            M=metrics(Jp(valid,:),Jn(valid,:));
            rr=row0();rr.case_id=c.id;rr.mode=mode;rr.window_half_width_mm=halfWidths(ih);
            rr.window_full_width_mm=2*halfWidths(ih);rr.sample_count=nnz(valid);
            rr.weight_energy=sum(q(valid).^2);rr.rank_conditioned=M.rank_conditioned;
            rr.lambda_min_scaled=M.lambda_min_scaled;rr.condition_scaled=M.condition_scaled;
            rows(end+1,1)=rr; %#ok<AGROW>
        end
    end
end
Scan=struct2table(rows);
Gate=summarize_gate(Scan);
out=fullfile(analysisDir,'output','03_window_aperture_causal_scan');
if ~exist(out,'dir'),mkdir(out);end
writetable(Scan,fullfile(out,'window_scan.csv'));
writetable(Gate,fullfile(out,'gate_summary.csv'));
save(fullfile(out,'window_scan.mat'),'Scan','Gate','halfWidths');
write_report(out,Gate);
disp(Gate);Result=struct('Scan',Scan,'Gate',Gate,'outputDir',out);
end

function M=metrics(Jp,Jn)
Q=orthb(Jn);Jps=Jp./max(vecnorm(Jp,2,1),eps);Jc=Jps-Q*(Q.'*Jps);
e=sort(real(eig((Jc.'*Jc+(Jc.'*Jc).')/2)),'ascend');
M=struct('rank_conditioned',rank(Jc,1e-10*norm(Jc,2)), ...
    'lambda_min_scaled',max(0,e(1)), ...
    'condition_scaled',e(end)/max(e(1),eps));
end
function Q=orthb(A)
[U,S,~]=svd(A,'econ');s=diag(S);r=sum(s>max(size(A))*eps(max(s)));Q=U(:,1:r);
end
function T=summarize_gate(S)
ids=unique(S.case_id,'stable');out=repmat(struct('case_id',"", ...
    'exact_first_full_rank_width_mm',NaN,'frozen_max_rank',NaN, ...
    'exact_information_monotonic_spearman',NaN),numel(ids),1);
for i=1:numel(ids)
    q=S(S.case_id==ids(i),:);e=q(q.mode=="exact",:);f=q(q.mode=="frozen",:);
    j=find(e.rank_conditioned==4,1,'first');if ~isempty(j),w=e.window_full_width_mm(j);else,w=NaN;end
    out(i)=struct('case_id',ids(i),'exact_first_full_rank_width_mm',w, ...
        'frozen_max_rank',max(f.rank_conditioned), ...
        'exact_information_monotonic_spearman', ...
        corr(e.window_full_width_mm,e.lambda_min_scaled,'Type','Spearman'));
end
T=struct2table(out);
end
function write_report(out,G)
fid=fopen(fullfile(out,'README.md'),'w');fprintf(fid,'# Window-aperture causal scan\n\n');
fprintf(fid,'The scan changes the physical window and therefore both phase aperture and sensitivity energy. It establishes a necessary causal trend but is not yet the energy-normalized control.\n\n');
for i=1:height(G)
fprintf(fid,'- %s: exact first full-rank width %.6g mm; frozen maximum rank %.0f; Spearman(width, information) %.6g.\n', ...
    G.case_id(i),G.exact_first_full_rank_width_mm(i),G.frozen_max_rank(i),G.exact_information_monotonic_spearman(i));
end
fclose(fid);
end
function r=row0()
r=struct('case_id',"",'mode',"",'window_half_width_mm',NaN, ...
    'window_full_width_mm',NaN,'sample_count',NaN,'weight_energy',NaN, ...
    'rank_conditioned',NaN,'lambda_min_scaled',NaN,'condition_scaled',NaN);
end
