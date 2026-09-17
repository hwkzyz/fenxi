function Result=Run_StructuredLayoutCoherenceAudit(nRandom)
%RUN_STRUCTUREDLAYOUTCOHERENCEAUDIT Rank layouts by worst EO-subspace coherence.
if nargin<1,nRandom=2000;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;eoGrid=6:30;
x=linspace(lib.domain(1),lib.domain(2),121).';fx=eval_gap_derivative(lib,.5,x);
layouts=base_layouts(cfg.alpha_k);rows=repmat(row0(),numel(layouts),1);
for i=1:numel(layouts)
    s=layout_score(layouts(i).alpha,eoGrid,x,fx,cfg.R_tip);
    rows(i)=pack(layouts(i).name,layouts(i).alpha,s);
end
rng(20260807,'twister');bestScore=Inf;bestAlpha=[];
for i=1:nRandom
    a=sort([0,2*pi*rand(1,4)]);
    sep=diff([a,a(1)+2*pi]);if min(sep)<.20,continue;end
    s=layout_score(a,eoGrid,x,fx,cfg.R_tip);
    if s.maxCoherence<bestScore,bestScore=s.maxCoherence;bestAlpha=a;bestDetail=s;end
end
rows(end+1)=pack("optimized5",bestAlpha,bestDetail);
Audit=sortrows(struct2table(rows),'max_subspace_coherence','ascend');
out=fullfile(root,'output','structured_layout_coherence');if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'layout_coherence.csv'));
save(fullfile(out,'layout_coherence.mat'),'Audit','eoGrid','x','nRandom');disp(Audit);
Result=struct('Audit',Audit,'bestAlpha',bestAlpha,'outputDir',out);
end

function L=base_layouts(current)
gold=(1+sqrt(5))/2;
L=[layout("current3",current),layout("uniform3",(0:2)*2*pi/3),...
 layout("golden3",sort(mod((0:2)*2*pi/gold,2*pi))),...
 layout("uniform4",(0:3)*2*pi/4),layout("uniform5",(0:4)*2*pi/5),...
 layout("golden5",sort(mod((0:4)*2*pi/gold,2*pi)))];
end

function q=layout(name,alpha),q=struct('name',string(name),'alpha',alpha);end

function s=layout_score(alpha,eoGrid,x,fx,R)
theta=reshape(x/R+alpha,[],1);w=repmat(abs(fx),numel(alpha),1);
Q=cell(numel(eoGrid),1);
for i=1:numel(eoGrid)
    X=w.*[sin(eoGrid(i)*theta),cos(eoGrid(i)*theta)];[Q{i},~]=qr(X,0);
end
best=-Inf;pair=[NaN NaN];
for i=1:numel(eoGrid)-1
    for j=i+1:numel(eoGrid)
        rho=norm(Q{i}'*Q{j},2);
        if rho>best,best=rho;pair=[eoGrid(i),eoGrid(j)];end
    end
end
s=struct('maxCoherence',best,'worstPair',pair,'minSeparation',min(diff([alpha,alpha(1)+2*pi])));
end

function r=pack(name,alpha,s)
r=struct('layout',name,'sensor_count',numel(alpha),'alpha_rad',string(mat2str(alpha,6)),...
    'max_subspace_coherence',s.maxCoherence,'worst_eo_pair',string(mat2str(s.worstPair)),...
    'minimum_angular_separation',s.minSeparation);
end

function r=row0()
r=struct('layout',"",'sensor_count',NaN,'alpha_rad',"",...
    'max_subspace_coherence',NaN,'worst_eo_pair',"",'minimum_angular_separation',NaN);
end
