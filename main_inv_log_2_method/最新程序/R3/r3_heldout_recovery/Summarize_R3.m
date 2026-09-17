function [Summary,Paired] = Summarize_R3(Detail,P)
%SUMMARIZE_R3 Aggregate method probabilities and paired rescue/penalty events.

root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'local_func'));

gList=unique(Detail.g_truth_mm).';aList=unique(Detail.A_true_mm).';
names=["fixed","adaptive","state_matched"];
rows=repmat(summary_row(),numel(gList)*numel(aList)*3,1);ir=0;
for g=gList
    for A=aList
        for name=names
            q=abs(Detail.g_truth_mm-g)<1e-12&Detail.A_true_mm==A&Detail.method==name;
            if ~any(q),continue;end
            ir=ir+1;n=nnz(q);pv=NaN;lo=NaN;hi=NaN;pf=NaN;pj=NaN;
            flo=NaN;fhi=NaN;jlo=NaN;jhi=NaN;
            if A>0
                k=nnz(Detail.vib_success(q));pv=k/n;[lo,hi]=r3_wilson_interval(k,n,P.wilsonAlpha);
                kj=nnz(Detail.joint_success(q));pj=kj/n;
                [jlo,jhi]=r3_wilson_interval(kj,n,P.wilsonAlpha);
            else
                validFP=~ismissing(Detail.false_positive(q));nf=nnz(validFP);
                if nf>0
                    kf=nnz(Detail.false_positive(q));pf=kf/nf;
                    [flo,fhi]=r3_wilson_interval(kf,nf,P.wilsonAlpha);
                end
            end
            rows(ir)=struct('method',name,'g_truth_mm',g,...
                'delta_gap_mm',g-P.gReferenceMm,'A_true_mm',A,'n',n,...
                'P_vib',pv,'P_vib_wilson_low',lo,'P_vib_wilson_high',hi,...
                'P_FP',pf,'P_FP_wilson_low',flo,'P_FP_wilson_high',fhi,...
                'P_joint',pj,'P_joint_wilson_low',jlo,'P_joint_wilson_high',jhi,...
                'median_amplitude_error_mm',...
                median(Detail.amplitude_error_mm(q),'omitnan'),'median_gap_error_mm',...
                median(Detail.gap_error_mm(q),'omitnan'));
        end
    end
end
Summary=struct2table(rows(1:ir));

prows=repmat(paired_row(),numel(gList)*numel(aList),1);ip=0;
for g=gList
    for A=aList(aList>0)
        q=abs(Detail.g_truth_mm-g)<1e-12&Detail.A_true_mm==A;
        D=Detail(q,:);F=D(D.method=="fixed",{'high_seed','vib_success'});
        Adata=D(D.method=="adaptive",{'high_seed','vib_success'});
        F.Properties.VariableNames{2}='fixed_success';Adata.Properties.VariableNames{2}='adaptive_success';
        J=innerjoin(F,Adata,'Keys','high_seed');if isempty(J),continue;end
        ip=ip+1;prows(ip)=struct('g_truth_mm',g,'delta_gap_mm',g-P.gReferenceMm,...
            'A_true_mm',A,'n',height(J),'P_rescue',mean(J.adaptive_success&~J.fixed_success),...
            'P_penalty',mean(~J.adaptive_success&J.fixed_success),...
            'paired_gain',mean(J.adaptive_success)-mean(J.fixed_success));
    end
end
Paired=struct2table(prows(1:ip));
end

function r=summary_row(),r=struct('method',"",'g_truth_mm',NaN,'delta_gap_mm',NaN,...
    'A_true_mm',NaN,'n',NaN,'P_vib',NaN,'P_vib_wilson_low',NaN,...
    'P_vib_wilson_high',NaN,'P_FP',NaN,'P_FP_wilson_low',NaN,...
    'P_FP_wilson_high',NaN,'P_joint',NaN,'P_joint_wilson_low',NaN,...
    'P_joint_wilson_high',NaN,...
    'median_amplitude_error_mm',NaN,'median_gap_error_mm',NaN);end
function r=paired_row(),r=struct('g_truth_mm',NaN,'delta_gap_mm',NaN,...
    'A_true_mm',NaN,'n',NaN,'P_rescue',NaN,'P_penalty',NaN,'paired_gain',NaN);end
