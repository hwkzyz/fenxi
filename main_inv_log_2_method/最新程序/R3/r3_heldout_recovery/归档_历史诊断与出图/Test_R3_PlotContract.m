function Test_R3_PlotContract()
%TEST_R3_PLOTCONTRACT Render both figure scripts from a complete mock table.

P=R3_Protocol();g=sort(P.gMainMm);A=P.amplitudeMainMm;
names=["fixed","adaptive","state_matched"];
method=strings(0,1);gt=[];at=[];pv=[];pf=[];pj=[];
for m=names
    for x=[0 A]
        for z=g
            method(end+1,1)=m; %#ok<AGROW>
            gt(end+1,1)=z;at(end+1,1)=x; %#ok<AGROW>
            pv(end+1,1)=double(x>0)*0.95; %#ok<AGROW>
            pf(end+1,1)=double(x==0)*0.05; %#ok<AGROW>
            pj(end+1,1)=double(x>0)*0.92; %#ok<AGROW>
        end
    end
end
Summary=table(method,gt,at,pv,pf,pj,'VariableNames',...
    {'method','g_truth_mm','A_true_mm','P_vib','P_FP','P_joint'});
[G,AA]=meshgrid(g,A);
Paired=table(G(:),AA(:),0.4*ones(numel(G),1),...
    'VariableNames',{'g_truth_mm','A_true_mm','P_rescue'});
Result=struct('P',P,'Summary',Summary,'Paired',Paired);
mockFile=fullfile(tempdir,'r3_plot_contract.mat');save(mockFile,'Result');
F1=Plot_R3_Fig3(mockFile);F2=Plot_R3_Supplementary(mockFile);
assert(all(structfun(@isfile,F1))&&all(structfun(@isfile,F2)));
cellfun(@delete,struct2cell(F1));cellfun(@delete,struct2cell(F2));
delete(mockFile);close all;disp('PLOT_CONTRACT_OK');
end
