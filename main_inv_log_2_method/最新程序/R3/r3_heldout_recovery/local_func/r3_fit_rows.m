function rows = r3_fit_rows(R3,C0,gTruth,Atrue,phi,sigmaV,highSeed,thresholds)
%R3_FIT_ROWS Generate one record and run all three comparator fits.

[~,highMap,truth]=r3_generate_high_case(R3,gTruth,Atrue,phi,sigmaV,highSeed);
Fits=r3_run_three_methods_global_final(R3,C0,highMap,gTruth);
names=["fixed","adaptive","state_matched"];
rows=repmat(r3_pack_result(Fits.fixed,truth,R3.P,NaN,C0),3,1);
for i=1:3
    name=names(i);threshold=NaN;
    if nargin>=8&&isstruct(thresholds)&&isfield(thresholds,char(name))
        threshold=thresholds.(char(name));
    end
    rows(i)=r3_pack_result(Fits.(char(name)),truth,R3.P,threshold,C0);
end
end
