function report = R5_V3_LeaveState(profileCells, lambda)
%R5_V3_LEAVESTATE Test slow-prior path against each independently selected window.
% This is a temporal consistency test, not a physical-gap validation.
arguments, profileCells cell; lambda (1,1) double {mustBeNonnegative} = 0; end
path=R5_Select_StatePath_SlowPrior(profileCells,lambda);
independent=nan(1,numel(profileCells));
for k=1:numel(profileCells), independent(k)=profileCells{k}.selected.stateValue; end
report=struct('path',path,'independentState',independent, ...
    'stateChanged',independent~=path.stateValue,'status','temporal_effective_state_diagnostic');
end
