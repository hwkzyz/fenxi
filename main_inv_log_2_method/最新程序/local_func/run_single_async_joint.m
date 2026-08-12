function R=run_single_async_joint(highMap,templateLib,cfg,staticState)
%RUN_SINGLE_ASYNC_JOINT  Single asynchronous component with continuous f.
R=run_single_frequency_joint(highMap,templateLib,cfg,staticState,'async');
R.frequency_constraint="continuous";
end
