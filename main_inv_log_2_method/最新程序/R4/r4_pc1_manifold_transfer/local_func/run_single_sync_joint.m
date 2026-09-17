function R=run_single_sync_joint(highMap,templateLib,cfg,staticState,eo)
%RUN_SINGLE_SYNC_JOINT  Single synchronous component with fixed EO.
cfg.singleSyncFrequencyHz=eo*(cfg.RPM_high/60);
R=run_single_frequency_joint(highMap,templateLib,cfg,staticState,'sync');
R.eo_id=eo; R.frequency_constraint="integer_EO";
end
