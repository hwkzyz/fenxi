function Restore_Config_Env_20241106(state)
%RESTORE_CONFIG_ENV_20241106 Restore environment values changed by a run.
for i = 1:numel(state.names)
    setenv(state.names{i}, state.values{i});
end
end
