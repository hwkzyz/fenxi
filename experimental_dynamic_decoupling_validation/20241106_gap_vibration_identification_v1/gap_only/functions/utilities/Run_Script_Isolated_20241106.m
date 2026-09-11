function Run_Script_Isolated_20241106(scriptFile)
%RUN_SCRIPT_ISOLATED_20241106 Run a legacy script in its own function workspace.
assert(isfile(scriptFile), 'Missing package script: %s', scriptFile);
run(scriptFile);
end
