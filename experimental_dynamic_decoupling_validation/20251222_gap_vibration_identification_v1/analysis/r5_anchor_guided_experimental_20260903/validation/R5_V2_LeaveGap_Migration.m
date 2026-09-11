function T = R5_V2_LeaveGap_Migration(S, outputDir)
%R5_V2_LEAVEGAP_MIGRATION Global held-state interpolation benchmark.
T = R5_Evaluate_StaticMigration(S,'leave_gap',outputDir);
end
