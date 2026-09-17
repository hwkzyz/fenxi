function modelDefs = get_inv_log_2_model_suite()
%get_inv_log_2_model_suite  Models used in closed-loop comparisons.

modelDefs = [get_baseline_model_def(); get_inv_log_2_model_def(); get_power_law_model_def()];
end
