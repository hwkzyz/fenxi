function modelDef = get_baseline_model_def()
%get_baseline_model_def  Return the baseline 1/g response model.

modelDef = table("baseline_1_over_g", "inv_g_linear", "", NaN, 1, ...
    'VariableNames', {'model_label','response_model','basis_name','parameter','basis_order'});
end
