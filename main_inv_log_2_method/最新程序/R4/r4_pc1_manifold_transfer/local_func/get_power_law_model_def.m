function modelDef = get_power_law_model_def()
%get_power_law_model_def  Return the comparison power-law response model.

modelDef = table("power_law_p025_order2", "power_law", "", 0.25, 2, ...
    'VariableNames', {'model_label','response_model','basis_name','parameter','basis_order'});
end
