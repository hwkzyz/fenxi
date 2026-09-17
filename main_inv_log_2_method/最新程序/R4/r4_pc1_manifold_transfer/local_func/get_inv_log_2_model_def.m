function modelDef = get_inv_log_2_model_def()
%get_inv_log_2_model_def  Return the current inv_log_2 response model.

modelDef = table("field_basis_inv_log_2", "field_basis", "inv_log_2", NaN, 3, ...
    'VariableNames', {'model_label','response_model','basis_name','parameter','basis_order'});
end
