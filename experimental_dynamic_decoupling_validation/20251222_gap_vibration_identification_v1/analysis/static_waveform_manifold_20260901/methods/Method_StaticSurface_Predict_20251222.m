function yhat = Method_StaticSurface_Predict_20251222(model, gap, varargin)
yhat=model.coeff(:,1)+model.coeff(:,2)./gap+model.coeff(:,3).*log(gap./model.g0);
end
