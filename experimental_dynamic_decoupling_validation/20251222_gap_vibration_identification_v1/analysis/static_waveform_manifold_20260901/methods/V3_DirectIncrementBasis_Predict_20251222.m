function deltaY = V3_DirectIncrementBasis_Predict_20251222(model, lowGap, targetGap, lowWave)
deltaY=model.coeff*[1;1/targetGap-1/lowGap;log(targetGap/lowGap)]; deltaY(~isfinite(lowWave))=NaN;
end
