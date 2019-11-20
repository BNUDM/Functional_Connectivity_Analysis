function [rCorr, pCorr] = FC_ROI( betaSeries_ROI1, betaSeries_ROI2 )
% [rCorr,pCorr]=FC_ROI( betaSeries_ROI1, betaSeries_ROI2 )
% Perform task-based functional connectivity analysis between ROIs.
% Input:
%   betaSeries_ROI1 - the betavalue series of ROI1, which is the output of DeconvGLM 
%   betaSeries_ROI2 - the betavalue series of ROI2, which is the output of DeconvGLM 
% Output:
%   rCorr - Pearson's Correlation Coefficient
%   pCorr - the P value

  roi1 = betaSeries_ROI1;
  roi2 = betaSeries_ROI2;
  [FC_r, FC_p] = corr(roi1', roi2');
  rCorr = FC_r;
  pCorr = FC_p;  
end
