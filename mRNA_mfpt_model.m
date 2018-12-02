function out = mRNA_mfpt_model(theta,y,z)
%JH 30/11/18
% 3D velocity jump process for mRNA transport
% don't use inputs y and z, included so function has same arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpt = summary_statistic_calculator_fpt(theta,5,1);
out = median(fpt);
