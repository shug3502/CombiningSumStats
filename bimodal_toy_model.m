function out = bimodal_toy_model(theta,y,z)
%JH 29/11/18
%bimodal toy model
% don't use inputs y and z, included so function has same arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out = sin(theta) + 0.1*randn(1,2);
