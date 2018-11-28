function out = gaussian_toy_model(theta,y,z)
%JH 27/11/18
% gaussian toy problem as in prangle 2017
% don't use inputs y and z, included so function has same arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out = [theta + 0.1*randn(1), randn(1)];
