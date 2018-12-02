function out = thesis_2D_gaussian_toy_model(theta,y,z)
%JH 28/11/18
%toy model for thesis
% don't use inputs y and z, included so function has same arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out = [1,2].*theta + [2,1].*randn(1,2);
