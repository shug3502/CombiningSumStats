function out = toy_model(theta,y,z)
%JH 16/5/17
%simple toy model, uniform on [0,a] where a is unknown
% don't use inputs y and z, included so function has same arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = theta*rand(1,10);
out = sort(x);
