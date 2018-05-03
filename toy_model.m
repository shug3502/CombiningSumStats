function out = toy_model(theta,y,z)
%JH 16/5/17
%simple toy model, uniform on [0,a] where a is unknown
% don't use inputs y and z
x = theta*rand(1,100);
%out = [randn(1,100),sort(x)];
out = [max(x),mean(x),7,prod(x/mean(x)),sum(x.^2),prod(log(x))];
