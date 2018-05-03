function d_hat = hellinger_knn_estimator(x_samples,y_samples,k)
%based on estimator in poczos et al 2012 and sutherland et al 2012
%compare to python code by dougal sutherland in skl-groups
%JH 19/12/17
%%%%%%%%%%%%%%%%%%%%%%%%
% x_samples has dim n_samples by n_sumstats
% y_samples has dim n_samples by s_sumstats

n_samples = size(x_samples,1);
d = size(x_samples,2); % i think this is what the dimensionality d refers to
alpha = -1/2; % method is for general alpha and beta, we set these for hellinger
beta = 1/2; % carfeul because these are also inbuilt function names
if size(y_samples,1)~= n_samples
    error('We assume same number of x and y samples, although estimator is valid in general');
end

[IDx,Dx] = knnsearch(x_samples,x_samples,'K',k);
[IDy,Dy] = knnsearch(y_samples,y_samples,'K',k);
rho_k = Dx(:,k); %distance of kth nearest neighbour in each case
nu_k = Dy(:,k); 

d_ab = sum(rho_k.^(-d*alpha).*nu_k.^(-d*beta));

ball_vol = (pi^(d/2)/gamma(d/2 + 1))^(-alpha-beta)*(gamma(k)^2)/(gamma(k-alpha)*gamma(k-beta)); 
correction_factor = ball_vol/(n_samples*((n_samples-1)^alpha)*(n_samples^beta)); 

d_hat = 1 - correction_factor*d_ab;
