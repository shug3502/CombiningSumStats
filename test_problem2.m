function ss = test_problem2(theta,t_end,recording_interval)
A = 10; %initial population level
if nargin<1
    phi = 0.1; %decay rate
    sigma = 0.01;
    t_end=100;
    recording_interval = 1;
else
    phi = theta(1);
    sigma = theta(2);
end
noise = 0; %magnitude of observation noise

time = 0;

%set up recording
%t_end = 100;
%recording_interval = 1; %now required as an input
n_t_pts = round(t_end/recording_interval);
A_rec = zeros(n_t_pts,1);
j=0;  %j=0;

while (A>0 && time<t_end)
  tau = -1/(phi*A)*log(rand(1));
  time = time+tau;
  while (time>recording_interval*j)
    j=j+1;
    A_rec(j)=A;
    %tau = tau-recording_interval;
  end
  %perform reactions
  A = A-1;
end
A_data  = A_rec + sqrt(noise)*randn(max(n_t_pts,numel(A_rec)),1);

%z = (A_data(1:t_end) - A_data(2:(1+t_end)))./A_data(2:(1+t_end))
%z1 = find(A_data<=0,1)-1; %the first time that we see extinction of the population
%if isempty(z1)
%    %warning('no such index found');
%    z1 = t_end;
%end
ss = A_data(1:n_t_pts)';
%ss(1:2:t_end+1) = randn(1,t_end/2+1); %don't currently corrupt data
ss = [ss, sigma*randn(1)];
%%z2 = sum(A_data);
%%ss = [z1, sigma*randn(1),z2];
%%ss = [z1, -z1, rand(1)];
%%ss = [ss, exp(ss/50), tanh(ss)];
