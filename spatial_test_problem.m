function ss = spatial_test_problem(theta, t_end, dt, model)
%% JH 2/8/16
%% spatial test problem for ABC dist weights
%%%%%%%%%%%%%%%%%%%

if nargin<1
    theta = 10^-1;
    t_end=10;
    dt = 1;
end

nbox = 8;
num_reactions = 2*nbox-2;

if nargin<4
%create model
model = sbiomodel('Spatial model');

formulas = cell(num_reactions,1); %store the reaction formulae
for j=2:nbox-1
    formulas{j} = sprintf('s%d -> s%d',j,j+1);
    formulas{nbox+j-1} = sprintf('s%d -> s%d',j,j-1);
end
formulas{1} = sprintf('s%d -> s%d',1,2);
formulas{nbox} = sprintf('s%d -> s%d',nbox,nbox-1); %edge cases

%enter reactions
% and set reactions to be mass action
%r = cell(num_reactions,1);
kl = cell(num_reactions,1);
for i=1:num_reactions
r = addreaction(model, formulas{i});
kl{i} = addkineticlaw(r, 'MassAction');
end
%add rate constants for each reaction
%set kinetic law constants for each kinetic law
p=cell(num_reactions,1);
for i=1:num_reactions
%p{i}  = addparameter(kl{i}, sprintf('c%d',i),  'Value', theta); %1
p{i}  = addparameter(model, sprintf('c%d',i), theta); 
kl{i}.ParameterVariableNames = {sprintf('c%d',i)};
end
% %specify initial conditions
for i=1:nbox/2
model.species(i).InitialAmount = 10;    % s1
end
% %display model reaction objects
% model.Reactions
% %ditto species
% model.Species
end

%set parameters
for i=1:num_reactions
set(model.Parameters(i),'Value',theta)
end
%get the active configuration set for the model
cs = getconfigset(model,'active');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now compare to using tau leaping
cs.StopTime = t_end;
cs.SolverType = 'ssa';
solver = cs.SolverOptions;
%solver.LogDecimation = 10;
[t_etl, x_etl] = sbiosimulate(model);
if max(t_etl)<t_end
    t_etl = [t_etl;t_end];
    x_etl = [x_etl;x_etl(end,:)];
end
ts = timeseries(x_etl,t_etl,'Name','spatial');
ts1 = resample(ts,[dt:dt:t_end],'zoh');% resample(ts,[0:dt:t_end],'zoh');
ss = reshape(ts1.Data,1,[]);
%data is in form of tend/dt data points about each box


