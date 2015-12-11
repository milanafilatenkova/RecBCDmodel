
%%  MC_RECBCD

% performs MCMC algorithm
% in 3D-parameter space of the model
% taking as an input 1) the initial point (input: theta_0)
% ... 2) the total number of steps (input: Number_iterations)
% ... 3) the maximum step size (input: epsilon)
% the output is the distribution of the parameters (output: Theta_PDF)
function [Theta, acceptance] = MC_RecBCDsimple(data, ChiPositions, Number_iterations,theta_0,epsilon)

% initialise PDF on the parameter space


Theta = zeros(3,Number_iterations);

theta = zeros(3,1);
theta_new = zeros(3,1);


% initialise maximum step size
epsilon_s = epsilon(1);
epsilon_chi = epsilon(2);
epsilon_m = epsilon(3);


% initialise the count of positive events
% when the proposal is accepted
N_positive = 0;

% initialise the start point of MCMC algorithm
theta = theta_0;

for i=1:Number_iterations
    
    %1) proposing innovation (for example a new x' can be drawn from uniform
    %distribution having its center in x)
    innov_s = unifrnd(-epsilon_s,epsilon_s,1);
    innov_chi = unifrnd(-epsilon_chi,epsilon_chi,1);
    innov_m = unifrnd(-epsilon_m,epsilon_m,1);

    
    innov = [innov_s,innov_chi,innov_m];
 
    theta_new = theta + innov;
    
    %2) Compute the ratio of probabilities p(theta_new)/p(theta)
    delta_log = RecBCDdeterministic([theta_new(1),...
            theta_new(2),theta_new(3)],data,ChiPositions)...
               - RecBCDdeterministic([theta(1),...
                   theta(2),theta(3)],data,ChiPositions);
    
    
    % Check if the probability ratio is bigger then a random number. If so
    %make a move into theta_new
    prob = min(1,exp(delta_log));
    
    if rand(1) < prob & theta_new(1) > 0 ...
            & theta_new(2) > 0 ...
            & theta_new(2) < 1 & theta_new(3) > 0 ...
            & theta_new(3) < 1  ...
            & theta_new(1) < 1 
       
        theta = theta_new;
        N_positive = N_positive +1;
   
        
    end 
    Theta(:,i) = theta;
    acceptance(i) =  N_positive/i;
    
   
end
%%




