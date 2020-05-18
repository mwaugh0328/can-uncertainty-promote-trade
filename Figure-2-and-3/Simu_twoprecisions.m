%{
Can Global Uncertainty Promote International Trade?
Authors: Isaac Baley, Laura Veldkamp, and Mike Waugh 
Data: May 2020

---------------------------------------------------------------------------
---------------------------------------------------------------------------
                  SIMULATION (behind Figures 2 and 3)
---------------------------------------------------------------------------
---------------------------------------------------------------------------

Simulates the model for the cases with perfect and imperfect information.  
%}


close all;
clear;
clc;

% Random number generator
rng(1)


%% Parameters

% % Structural parameters
param.theta = 0.75;
param.sigma = 1 - param.theta;
param.tau   = 0;            % Iceberg cost
param.m_x   = 0;            % mean aggregate productivity country x;
param.m_y   = 0;            % mean aggregate productivity country y;
param.s_x   = 1;            % st.d. aggregate productivity country x;
param.s_y   = 1;            % st.d. aggregate productivity country y;
param.sig_x = sqrt(2);      % st.d. idiosincratic productivity country x;
param.sig_y = sqrt(2);      % st.d. idiosincratic productivity country y;

% % Grid parameters
% Country x
N = 11;
param.N_mu_x        = N;                                                   % grid size for domestic productivity
param.N_post_mu_y   = N;                                                   % grid size for posterior mean of foreign productivity
param.mu_x_min      = param.m_x - 3*param.s_x; 
param.mu_x_max      = param.m_x + 3*param.s_x;  
param.post_mu_y_min = param.m_y - 3*param.s_y;
param.post_mu_y_max = param.m_y + 3*param.s_y;

% Country y
N = 11;
param.N_mu_y        = N;                                                   % grid size for foreign productivity
param.N_post_mu_x   = N;                                                   % grid size for posterior mean of domestic productivity
param.mu_y_min      = param.m_y - 3*param.s_y;  
param.mu_y_max      = param.m_y + 3*param.s_y;  
param.post_mu_x_min = param.m_x - 3*param.s_x;
param.post_mu_x_max = param.m_x + 3*param.s_x;

% % Algorithm parameters
param.relax  = 1 - 0.05;                                                   % relaxation parameter for algorithm 
param.tol    = 10^(-5);                                                    % tolerance for iterations
param.N_quad = 10;                                                         % quadrature nodes

% Sets simulation size
param.T_sim  = 20;


%% Generates grids
% Linear splines at uniform nodes
fspace_x = fundef({'spli', nodeunif(param.N_mu_x,param.mu_x_min,param.mu_x_max),0,1},...     % Functional space for country x
    {'spli', nodeunif(param.N_post_mu_y,param.post_mu_y_min,param.post_mu_y_max),0,1});
fspace_y = fundef({'spli', nodeunif(param.N_mu_y,param.mu_y_min,param.mu_y_max),0,1},...     % Functional space for country y
    {'spli', nodeunif(param.N_post_mu_x,param.post_mu_x_min,param.post_mu_x_max),0,1});

state_x        = gridmake(funnode(fspace_x));  
grid.mu_x      = state_x(1:param.N_mu_x,1);
grid.post_mu_y = state_x(1:param.N_mu_x:length(state_x),2);
state_y        = gridmake(funnode(fspace_y));         
grid.mu_y      = state_y(1:param.N_mu_y,1);
grid.post_mu_x = state_y(1:param.N_mu_y:length(state_y),2);
state          = gridmake(funnode(fspace_x),funnode(fspace_y));

%% Perfect Info Guess
[q_x,q_y,Psi_PI,Gamma_PI] = PI(param,grid);

q_x      = reshape(q_x,[param.N_mu_x*param.N_mu_y 1]);
Psi_PI   = reshape(Psi_PI,[param.N_mu_x*param.N_mu_y 1]);
Gamma_PI = reshape(Gamma_PI,[param.N_mu_x*param.N_mu_y 1]);
coef     = guess(param,Psi_PI,q_x,Gamma_PI);

%% Imperfect Info

% Draws endowment shocks (same in both information scenarios)
mu_simu      = mvnrnd([param.m_x param.m_y], ...
                [param.s_x^2 0; 0 param.s_y^2],param.T_sim);
mu_simu(:,1) = max(min(mu_simu(:,1),param.mu_x_max),param.mu_y_min);
mu_simu(:,2) = max(min(mu_simu(:,2),param.mu_y_max),param.mu_y_min);

mu_weights = normpdf(grid.mu_x,param.m_x,param.s_x);
mu_weights = mu_weights/sum(mu_weights);

% Sets two levels of singal precision
N_eta   = 2;
eta_min = 0.0000000001; %perfect info
eta_max = 100000000000; %imperfect info
eta_vec = linspace(eta_min,eta_max,N_eta);

% Inizialization
Ep_eta_r             = zeros(param.N_mu_x,param.N_mu_y,N_eta);
Vp_eta_r             = zeros(param.N_mu_x,param.N_mu_y,N_eta);
Cp_eta_r             = zeros(param.N_mu_x,param.N_mu_y,N_eta);
Exports_eta_r        = zeros(param.N_mu_x,param.N_mu_y,N_eta);
meanEp_eta           = zeros(1,N_eta);
meanVp_eta           = zeros(1,N_eta);
meanCp_eta           = zeros(1,N_eta);
meanExports_eta      = zeros(1,N_eta);
meanxEp_eta          = zeros(param.N_mu_y,N_eta);
meanxVp_eta          = zeros(param.N_mu_y,N_eta);
meanxCp_eta          = zeros(param.N_mu_y,N_eta);
meanxExports_eta     = zeros(param.N_mu_y,N_eta);
meanyEp_eta          = zeros(param.N_mu_x,N_eta);
meanyVp_eta          = zeros(param.N_mu_x,N_eta);
meanyCp_eta          = zeros(param.N_mu_x,N_eta);
meanyExports_eta     = zeros(param.N_mu_x,N_eta);
post_weights_eta     = zeros(param.N_mu_x,N_eta);

%One loop for each level of signal precision
for j = 1:N_eta 
    
    disp('Iteration')
    disp(j)
    
    %Set signal noises, uncertainties
    param.p_eta_y    = eta_vec(j);                                         % Foreign's perception of domestic noise
    param.eta_y      = eta_vec(j);                                         % Domestic signal noise
    param.p_eta_x    = eta_vec(j);                                         % Domestic's perception of foreign noise 
    param.eta_x      = eta_vec(j);                                         % Foreign signal noise
    param.post_s_y   = 1/sqrt(1/param.s_y^2 + 1/param.eta_y^2);            % Posterior uncertainty about foreign productivity
    param.post_s_x   = 1/sqrt(1/param.s_x^2 + 1/param.eta_x^2);            % Posterior uncertainty about domestic productivity
    param.second_s_x = 1/(param.p_eta_x/param.s_x^2 + 1/param.p_eta_x);    % Second order uncertainty about domestic productivity
    param.second_s_y = 1/(param.p_eta_y/param.s_y^2 + 1/param.p_eta_y);    % Second order uncertainty about foreign productivity
    
    %Draw signals using corresponding precisions
    signal_simu = mvnrnd(mu_simu,[param.eta_x^2 0;0 param.eta_y^2]);       % Generate signals centered at true realization
    post_var    = [1/(param.s_x^(-2)/param.eta_x^(-2) + 1),...             % Posterior variances
                   1/(param.s_y^(-2)/param.eta_y^(-2) + 1)];
    post_simu   = [(1-post_var(1))*param.m_x+post_var(1)*signal_simu(:,1),... % Second state, posterior mean of foreign productivity
                   (1-post_var(2))*param.m_y+post_var(2)*signal_simu(:,2)];
    post_simu(:,1) = max(min(post_simu(:,1),param.mu_x_max),param.mu_x_min);
    post_simu(:,2) = max(min(post_simu(:,2),param.mu_y_max),param.mu_y_min);
    
    %Calls modelsolve
    [coef,Ep, Vp, Cp, Exports] = modelsolve(coef,param,state_x,state_y,fspace_x,fspace_y);   
    
    post_weights_eta(:,j) = normpdf(grid.mu_y,param.m_y,param.s_y);
    post_weights_eta(:,j) = post_weights_eta(:,j)/sum(post_weights_eta(:,j));
    
    % Expected terms of trade q
    Ep_eta_r(:,:,j)   = reshape(Ep,param.N_mu_x,param.N_post_mu_y);
    meanyEp_eta(:,j)  = Ep_eta_r(:,:,j)*post_weights_eta(:,j);
    meanxEp_eta(:,j)  = mu_weights'*Ep_eta_r(:,:,j);
    meanEp_eta(:,j)   = mu_weights'*Ep_eta_r(:,:,j)*post_weights_eta(:,j);
    
    % Variance of terms of trade q
    Vp_eta_r(:,:,j)   = reshape(Vp,param.N_mu_x,param.N_post_mu_y);
    meanyVp_eta(:,j)  = Vp_eta_r(:,:,j)*post_weights_eta(:,j);
    meanxVp_eta(:,j)  = mu_weights'*Vp_eta_r(:,:,j);
    meanVp_eta(:,j)   = mu_weights'*Vp_eta_r(:,:,j)*post_weights_eta(:,j);
    
    % Coeff of variation of trade q
    Cp_eta_r(:,:,j)   = reshape(Cp,param.N_mu_x,param.N_post_mu_y);
    meanyCp_eta(:,j)  = Cp_eta_r(:,:,j)*post_weights_eta(:,j);
    meanxCp_eta(:,j)  = mu_weights'*Cp_eta_r(:,:,j);
    meanCp_eta(:,j)   = mu_weights'*Cp_eta_r(:,:,j)*post_weights_eta(:,j);
    
    % Total exports
    Exports_eta_r(:,:,j)   = reshape(Exports,param.N_mu_x,param.N_post_mu_y);
    meanyExports_eta(:,j)  = Exports_eta_r(:,:,j)*post_weights_eta(:,j);
    meanxExports_eta(:,j)  = mu_weights'*Exports_eta_r(:,:,j);
    meanExports_eta(:,j)   = mu_weights'*Exports_eta_r(:,:,j)*post_weights_eta(:,j); 
    
    [T_x_simu(:,j),T_y_simu(:,j),q_simu(:,j)]=simulation(param,fspace_x,fspace_y,coef,mu_simu,post_simu);
end

save results_simu

