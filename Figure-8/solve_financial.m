
%{
Can Global Uncertainty Promote International Trade?
Authors: Isaac Baley, Laura Veldkamp, and Mike Waugh 
Data: May 2020

---------------------------------------------------------------------------
---------------------------------------------------------------------------
        Solves policies with two types of agents: 

%     1)  A fraction alpha of agents operates udner complete markets (perfect information) 
%     2)  A fraction 1-alpha financial autarchy (imperfect info with precision zero, posterior are equal to the prior). 
% 
      (Behind Figure 8) 
---------------------------------------------------------------------------
---------------------------------------------------------------------------

%}

clear all
close all
clc

%% Parameters

% Vector of alpha's (fraction of agents with complete markets)
n_alpha = 5;
alpha_vec=linspace(0.1,0.9,n_alpha); 

%structural
param.theta = 0.3;             % Low elasticity of substitution
param.sigma = 1 - param.theta; %Baseline CES-like model  
param.tau   = 0;               % Iceberg cost
param.m_x   = 0;               % mean aggregate productivity country x;
param.m_y   = 0;               % mean aggregate productivity country y;
param.s_x   = 1;               % st.d. aggregate productivity country x;
param.s_y   = 1;               % st.d. aggregate productivity country y;
param.sig_x = sqrt(2);      % st.d. idiosincratic productivity country x;
param.sig_y = sqrt(2);      % st.d. idiosincratic productivity country y;

% Grid parameters
% Country x
N = 5;
param.N_mu_x        = N;	% grid size for domestic productivity
param.N_post_mu_y   = N;	% grid size for posterior mean of foreign productivity
param.mu_x_min      = param.m_x - 3*param.s_x; 
param.mu_x_max      = param.m_x + 3*param.s_x;  
param.post_mu_y_min = param.m_y - 3*param.s_y;
param.post_mu_y_max = param.m_y + 3*param.s_y;

% Country y
N = 5;
param.N_mu_y        = N;	% grid size for domestic productivity
param.N_post_mu_x   = N;    % grid size for posterior mean of foreign productivity
param.mu_y_min      = param.m_y - 3*param.s_y;  
param.mu_y_max      = param.m_y + 3*param.s_y;  
param.post_mu_x_min = param.m_x - 3*param.s_x;
param.post_mu_x_max = param.m_x + 3*param.s_x;

% % Algorithm parameters
param.relax  = 1 - 0.05;                                                   % relaxation parameter for algorithm 
param.tol    = 10^(-5);                                                    % tolerance for iterations
param.N_quad = 10;                                                         % quadrature points

% Simulation param
param.T_sim = 1000;                                                       % periods to simulate

%% Set up functional spaces

% Functional spaces
fspace_x = fundef({'spli', nodeunif(param.N_mu_x,param.mu_x_min,param.mu_x_max),0,1});
% Functional space for X country
fspace_y = fundef({'spli', nodeunif(param.N_mu_y,param.mu_y_min,param.mu_y_max),0,1});
% States (for policies)
state_x        = gridmake(funnode(fspace_x));  
state_y        = gridmake(funnode(fspace_y));   
state          = gridmake(funnode(fspace_x),funnode(fspace_y));     

%Initial coefficients
c_g1     = funfitxy(fspace_x,state_x,0.5*ones(param.N_mu_x,1));
c_g2     = c_g1;
c_h1     = c_g1;
c_h2     = c_g1;
coef_ini = [c_g1,c_g2,c_h1,c_h2];


%% Run and fix simulation for aggregate productivities...

mu_simu= mvnrnd([param.m_x param.m_y], ...
                [param.s_x^2 0; 0 param.s_y^2],param.T_sim);
mu_simu(:,1) = max(min(mu_simu(:,1),param.m_x+4*param.s_x),param.m_x-4*param.s_x);
mu_simu(:,2) = max(min(mu_simu(:,2),param.m_y+4*param.s_y),param.m_y-4*param.s_y);

for i=1:n_alpha  % Solves model for each level of alpha (fraction of agents with complete markets)
i
alpha = alpha_vec(i);

% Solve for policies and coefficients
coef_alpha = financial_contracts(param,fspace_x,fspace_y,coef_ini,alpha);

% Run simulation
[stats_fund(:,:,i), stats_dom(:,:,i), stats_for(:,:,i),    ...
  corr_fund(:,:,i), corr_dom(:,:,i) ,  corr_for(:,:,i) ]  ...
= simu_financial(param,fspace_x,fspace_y,coef_alpha,mu_simu,alpha);
end

save results_financial



