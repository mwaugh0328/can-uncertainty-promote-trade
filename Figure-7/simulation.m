%% Simulation, returns trade and utility correlation. 

function [corr_T,corr_C]=simulation(param,fspace_x,fspace_y,coef,mu_simu,post_simu)

c_g1 = coef(:,1);       
c_g2 = coef(:,2);  
c_h1 = coef(:,3);  
c_h2 = coef(:,4);  

% Evaluate policies at simulated states
g1_simu    = funeval(c_g1,fspace_x,[mu_simu(:,1) post_simu(:,2)]);
g2_simu    = funeval(c_g2,fspace_x,[mu_simu(:,1) post_simu(:,2)]);
Psi_1_simu = 1./(1+(1+param.tau)^(param.theta/(1-param.theta)).*(g1_simu./g2_simu).^(1/(1-param.theta)));
h1_simu    = funeval(c_h1,fspace_y,[mu_simu(:,2) post_simu(:,1)]);
h2_simu    = funeval(c_h2,fspace_y,[mu_simu(:,2) post_simu(:,1)]);
Gam_1_simu = 1./(1+(1+param.tau)^(param.theta/(1-param.theta)).*(h1_simu./h2_simu).^(1/(1-param.theta)));

% Fundamentals
f_x_simu = exp(mu_simu(:,1)+0.5*param.sig_x^2);
f_y_simu = exp(mu_simu(:,2)+0.5*param.sig_y^2);
f_simu   = f_x_simu./f_y_simu;    

% Aggregate variables
q_simu    = f_simu.*Psi_1_simu./Gam_1_simu;

% X country
T_x_simu  = f_x_simu.*    Psi_1_simu  ;                                                                                  
C_x_simu  = f_x_simu.* (1-Psi_1_simu) ;                                               
C_y_simu  = T_x_simu./ ((1+param.tau).*q_simu); 
C_simu    = (C_x_simu.^param.theta + C_y_simu.^param.theta).^(1/param.theta);   

% Y country
T_y_simu  = f_y_simu.*  Gam_1_simu ; 
Cf_y_simu = f_y_simu.* (1-Gam_1_simu);
Cf_x_simu = T_y_simu.* q_simu/(1+param.tau);
Cf_simu = (Cf_x_simu.^param.theta + Cf_y_simu.^param.theta).^(1/param.theta);  

corr_T = corr(T_x_simu,T_y_simu);
corr_C = corr(C_simu,Cf_simu);