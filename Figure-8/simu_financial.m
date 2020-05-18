function [stats_fund, stats_dom, stats_for,  ...
          corr_fund , corr_dom , corr_for   ...
          ]=simu_financial(param,fspace_x,fspace_y,coef,mu_simu,alpha)

c_g1 = coef(:,1);
c_g2 = coef(:,2);
c_h1 = coef(:,3);
c_h2 = coef(:,4);

% Evaluate policies at simulated states
g1_simu    = funeval(c_g1,fspace_x,mu_simu(:,1));
g2_simu    = funeval(c_g2,fspace_x,mu_simu(:,1));
Psi_1_simu = 1./(1+(1+param.tau)^(param.theta/(1-param.theta)).*(g1_simu./g2_simu).^(1/(1-param.theta)));
h1_simu    = funeval(c_h1,fspace_y,mu_simu(:,2));
h2_simu    = funeval(c_h2,fspace_y,mu_simu(:,2));
Gam_1_simu = 1./(1+(1+param.tau)^(param.theta/(1-param.theta)).*(h1_simu./h2_simu).^(1/(1-param.theta)));

% Fundamentals
f_x_simu = exp(mu_simu(:,1)+0.5*param.sig_x^2);
f_y_simu = exp(mu_simu(:,2)+0.5*param.sig_y^2);
f_simu   = f_x_simu./f_y_simu;    

% Aggregate variables
  options=optimset('fsolve');
  options=optimset(options,'display','off'); 

q_simu = fsolve(@(x) ...
                 x - f_simu.*...
                 (alpha./(1+x.^(param.theta/(1-param.theta)))+(1-alpha).*Psi_1_simu)./(alpha./(1+x.^(-param.theta/(1-param.theta)))+(1-alpha).*Gam_1_simu),ones(param.T_sim,1),options );

% X country
      T_x_simu  = (1-alpha).*f_x_simu.* Psi_1_simu  ;  
 Omega_x_simu   = 1./(1+q_simu.^(param.theta/(1-param.theta)));
tilde_T_x_simu  = alpha.*f_x_simu.*Omega_x_simu;
 T_x_total_simu = T_x_simu  + tilde_T_x_simu; 
       share_x  = alpha.*Omega_x_simu + (1-alpha).*Psi_1_simu;
       
% Y country
      T_y_simu  = (1-alpha).*f_y_simu.* Gam_1_simu ; 
 Omega_y_simu   = 1./(1+q_simu.^(-param.theta/(1-param.theta)));
tilde_T_y_simu  = alpha*f_y_simu.*Omega_y_simu ;
 T_y_total_simu = T_y_simu  + tilde_T_y_simu; 
        share_y = alpha.*Omega_y_simu + (1-alpha).*Gam_1_simu;

% Results: Fundamentals, Domestic, Foreign
results_fund   = [f_x_simu   f_y_simu    f_simu    q_simu];
results_dom    = [Psi_1_simu Omega_x_simu share_x  T_x_simu  tilde_T_x_simu  T_x_total_simu];  
results_for    = [Gam_1_simu Omega_y_simu share_y  T_y_simu  tilde_T_y_simu  T_y_total_simu]; 

% Statistics (Times series)
stats_fund  = [mean(results_fund); std(results_fund) ; std(results_fund) ./mean(results_fund) ]'; 
stats_dom   = [mean(results_dom);  std(results_dom)  ; std(results_dom)  ./mean(results_dom)  ]';
stats_for   = [mean(results_for);  std(results_for)  ; std(results_for)  ./mean(results_for)  ]'; 

% Correlation matrix
 corr_fund = corr(results_fund);
 corr_dom  = corr(results_dom);
 corr_for  = corr(results_for);