function [q_x,q_y,Psi_PI,Gamma_PI] = PI(param,grid)

f_x       = exp(grid.mu_x + 0.5*param.sig_x^2);
f_y       = exp(grid.mu_y + 0.5*param.sig_y^2) ; 
f         = gridmake(f_x,f_y);                                                     
f         = f(:,1)./f(:,2); 
Psii_PI   = @(q) 1./(1 + ...
    ((1 + param.tau).*q).^(param.theta/(1-param.theta)));
Gammaa_PI = @(q) 1./(1 + ...
    ((1 + param.tau)./q).^(param.theta/(1-param.theta)));     
q         = fsolve(@(x) x-f.*Psii_PI(x)./Gammaa_PI(x), ...
    ones(param.N_mu_x*param.N_mu_y,1));
q_x       = reshape(q, [param.N_mu_x param.N_mu_y]);
q_y       = 1./q_x;
Psi_PI    = Psii_PI(q_x);
Gamma_PI  = Gammaa_PI(q_x);