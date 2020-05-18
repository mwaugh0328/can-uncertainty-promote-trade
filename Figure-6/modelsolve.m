function [coef,Ep, Vp, Cp, Exports] = modelsolve(coef,param,state_x,state_y,fspace_x,fspace_y)      
c_g1_new   = coef(:,1);
c_g2_new   = coef(:,2);
c_h1_new   = coef(:,3);
c_h2_new   = coef(:,4);

exp_x_1    = zeros(length(state_x),1);
exp_x_2    = zeros(length(state_x),1);
exp_qx     = zeros(length(state_x),1);
var_qx     = zeros(length(state_x),1);
cv_qx      = zeros(length(state_x),1);
exp_px     = zeros(length(state_x),1);
var_px     = zeros(length(state_x),1);
cv_px      = zeros(length(state_x),1);

exp_y_1    = zeros(length(state_y),1);
exp_y_2    = zeros(length(state_y),1);
exp_qy     = zeros(length(state_y),1);
var_qy     = zeros(length(state_y),1);
cv_qy      = zeros(length(state_y),1);
exp_py     = zeros(length(state_y),1);
var_py     = zeros(length(state_y),1);
cv_py      = zeros(length(state_y),1);

        
% Fixed point algorithm
diff = 100;
i    = 0;
        
while(diff>param.tol && i<1000)
c_g1 = c_g1_new;
c_g2 = c_g2_new;
c_h1 = c_h1_new;
c_h2 = c_h2_new;
            
%%%%%%%%%%%%%%%%%%%%%%      X COUNTRY %%%%%%%%%%%%%%%%%%%%%%%%%           
% State: a = (mu_x,post_mu_y)
            
% Evaluate Psi(a)
Psi_1 = 1./(1 + (1 + param.tau)^(param.theta/(1 - param.theta)) ... 
    .*(c_g1./c_g2).^(1/(1 - param.theta)));
for a = 1:length(state_x)
    % Nodes and weights
    [x_nodes, x_weights] = qnwnorm([param.N_quad param.N_quad], ...
    	[state_x(a,2) (param.s_x^(-2)/param.p_eta_x^(-2)*param.m_x + state_x(a,1))/(param.s_x^(-2)/param.p_eta_x^(-2)+1)], ...
        [param.post_s_y^2 0;
        0 param.second_s_x^2]) ;
    % Evaluate Gam at integration nodes and construct q
    h1_x    = funeval(c_h1,fspace_y,x_nodes); h1_x = real(h1_x);
    h2_x    = funeval(c_h2,fspace_y,x_nodes); h2_x = real(h2_x);
    Gam_1_x = 1./(1 + (1 + param.tau)^(param.theta/(1 - param.theta)) ...
        .*(h1_x./h2_x).^(1/(1 - param.theta))); Gam_1_x = real(Gam_1_x);
    q_x     = exp(state_x(a,1) - x_nodes(:,1)).*exp(0.5*(param.sig_x^2 - param.sig_y^2)).*Psi_1(a)./Gam_1_x; q_x = real(q_x);
    p_x     = 1./q_x;
    % Conditional expectations
    Psi_2        = 1/(1 - param.sigma)*( (1 - Psi_1(a)).^param.theta ...
        + (Psi_1(a)./((1 + param.tau).*q_x)).^param.theta ).^((1 - param.sigma)/param.theta);
    exp_x_1(a,1) = x_weights'*Psi_2.^((1 - param.sigma)*(1 - param.sigma - param.theta));
    exp_x_2(a,1) = x_weights'*(Psi_2.^((1 - param.sigma)*(1 - param.sigma - param.theta)).*q_x.^(-param.theta));
    exp_qx(a,1)  = x_weights'*q_x;
    var_qx(a,1)  = x_weights'*(q_x.^2) - exp_qx(a,1).^2;
    cv_qx(a,1)   = sqrt(var_qx(a,1))./exp_qx(a,1);
    exp_px(a,1)  = x_weights'*p_x;
    var_px(a,1)  = x_weights'*(p_x.^2) - exp_px(a,1).^2;
    cv_px(a,1)   = sqrt(var_px(a,1))./exp_px(a,1);
end
    % Update basis coefficients
    c_g1_new = (1 - param.relax)*funfitxy(fspace_x,state_x,exp_x_1) + param.relax*c_g1;
    c_g2_new = (1 - param.relax)*funfitxy(fspace_x,state_x,exp_x_2) + param.relax*c_g2;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%      Y COUNTRY    %%%%%%%%%%%%%%%%%%
% State: a=(mu_y,post_mu_x)
            
% Evaluate Gam(a)
Gam_1 =  1./(1 + (1 + param.tau)^(param.theta/(1-param.theta)) ...
    .*(c_h1./c_h2).^(1/(1 - param.theta)));
for a = 1:length(state_y)
    [y_nodes, y_weights] = qnwnorm([param.N_quad param.N_quad], ...
    	[state_y(a,2) (param.s_y^(-2)/param.p_eta_y^(-2)*param.m_y + state_y(a,1))/(param.s_y^(-2)/param.p_eta_y^(-2) + 1)], ...
        [param.post_s_x^2 0;
        0 param.second_s_y^2]);
    % Evaluate Psi(integration nodes) and construct q
    g1_y    = funeval(c_g1,fspace_x,y_nodes); g1_y = real(g1_y);
    g2_y    = funeval(c_g2,fspace_x,y_nodes); g2_y = real(g2_y);
    Psi_1_y = 1./(1 + (1 + param.tau)^(param.theta/(1 - param.theta)) ...
        .*(g1_y./g2_y).^(1/(1 - param.theta))); Psi_1_y = real(Psi_1_y);                
    q_y     = exp(y_nodes(:,1) - state_y(a,1)).*exp(0.5*(param.sig_x^2 - param.sig_y^2)).* Psi_1_y./Gam_1(a); q_y = real(q_y);
    p_y = 1./q_y;
    % Conditional expectations
    Gam_2        = 1/(1 - param.sigma)*( (1 - Gam_1(a)).^param.theta ...
        + (Gam_1(a).*q_y./(1 + param.tau)).^param.theta ).^((1 - param.sigma)/param.theta);
    exp_y_1(a,1) = y_weights'* Gam_2.^((1 - param.sigma)*(1 - param.sigma - param.theta));
    exp_y_2(a,1) = y_weights'*(Gam_2.^((1 - param.sigma)*(1 - param.sigma - param.theta)).*q_y.^param.theta);
    exp_qy(a,1)  = y_weights'*q_y;
    var_qy(a,1)  = y_weights'*q_y.^2 - exp_qy(a,1).^2 ;
    cv_qy(a,1)   = sqrt(var_qy(a,1))./ exp_qy(a,1);
    exp_py(a,1)  = y_weights'*p_y;
    var_py(a,1)  = y_weights'*p_y.^2 - exp_py(a,1).^2 ;
    cv_py(a,1)   = sqrt(var_py(a,1))./ exp_py(a,1);
end
    % Update basis coefficients
    c_h1_new = (1 - param.relax)*funfitxy(fspace_y,state_y,exp_y_1) + param.relax*c_h1;
    c_h2_new = (1 - param.relax)*funfitxy(fspace_y,state_y,exp_y_2) + param.relax*c_h2;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%      CHECK CONVERGENCE    %%%%%%%%%%%%%%%%%%%%%
            
g1_new    = c_g1_new;
g2_new    = c_g2_new;
Psi_1_new =  1./(1+(1 + param.tau)^(param.theta/(1 - param.theta)) ...
    .*(g1_new./g2_new).^(1/(1 - param.theta)));
h1_new    = c_h1_new;
h2_new    = c_h2_new;
Gam_1_new =  1./(1 + (1 + param.tau)^(param.theta/(1 - param.theta)) ...
    .*(h1_new./h2_new).^(1/(1 - param.theta)));

diff = max(abs(Psi_1_new - Psi_1) + abs(Gam_1_new - Gam_1));
i    = i+1;       
end
        
coef    = [c_g1,c_g2,c_h1,c_h2];      
Ep      = exp_px;
Vp      = var_px;
Cp      = cv_px;
Exports = exp(state_x(:,1)+ 0.5*param.sig_x^2).*Psi_1;

