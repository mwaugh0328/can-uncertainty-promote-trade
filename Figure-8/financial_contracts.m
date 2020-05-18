function coef = financial_contracts(param,fspace_x,fspace_y,inp1,inp2)

% Set optimization options
options=optimset('fsolve');
options=optimset(options,'display','off');

% Read inputs (coefficients and alpha (fraction of agents with complete markets)
coef       = inp1;
alpha      = inp2;


% Initialize vectors
c_g1_new   = coef(:,1);
c_g2_new   = coef(:,2);
c_h1_new   = coef(:,3);
c_h2_new   = coef(:,4);

% For simplicity, assign parameter names
theta          = param.theta ;                                              % Aggregator parameter for composite good
tau            = param.tau;                                                 % Iceberg cost
m_x            = param.m_x ;                                                % mean  of domestic aggregate productivity
m_y            = param.m_y ;                                                % mean  of foreign  aggregate productivity
s_x            = param.s_x ;                                                % stdev of domestic aggregate productivity
s_y            = param.s_y ;                                                % stdev of foreign  aggregate productivity

% Algorithm parameters
tol            = param.tol; 
N_quad         = param.N_quad;  
relax          = param.relax;

state_x    = gridmake(funnode(fspace_x));
state_y    = gridmake(funnode(fspace_y));
exp_x_1    = zeros(length(state_x),1);
exp_x_2    = zeros(length(state_x),1);
exp_y_1    = zeros(length(state_y),1);
exp_y_2    = zeros(length(state_y),1);

% Fixed point algorithm
diff  = 100;
i=0;

while(diff>tol && i<1000)
    c_g1 = c_g1_new;
    c_g2 = c_g2_new;
    c_h1 = c_h1_new;
    c_h2 = c_h2_new;
    
    %%%%%%%%%%%%%%%%%%%%%%      X COUNTRY %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % State: a = mu_x
    
    % Evaluate Psi(a)
    g1    = funeval(c_g1,fspace_x,state_x);
    g2    = funeval(c_g2,fspace_x,state_x);
    Psi_1 = 1./(1+(1+tau)^(theta/(1-theta)).*(g1./g2).^(1/(1-theta)));
    
    % Nodes and weights
    [x_nodes, x_weights] = qnwnorm(N_quad, m_y, s_y^2) ;
    
    % Evaluate Gam at integration nodes and construct q
    h1_x    = funeval(c_h1,fspace_y,x_nodes);
    h2_x    = funeval(c_h2,fspace_y,x_nodes);
    Gam_1_x =  1./(1+(1+tau)^(theta/(1-theta)).*(h1_x./h2_x).^(1/(1-theta)));
    
    for a=1:length(state_x)
        %Here we solve for the equilibrium price that includes both types of agents
        q_x = fsolve(@(x) ...
            x - exp(state_x(a)-x_nodes).*...
            (alpha./(1+x.^(theta/(1-theta)))+(1-alpha).*Psi_1(a))./(alpha./(1+x.^(-theta/(1-theta)))+(1-alpha).*Gam_1_x),ones(length(x_nodes),1),options);
        % Conditional expectations        
        %Psi_2   = menu_fun('psi 2',param,[],[],[],Psi_1(a),q_x);
        
        Psi_2 = ((1-Psi_1(a)).^param.theta + (Psi_1(a)./q_x).^param.theta).^(1/param.theta); 
        exp_x_1(a,1) = x_weights'* Psi_2.^(1-theta);
        exp_x_2(a,1) = x_weights'*(Psi_2.^(1-theta).*q_x.^(-theta));
    end
    % Update basis coefficients
    c_g1_new = (1-relax)*funfitxy(fspace_x,state_x,exp_x_1)+relax*c_g1;
    c_g2_new = (1-relax)*funfitxy(fspace_x,state_x,exp_x_2)+relax*c_g2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%      Y COUNTRY    %%%%%%%%%%%%%%%%%%
    
    % State: a=mu_y
    
    % Evaluate Gam(a)
    h1    = funeval(c_h1,fspace_y,state_y);
    h2    = funeval(c_h2,fspace_y,state_y);
    Gam_1 =  1./(1+(1+tau)^(theta/(1-theta)).*(h1./h2).^(1/(1-theta)));
    [y_nodes, y_weights] = qnwnorm(N_quad, m_x, s_x^2) ;
    
    % Evaluate Psi(integration nodes) and construct q
    g1_y    = funeval(c_g1,fspace_x,y_nodes);
    g2_y    = funeval(c_g2,fspace_x,y_nodes);
    Psi_1_y =  1./(1+(1+tau)^(theta/(1-theta)).*(g1_y./g2_y).^(1/(1-theta)));
    
    for a=1:length(state_y)
        %Here we solve for the equilibrium price that includes both types of agents
        q_y = fsolve(@(x) ...
            x - exp(y_nodes-state_y(a)).*...
            (alpha./(1+x.^(theta/(1-theta)))+(1-alpha).*Psi_1_y)./(alpha./(1+x.^(-theta/(1-theta)))+(1-alpha).*Gam_1(a)),ones(length(x_nodes),1),options );
        % Conditional expectations
        %Gam_2   = menu_fun('psi 2',param,[],[],[],Gam_1(a),1./q_y);
        Gam_2   = ((1-Gam_1(a)).^param.theta + (Gam_1(a).*q_y).^param.theta).^(1/param.theta); 
        exp_y_1(a,1) = y_weights'* Gam_2.^(1-theta);
        exp_y_2(a,1) = y_weights'*(Gam_2.^(1-theta).*q_y.^theta);
    end
   
        
    % Update basis coefficients
    c_h1_new = (1-relax)*funfitxy(fspace_y,state_y,exp_y_1)+relax*c_h1;
    c_h2_new = (1-relax)*funfitxy(fspace_y,state_y,exp_y_2)+relax*c_h2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%      CHECK CONVERGENCE    %%%%%%%%%%%%%%%%%%%%%
    
    g1_new    = funeval(c_g1_new,fspace_x,state_x);
    g2_new    = funeval(c_g2_new,fspace_x,state_x);
    Psi_1_new =  1./(1+(1+tau)^(theta/(1-theta)).*(g1_new./g2_new).^(1/(1-theta)));
    h1_new    = funeval(c_h1_new,fspace_y,state_y);
    h2_new    = funeval(c_h2_new,fspace_y,state_y);
    Gam_1_new =  1./(1+(1+tau)^(theta/(1-theta)).*(h1_new./h2_new).^(1/(1-theta)));
    
    diff = max(abs(Psi_1_new-Psi_1)+abs(Gam_1_new-Gam_1));
    i=i+1;
    
end

coef = [c_g1,c_g2,c_h1,c_h2];