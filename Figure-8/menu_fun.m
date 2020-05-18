%Function menu for main.m

function [out1 out2 out3 out4 out5 out6 out7 out8 out9 out10 out11 out12 out13 out14 out15 out16 out17 out18 out19 out20 out21] = menu_fun(flag,param,grid,fspace_x,fspace_y,inp1,inp2,inp3,inp4)

%out  = outputs
%inp  = inputs
%flag = chooses function to evaluate

switch flag
    
    case 'fundamentals long (mu_x,mu_y)'                                    % Uses true foreign productivity
        out2 = exp(grid.mu_x + 0.5*param.sig_x^2) ;                         % Domestic fundamental f_x
        out3 = exp(grid.mu_y + 0.5*param.sig_y^2) ;                         % Foreign fundamental f_y
        out1 = gridmake(out2,out3);                                         % organized as [f_x,f_y]
        out1 = out1(:,1)./out1(:,2);                                        % f_x/f_y
        
    case 'fundamentals long (mu_y,mu_x)'                                    % Uses true foreign productivity
        out2 = exp(grid.mu_x + 0.5*param.sig_x^2) ;                         % Domestic fundamental f_x
        out3 = exp(grid.mu_y + 0.5*param.sig_y^2) ;                         % Foreign fundamental f_y
        out1 = gridmake(out3,out2);                                         % organized as [f_y,f_x]
        out1 = out1(:,2)./out1(:,1);                                        % f_x/f_y
        
    case 'index_pi'
        out1  = 1./(1 + (inp1/param.beta).^(param.theta/(1-param.theta))) ;
        
    case 'perfect info with bias'
        
        q_x = fsolve(@(x)                                           ...     % Price with state (mu_x,mu_y)
            x-menu_fun('fundamentals long (mu_x,mu_y)',param,grid)  ...
            .*menu_fun('index_pi',param,[],[],[],x)                 ...
            ./menu_fun('index_pi',param,[],[],[],(1./x)),           ...
            ones(param.N_mu_x*param.N_mu_y,1));
        
        q_y = fsolve(@(x)                                           ...     % Price with state (mu_y,mu_x)
            x-menu_fun('fundamentals long (mu_y,mu_x)',param,grid) ...
            .*menu_fun('index_pi',param,[],[],[],x)                ...
            ./menu_fun('index_pi',param,[],[],[],(1./x)),          ...
            ones(param.N_mu_x*param.N_mu_y,1));
        
        out1 = menu_fun('index_pi',param,[],[],[],q_x) ;                    % Psi_pi
        out2 = menu_fun('index_pi',param,[],[],[],1./q_y) ;                 % Gam_pi
        out3 = q_x;                                                         % Price q_x=q(\mu_x,\mu_y)
        out4 = q_y;                                                         % Price q_y=q(\mu_y,\mu_x)
        
    case 'psi 2'
        param.beta = 1;
        out1 = ((1-inp1).^param.theta + ...                                 %inp1 = Psi or Gam
            (inp1.*param.beta./inp2).^param.theta).^(1/param.theta);        %inp2 = q   or 1/q
        
    case 'perfect info coef'
        param_names;
        state_x = gridmake(funnode(fspace_x));
        state_y = gridmake(funnode(fspace_y));
        Psi_input = inp1;
        Gam_input = inp2;
        qx_input  = inp3;
        qy_input  = inp4;
                
        %1)Construct arguments of conditional expectations
        Psi_2_input   = menu_fun('psi 2',param,[],[],[],Psi_input,qx_input);
        Gam_2_input   = menu_fun('psi 2',param,[],[],[],Gam_input,qy_input);
        
        %2) Conditional expectations (with perfect info, same as functions)
        exp_x_1_input = Psi_2_input.^(1-theta);
        exp_x_2_input = exp_x_1_input.*qx_input.^(-theta);
        exp_y_1_input = Gam_2_input.^(1-theta);
        exp_y_2_input = exp_y_1_input.*qy_input.^(theta);
        
        %3) Compute basis coefficients (initial coefficient guess)
        c_g1 = funfitxy(fspace_x,state_x,exp_x_1_input);
        c_g2 = funfitxy(fspace_x,state_x,exp_x_2_input);
        c_h1 = funfitxy(fspace_y,state_y,exp_y_1_input);
        c_h2 = funfitxy(fspace_y,state_y,exp_y_2_input);
        out1 = [c_g1,c_g2,c_h1,c_h2];
        
        %4) Check
        g1    = funeval(c_g1,fspace_x,state_x);
        g2    = funeval(c_g2,fspace_x,state_x);
        Psi_check = 1./(1+(1+tau)^(theta/(1-theta)).*(g1./g2).^(1/(1-theta)));
        h1    = funeval(c_h1,fspace_y,state_y);
        h2    = funeval(c_h2,fspace_y,state_y);
        Gam_check = 1./(1+(1+tau)^(theta/(1-theta)).*(h1./h2).^(1/(1-theta)));
        
        
    case 'signal'
        
        coef       = inp1;
        c_g1_new   = coef(:,1);
        c_g2_new   = coef(:,2);
        c_h1_new   = coef(:,3);
        c_h2_new   = coef(:,4);
        param_names;
        state_x    = gridmake(funnode(fspace_x));
        state_y    = gridmake(funnode(fspace_y));
        exp_x_1    = zeros(length(state_x),1);
        exp_x_2    = zeros(length(state_x),1);
        exp_qx     = zeros(length(state_x),1);
        var_qx     = zeros(length(state_x),1);
        cv_qx      = zeros(length(state_x),1);
        skew_qx    = zeros(length(state_x),1);
        exp_invqx  = zeros(length(state_x),1);
        var_invqx  = zeros(length(state_x),1);
        cv_invqx   = zeros(length(state_x),1);
        exp_y_1    = zeros(length(state_y),1);
        exp_y_2    = zeros(length(state_y),1);
        exp_qy     = zeros(length(state_y),1);
        var_qy     = zeros(length(state_y),1);
        cv_qy      = zeros(length(state_y),1);
        skew_qy    = zeros(length(state_y),1);
        
        % 	PAU UPDATE
        
        I_x  = zeros(length(state_x),1);
        II_1_x  = zeros(length(state_x),1);
        II_2_x  = zeros(length(state_x),1);
        II_x  = zeros(length(state_x),1);
        III_x  = zeros(length(state_x),1);
       
        I_y  = zeros(length(state_y),1);
        II_1_y  = zeros(length(state_y),1);
        II_2_y  = zeros(length(state_y),1);
        II_y  = zeros(length(state_y),1);
        III_y  = zeros(length(state_y),1);
        
        
        % Fixed point algorithm
        diff  = 100;
        i=0;
        
        while(diff>tol && i<1000)
            c_g1 = c_g1_new;
            c_g2 = c_g2_new;
            c_h1 = c_h1_new;
            c_h2 = c_h2_new;
            
            %%%%%%%%%%%%%%%%%%%%%%      X COUNTRY %%%%%%%%%%%%%%%%%%%%%%%%%
            
            % State: a = (mu_x,post_mu_y)
            
            % Evaluate Psi(a)
            g1    = funeval(c_g1,fspace_x,state_x);
            g2    = funeval(c_g2,fspace_x,state_x);
            Psi_1 = 1./(1+(1+tau)^(theta/(1-theta)).*(g1./g2).^(1/(1-theta)));
            
            for a=1:length(state_x)
                % Nodes and weights
                [x_nodes x_weights] = qnwnorm([N_quad N_quad], ...
                    [state_x(a,2) (s_x^(-2)/p_eta_x^(-2)*m_x + state_x(a,1))/(s_x^(-2)/p_eta_x^(-2)+1)], ...
                    [post_s_y^2          0      ;
                    0      second_s_x^2]) ;
                % Evaluate Gam at integration nodes and construct q
                h1_x    = funeval(c_h1,fspace_y,x_nodes);
                h2_x    = funeval(c_h2,fspace_y,x_nodes);
                Gam_1_x =  1./(1+(1+tau)^(theta/(1-theta)).*(h1_x./h2_x).^(1/(1-theta)));
                q_x     = exp(state_x(a,1)-x_nodes(:,1)).*exp(0.5*(sig_x^2-sig_y^2)).* Psi_1(a)./Gam_1_x;
                % Conditional expectations
                Psi_2   = menu_fun('psi 2',param,[],[],[],Psi_1(a),q_x);
                exp_x_1(a,1) = x_weights'* Psi_2.^(1-theta);
                exp_x_2(a,1) = x_weights'*(Psi_2.^(1-theta).*q_x.^(-theta));
                exp_qx(a,1)  = x_weights'*q_x;
                var_qx(a,1)  = x_weights'*(q_x.^2) - exp_qx(a,1).^2;
                cv_qx(a,1)   = sqrt(var_qx(a,1))./ exp_qx(a,1);
                skew_qx(a,1) = x_weights'*(((q_x-exp_qx(a,1))./sqrt(var_qx(a,1))).^3);
                
                exp_invqx(a,1)  = x_weights'*(1./q_x);
                var_invqx(a,1)  = x_weights'*((1./q_x).^2) - exp_invqx(a,1).^2 ;
                cv_invqx(a,1)   = sqrt(var_invqx(a,1))./ exp_invqx(a,1);
                
                
            % PAU UPDATE (Country X)
                
                C_x = ((1-Psi_1(a)).^theta + (Psi_1(a)./q_x).^theta).^(1/theta);
                I_x(a,1) = x_weights'*(q_x.^(-theta));
                III_x(a,1) = x_weights'*(C_x.^(1-theta));
                II_1_x(a,1) = x_weights'*(C_x.^(1-theta).*q_x.^(-theta));
                II_2_x(a,1) = III_x(a,1).*I_x(a,1);
                II_x(a,1)  = II_1_x(a,1) - II_2_x(a,1);
                
                
                
                
            
            end
            % Update basis coefficients
            c_g1_new = (1-relax)*funfitxy(fspace_x,state_x,exp_x_1)+relax*c_g1;
            c_g2_new = (1-relax)*funfitxy(fspace_x,state_x,exp_x_2)+relax*c_g2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%      Y COUNTRY    %%%%%%%%%%%%%%%%%%
            
            % State: a=(mu_y,post_mu_x)
            
            % Evaluate Gam(a)
            h1    = funeval(c_h1,fspace_y,state_y);
            h2    = funeval(c_h2,fspace_y,state_y);
            Gam_1 =  1./(1+(1+tau)^(theta/(1-theta)).*(h1./h2).^(1/(1-theta)));
            
            for a=1:length(state_y)
                [y_nodes y_weights] = qnwnorm([N_quad N_quad], ...
                    [state_y(a,2) (s_y^(-2)/p_eta_y^(-2)*m_y + state_y(a,1))/(s_y^(-2)/p_eta_y^(-2)+1)], ...
                    [post_s_x^2          0      ;
                    0      second_s_y^2]);
                % Evaluate Psi(integration nodes) and construct q
                g1_y    = funeval(c_g1,fspace_x,y_nodes);
                g2_y    = funeval(c_g2,fspace_x,y_nodes);
                Psi_1_y =  1./(1+(1+tau)^(theta/(1-theta)).*(g1_y./g2_y).^(1/(1-theta)));
                q_y     = exp(y_nodes(:,1) - state_y(a,1)).*exp(0.5*(sig_x^2-sig_y^2)).* Psi_1_y./Gam_1(a);
                % Conditional expectations
                Gam_2   = menu_fun('psi 2',param,[],[],[],Gam_1(a),1./q_y);
                exp_y_1(a,1) = y_weights'* Gam_2.^(1-theta);
                exp_y_2(a,1) = y_weights'*(Gam_2.^(1-theta).*q_y.^theta);
                exp_qy(a,1)   = y_weights'*q_y;
                var_qy(a,1)   = y_weights'*q_y.^2 - exp_qy(a,1).^2 ;
                cv_qy(a,1)    = sqrt(var_qy(a,1))./ exp_qy(a,1);
                skew_qy(a,1) = y_weights'*(((q_y-exp_qy(a,1))./sqrt(var_qy(a,1))).^3);
                
                
                 % PAU UPDATE (Country Y)
                
                C_y = ((1-Gam_1(a)).^theta + (Gam_1(a).*q_y).^theta).^(1/theta);
                I_y(a,1) = y_weights'*(q_y.^(theta));
                III_y(a,1) = y_weights'*(C_y.^(1-theta));
                II_1_y(a,1) = y_weights'*(C_y.^(1-theta).*q_y.^(theta));
                II_2_y(a,1) = III_y(a,1).*I_y(a,1);
                II_y(a,1)  = II_1_y(a,1) - II_2_y(a,1);
                
          
                
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
            
            %diff1 = max(abs(c_g1-c_g1_new))+max(abs(c_g2-c_g2_new))...
            %+max(abs(c_h1-c_h1_new))+max(abs(c_h2-c_h2_new));
            i=i+1;
            
        end
        
        coef = [c_g1,c_g2,c_h1,c_h2];
        out1 = Psi_1;
        out2 = Gam_1;
        out3 = coef;
        out4 = exp_qx; %conditional expectations of price q
        out5 = var_qx; %conditional variance of price q
        out6 = cv_qx;  %conditional coefficient of variation of price q
        out7 = exp_invqx; %conditional expectations of price q
        out8 = var_invqx; %conditional variance of price q
        out9 = cv_invqx;  %conditional coefficient of variation of price q
        % PAU UPDATE
        out10 = I_x;
        out11 = II_1_x;
        out12 = II_2_x;
        out13 = III_x;
        out14 = II_x;
        out15 = skew_qx;
        
        out16 = I_y;
        out17 = II_1_y;
        out18 = II_2_y;
        out19 = III_y;
        out20 = II_y;
        out21 = skew_qy;
        
    case 'fixed belief'
 
        
        coef = inp1;
        c_g1_new = coef(:,1);
        c_g2_new = coef(:,2);
      	param_names;
        state_x = gridmake(funnode(fspace_x));
        state_y    = gridmake(funnode(fspace_y));
        exp_x_1 = zeros(length(state_x),1);
        exp_x_2 = zeros(length(state_x),1);

        
        exp_qx     = zeros(length(state_x),1);
        var_qx     = zeros(length(state_x),1);
        cv_qx      = zeros(length(state_x),1);
        exp_invqx  = zeros(length(state_x),1);
        var_invqx  = zeros(length(state_x),1);
        cv_invqx   = zeros(length(state_x),1);
        exp_y_1    = zeros(length(state_y),1);
        exp_y_2    = zeros(length(state_y),1);
        exp_qy     = zeros(length(state_y),1);
        var_qy     = zeros(length(state_y),1);
        cv_qy      = zeros(length(state_y),1);
        skew_qx    = zeros(length(state_x),1);
        skew_qy    = zeros(length(state_y),1);
           
        
        I_x  = zeros(length(state_x),1);
        II_1_x  = zeros(length(state_x),1);
        II_2_x  = zeros(length(state_x),1);
        II_x  = zeros(length(state_x),1);
        III_x  = zeros(length(state_x),1);
        
        I_y  = zeros(length(state_y),1);
        II_1_y  = zeros(length(state_y),1);
        II_2_y  = zeros(length(state_y),1);
        II_y  = zeros(length(state_y),1);
        III_y  = zeros(length(state_y),1);
       
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%  Y COUNTRY  FIXED  %%%%%%%%%%%%%%%%%%
                        
       % Evaluate Gam(a)
       c_g1_fixed = coef(:,1);
       c_g2_fixed = coef(:,2);
       c_h1 = coef(:,3);
       c_h2 = coef(:,4);
%       state_y = gridmake(funnode(fspace_y));
       h1    = funeval(c_h1,fspace_y,state_y);
       h2    = funeval(c_h2,fspace_y,state_y);
       Gam_1 =  1./(1+(1+tau)^(theta/(1-theta)).*(h1./h2).^(1/(1-theta)));
       
            for a=1:length(state_y)
                [y_nodes y_weights] = qnwnorm([N_quad N_quad], ...
                    [state_y(a,2) (s_y^(-2)/p_eta_y^(-2)*m_y + state_y(a,1))/(s_y^(-2)/p_eta_y^(-2)+1)], ...
                    [post_s_x^2          0      ;
                    0      second_s_y^2]);
                % Evaluate Psi(integration nodes) and construct q
                g1_y_fixed    = funeval(c_g1_fixed,fspace_x,y_nodes);
                g2_y_fixed    = funeval(c_g2_fixed,fspace_x,y_nodes);
                Psi_1_y =  1./(1+(1+tau)^(theta/(1-theta)).*(g1_y_fixed./g2_y_fixed).^(1/(1-theta)));
                q_y     = exp(y_nodes(:,1) - state_y(a,1)).*exp(0.5*(sig_x^2-sig_y^2)).* Psi_1_y./Gam_1(a);
                 % Conditional expectations
                 Gam_2   = menu_fun('psi 2',param,[],[],[],Gam_1(a),1./q_y);
                 exp_y_1(a,1) = y_weights'* Gam_2.^(1-theta);
                 exp_y_2(a,1) = y_weights'*(Gam_2.^(1-theta).*q_y.^theta);
                 exp_qy(a,1)   = y_weights'*q_y;
                 var_qy(a,1)   = y_weights'*q_y.^2 - exp_qy(a,1).^2 ;
                 cv_qy(a,1)    = sqrt(var_qy(a,1))./ exp_qy(a,1);
                 skew_qy(a,1) = y_weights'*(((q_y-exp_qy(a,1))./sqrt(var_qy(a,1))).^3);
                 
                 
                 % PAU UPDATE (Country Y)
                
                C_y = ((1-Gam_1(a)).^theta + (Gam_1(a).*q_y).^theta).^(1/theta);
                I_y(a,1) = y_weights'*(q_y.^(theta));
                III_y(a,1) = y_weights'*(C_y.^(1-theta));
                II_1_y(a,1) = y_weights'*(C_y.^(1-theta).*q_y.^(theta));
                II_2_y(a,1) = III_y(a,1).*I_y(a,1);
                II_y(a,1)  = II_1_y(a,1) - II_2_y(a,1);
                
          
                
            end
       
       
        %%%%%%%%%%%%%%%%%%%%%%   X COUNTRY ITERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%

        
        % Fixed point algorithm
        diff  = 100;
        i=0;
        
        
        while(diff>tol && i<1000)
            c_g1 = c_g1_new;
            c_g2 = c_g2_new;


            % State: a = (mu_x,post_mu_y)
            
            % Evaluate Psi(a)
            g1    = funeval(c_g1,fspace_x,state_x);
            g2    = funeval(c_g2,fspace_x,state_x);
            Psi_1 = 1./(1+(1+tau)^(theta/(1-theta)).*(g1./g2).^(1/(1-theta)));
            
            for a=1:length(state_x)
                % Nodes and weights
                [x_nodes x_weights] = qnwnorm([N_quad N_quad], ...
                    [state_x(a,2) (s_x^(-2)/p_eta_x^(-2)*m_x + state_x(a,1))/(s_x^(-2)/p_eta_x^(-2)+1)], ...
                    [post_s_y^2          0      ;
                    0      second_s_x^2]) ;
                % Evaluate Gam at integration nodes and construct q
                h1_x    = funeval(c_h1,fspace_y,x_nodes);
                h2_x    = funeval(c_h2,fspace_y,x_nodes);
                Gam_1_x =  1./(1+(1+tau)^(theta/(1-theta)).*(h1_x./h2_x).^(1/(1-theta)));
                q_x     = exp(state_x(a,1)-x_nodes(:,1)).*exp(0.5*(sig_x^2-sig_y^2)).* Psi_1(a)./Gam_1_x;
                % Conditional expectations
                Psi_2   = menu_fun('psi 2',param,[],[],[],Psi_1(a),q_x);
                exp_x_1(a,1) = x_weights'* Psi_2.^(1-theta);
                exp_x_2(a,1) = x_weights'*(Psi_2.^(1-theta).*q_x.^(-theta));
                exp_qx(a,1)  = x_weights'*q_x;
                var_qx(a,1)  = x_weights'*(q_x.^2) - exp_qx(a,1).^2 ;
                cv_qx(a,1)   = sqrt(var_qx(a,1))./ exp_qx(a,1);
                skew_qx(a,1) = x_weights'*(((q_x-exp_qx(a,1))./sqrt(var_qx(a,1))).^3);
                
                exp_invqx(a,1)  = x_weights'*(1./q_x);
                var_invqx(a,1)  = x_weights'*((1./q_x).^2) - exp_invqx(a,1).^2 ;
                cv_invqx(a,1)   = sqrt(var_invqx(a,1))./ exp_invqx(a,1);
   
                % PAU UPDATE (Country X)
                
                C_x = ((1-Psi_1(a)).^theta + (Psi_1(a)./q_x).^theta).^(1/theta);
                I_x(a,1) = x_weights'*(q_x.^(-theta));
                III_x(a,1) = x_weights'*(C_x.^(1-theta));
                II_1_x(a,1) = x_weights'*(C_x.^(1-theta).*q_x.^(-theta));
                II_2_x(a,1) = III_x(a,1).*I_x(a,1);
                II_x(a,1)  = II_1_x(a,1) - II_2_x(a,1);
                
            end
            % Update basis coefficients
            c_g1_new = (1-relax)*funfitxy(fspace_x,state_x,exp_x_1)+relax*c_g1;
            c_g2_new = (1-relax)*funfitxy(fspace_x,state_x,exp_x_2)+relax*c_g2;
                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%      CHECK CONVERGENCE    %%%%%%%%%%%%%%%%%%%%%
            
            g1_new    = funeval(c_g1_new,fspace_x,state_x);
            g2_new    = funeval(c_g2_new,fspace_x,state_x);
            Psi_1_new =  1./(1+(1+tau)^(theta/(1-theta)).*(g1_new./g2_new).^(1/(1-theta)));
             
            diff = max(abs(Psi_1_new-Psi_1));
            i=i+1;
            
        end
        
        coef = [c_g1,c_g2,c_h1,c_h2];
        out1 = Psi_1;
        out2 = Gam_1;
        out3 = coef;
        out4 = exp_qx; %conditional expectations of price q
        out5 = var_qx; %conditional variance of price q
        out6 = cv_qx;  %conditional coefficient of variation of price q
        out7 = exp_invqx; %conditional expectations of price q
        out8 = var_invqx; %conditional variance of price q
        out9 = cv_invqx;  %conditional coefficient of variation of price q
        % PAU UPDATE
        out10 = I_x;
        out11 = II_1_x;
        out12 = II_2_x;
        out13 = III_x;
        out14 = II_x;
        out15 = skew_qx;
                
        out16 = I_y;
        out17 = II_1_y;
        out18 = II_2_y;
        out19 = III_y;
        out20 = II_y;
        out21 = skew_qy;
        
        
        
        case 'financial contracts'
        
            
        options=optimset('fsolve');
        options=optimset(options,'display','off'); 
        %options=optimset(options,'LargeScale','off'); 
        %options=optimset(options,'Jacobian','off');
        %options=optimset(options,'TolFun',1e-6);
          
        coef       = inp1;
        alpha      = inp2;
       
        
        c_g1_new   = coef(:,1);
        c_g2_new   = coef(:,2);
        c_h1_new   = coef(:,3);
        c_h2_new   = coef(:,4);
        param_names;
        state_x    = gridmake(funnode(fspace_x));
        state_y    = gridmake(funnode(fspace_y));
        exp_x_1    = zeros(length(state_x),1);
        exp_x_2    = zeros(length(state_x),1);
        exp_qx     = zeros(length(state_x),1);
        var_qx     = zeros(length(state_x),1);
        cv_qx      = zeros(length(state_x),1);
        skew_qx    = zeros(length(state_x),1);
        exp_invqx  = zeros(length(state_x),1);
        var_invqx  = zeros(length(state_x),1);
        cv_invqx   = zeros(length(state_x),1);
        exp_y_1    = zeros(length(state_y),1);
        exp_y_2    = zeros(length(state_y),1);
        exp_qy     = zeros(length(state_y),1);
        var_qy     = zeros(length(state_y),1);
        cv_qy      = zeros(length(state_y),1);
        skew_qy    = zeros(length(state_y),1);
        
        
        
        
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
            [x_nodes x_weights] = qnwnorm(N_quad, m_y, s_y^2) ;
            % Evaluate Gam at integration nodes and construct q
            h1_x    = funeval(c_h1,fspace_y,x_nodes);
            h2_x    = funeval(c_h2,fspace_y,x_nodes);
            Gam_1_x =  1./(1+(1+tau)^(theta/(1-theta)).*(h1_x./h2_x).^(1/(1-theta)));
            
            for a=1:length(state_x)
            %Here we solve for the equilibrium price that includes bothtypes of agents
                
                 q_x = fsolve(@(x) ...
                 x - exp(state_x(a)-x_nodes).*...
                 (alpha./(1+x.^(theta/(1-theta)))+(1-alpha).*Psi_1(a))./(alpha./(1+x.^(-theta/(1-theta)))+(1-alpha).*Gam_1_x),ones(length(x_nodes),1),options);
            % Conditional expectations
                Psi_2   = menu_fun('psi 2',param,[],[],[],Psi_1(a),q_x);
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
            [y_nodes y_weights] = qnwnorm(N_quad, m_x, s_x^2) ;
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
              Gam_2   = menu_fun('psi 2',param,[],[],[],Gam_1(a),1./q_y);
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
            
            %diff1 = max(abs(c_g1-c_g1_new))+max(abs(c_g2-c_g2_new))...
            %+max(abs(c_h1-c_h1_new))+max(abs(c_h2-c_h2_new));
            i=i+1;
            
        end
        
        coef = [c_g1,c_g2,c_h1,c_h2];
        out1 = Psi_1;
        out2 = Gam_1;
        out3 = coef;
        %out4 = exp_qx; %conditional expectations of price q
        %out5 = var_qx; %conditional variance of price q
        %out6 = cv_qx;  %conditional coefficient of variation of price q
        %out7 = exp_invqx; %conditional expectations of price q
        %out8 = var_invqx; %conditional variance of price q
        %out9 = cv_invqx;  %conditional coefficient of variation of price q
  
        
end
