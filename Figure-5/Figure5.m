%{
Can Global Uncertainty Promote International Trade?
Authors: Isaac Baley, Laura Veldkamp, and Mike Waugh 
Data: May 2020

---------------------------------------------------------------------------
---------------------------------------------------------------------------
                             FIGURES 5 
          Plots average foreign beliefs for low and high realizations 
          of domestic endowment, for various levels of signal noise. 

---------------------------------------------------------------------------
---------------------------------------------------------------------------

%}

close all
clear all
clc

% Parameters for domestic Country X
param.m_x   = 0;            % mean aggregate productivity country x;
param.sig_x = sqrt(2);      % st.d. idiosincratic productivity country x;
param.s_x   = 1;            % st.d. aggregate productivity country x;
param.mu_x_min      = param.m_x - 2*param.s_x; 
param.mu_x_max      = param.m_x + 2*param.s_x;  

% Precision grid
N_eta   = 21;
eta_min = 0.0001;
eta_max = 3;
eta_vec = linspace(eta_min,eta_max,N_eta);

%% Create posterior beliefs of foreign Country Y
mu_tilde_h = param.mu_x_max;
mu_hat_h   = (param.sig_x^(-2)*param.m_x + eta_vec.^(-2).*mu_tilde_h)./(param.sig_x^(-2) + eta_vec.^(-2));
mu_tilde_l = param.mu_x_min;
mu_hat_l   = (param.sig_x^(-2)*param.m_x + eta_vec.^(-2).*mu_tilde_l)./(param.sig_x^(-2) + eta_vec.^(-2));

%% Plot
font_l = 20;
font_s = 18;
set(0,'DefaultTextinterpreter','latex')

figure(1);
a1 = plot(eta_vec,param.mu_x_min*ones(length(eta_vec),1),'-','color',[0.5 0.5 0.5],'linewidth',4);
hold on
a2 = plot(eta_vec,param.mu_x_max*ones(length(eta_vec),1),'color',[0.7 0.7 0.7],'linewidth',4);
a3 = plot(eta_vec,param.m_x*ones(length(eta_vec),1),'k','linewidth',4);
a4 = plot(eta_vec,mu_hat_h,':','color',[0.7 0.7 0.7],'linewidth',4);
a5 = plot(eta_vec,mu_hat_l,':','color',[0.5 0.5 0.5],'linewidth',4);
set(gca,'box','on','fontname','times','fontsize',font_s)
xlabel('Uncertainty','fontname','times','fontsize',font_s);
ylim([-3 3])

% Create textbox
annotation('textbox',...
    [0.424023755472674 0.813445378151262 0.20367204475434 0.0474789915966429],...
    'String','High endowment \mu_x^{high}','FontSize',18,'FontName',...
    'Times New Roman','FitBoxToText','off','EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.428896465055942 0.168067226890756 0.209014999189234 0.0596638655462189],...
    'String','Low endowment \mu_x^{low}','FontSize',18,'FontName',...
    'Times New Roman','FitBoxToText','off','EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.449327793092262 0.53109243697479 0.136500000000003 0.0480392156862781],...
    'String','Prior m_x','FontSize',18,'FontName',...
    'Times New Roman','FitBoxToText','off','EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.638905705225714 0.305882352941177 0.212690476190477 0.0512605042016837],...
    'String','Foreign belief \mu_x|signal','FontSize',18,'FontName',...
    'Times New Roman','FitBoxToText','off','EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.640060082981357 0.66386554621849 0.213824263038549 0.0546218487394988],...
    'String','Foreign belief \mu_x|signal','FontSize',18,'FontName',...
    'Times New Roman','FitBoxToText','off','EdgeColor',[1 1 1]);

% Create arrows
annotation('arrow',[0.7105561861521 0.677639046538025],...
    [0.67463025210084 0.615126050420168]);
annotation('arrow',[0.713961407491487 0.677639046538025],...
    [0.353621848739499 0.418487394957983]);


saveas(gcf,'fig5_foreignbeliefs.eps','epsc')

