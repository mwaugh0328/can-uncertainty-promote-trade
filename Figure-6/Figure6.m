%{
Can Global Uncertainty Promote International Trade?
Authors: Isaac Baley, Laura Veldkamp, and Mike Waugh 
Data: May 2020
---------------------------------------------------------------------------
---------------------------------------------------------------------------
            FIGURE 6: State-dependent Responses to Uncertainty

Reads results  "results_statedep.mat" and plots Figures 6.  
---------------------------------------------------------------------------
---------------------------------------------------------------------------
%}

load('results_statedep.mat')
 
% Picks low and high domestic endowment to plot state-dependent outcomes
low  = ceil(N*0.3); 
high = ceil(N*0.7); 

% Fixed foreign endowment at the average
med  = ceil(N*0.5); 

% Sets fonts
font_l = 22;
font_s = 18;
set(0,'DefaultTextinterpreter','latex')

%% Three plots

figure(1)

%Plots average terms of trade
subplot(1,3,1)
h1 = plot(eta_vec,reshape(Ep_eta_r(low,med,:), [1 N_eta])./Ep_eta_r(low,med,1),'k','linewidth',6);    % low mu_x
hold on
h2 = plot(eta_vec,reshape(Ep_eta_r(high,med,:), [1 N_eta])./Ep_eta_r(high,med,1),':k','linewidth',6); % high mu_x
title('\textbf{Expected Terms of Trade}','fontsize',18);
set(gca,'box','on','fontname','times','fontsize',font_l)
xlabel('Uncertainty');
ylim([0.99 1.6])
h_leg=legend([h1 h2],'low \mu_{x}', 'high \mu_{x}');
set(h_leg,'box','on','location','Northwest','fontsize',25);
legend boxoff

%Plots volatility of terms of trade
subplot(1,3,2)
plot(eta_vec,reshape(Vp_eta_r(low,med,:), [1 N_eta]),'k','linewidth',6);    % low mu_x
hold on
plot(eta_vec,reshape(Vp_eta_r(high,med,:), [1 N_eta]),':k','linewidth',6);  % high mu_x
title('\textbf{Volatility of Terms of Trade}','fontsize',18);
set(gca,'box','on','fontname','times','fontsize',font_l)
xlabel('Uncertainty');
ylim([0 11])

%Plots export volume
subplot(1,3,3)
plot(eta_vec,reshape(Exports_eta_r(low,med,:), [1 N_eta])./Exports_eta_r(low,med,1),'k','linewidth',6); % low mu_x
hold on
plot(eta_vec,reshape(Exports_eta_r(high,med,:), [1 N_eta])./Exports_eta_r(high,med,1),':k','linewidth',6); % high mu_x
title('\textbf{Export Volume}','fontsize',18);
set(gca,'box','on','fontname','times','fontsize',font_l)
xlabel('Uncertainty');
ylim([0.99 1.05])

%saveas(gcf,'fig6_statedep.eps','epsc')

