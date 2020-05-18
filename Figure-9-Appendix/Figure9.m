%{
Can Global Uncertainty Promote International Trade?
Authors: Isaac Baley, Laura Veldkamp, and Mike Waugh 
Data: May 2020

---------------------------------------------------------------------------
---------------------------------------------------------------------------
            FIGURE 9 Appendix D: Effect of Uncertainty in Trade 
                 for Different Levels of Risk Aversion 

Reads three matrices of results and plots Figures 9:
  - "results_ra_bench.mat",   risk aversion sigma = 1-theta
  - "results_ra_neutral.mat", risk aversion sigma = 0
  - "results_ra_averse.mat", risk aversion sigma = 1.5 

---------------------------------------------------------------------------
---------------------------------------------------------------------------
%}

% Evaluate at average state
med  = ceil(N*0.5);

% Loads results for sigma = 1-theta
load('results_ra_bench.mat','Exports_eta_r','Export_eta_r_1','N','N_eta')
A1   = Exports_eta_r(med,med,:);
B1   = Export_eta_r_1(med,med,:);

% Loads results for sigma = 0
load('results_ra_neutral.mat','Exports_eta_r','Export_eta_r_1')
A2   = Exports_eta_r(med,med,:);
B2   = Export_eta_r_1(med,med,:);

% Loads results for sigma = 1.5
load('results_ra_averse.mat','Exports_eta_r','Export_eta_r_1')
A3   = Exports_eta_r(med,med,:);
B3   = Export_eta_r_1(med,med,:);

% Set fonts
font_l = 22;
font_s = 20;
set(0,'DefaultTextinterpreter','latex')

% Plot figure
figure;
subplot(1,2,1)
plot(reshape(A1,[1 N_eta]),'k','linewidth',4);
hold on
plot(reshape(A2,[1 N_eta]),'--k','linewidth',3);
plot(reshape(real(A3),[1 N_eta]),'--','Color',[0.5 0.5 0.5],'linewidth',3);
title('\textbf{High Substitution}','fontname','times','fontsize',18)
set(gca,'box','on','fontname','times','fontsize',font_l)
legend('\sigma = 1 - \theta = 0.2','\sigma = 0','\sigma = 1.5')
legend boxoff
ylabel('Exports, $T_x$','fontname','times','fontsize',font_l)
xlabel('Uncertainty','fontname','times','fontsize',font_l)
axis tight

subplot(1,2,2)
plot(reshape(B1,[1 N_eta]),'k','linewidth',4);
hold on
plot(reshape(B2,[1 N_eta]),'--k','linewidth',3);
plot(reshape(real(B3),[1 N_eta]),'--','Color',[0.5 0.5 0.5],'linewidth',3);
title('\textbf{Low Substitution}','fontname','times','fontsize',18)
set(gca,'box','on','fontname','times','fontsize',font_l)
legend('\sigma = 1 - \theta = 0.7','\sigma = 0','\sigma = 1.5')
legend boxoff
ylabel('Exports, $T_x$','fontname','times','fontsize',font_l)
xlabel('Uncertainty','fontname','times','fontsize',font_l)
axis tight

saveas(gcf,'fig9_riskaversion.eps','epsc')