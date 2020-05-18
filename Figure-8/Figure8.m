%{
Can Global Uncertainty Promote International Trade?
Authors: Isaac Baley, Laura Veldkamp, and Mike Waugh 
Data: May 2020

---------------------------------------------------------------------------
---------------------------------------------------------------------------
            FIGURE 8: Completing the Market Reduces Exports

Reads results  "results_financial.mat" and plots Figures 8.  
---------------------------------------------------------------------------
---------------------------------------------------------------------------
%}

load('results_financial.mat')

font_l = 22;
font_s = 20;
set(0,'DefaultTextinterpreter','latex')

stats_dom = real(stats_dom); 

figure
subplot(1,2,1)
plot(alpha_vec,squeeze(stats_dom(6,1,:)),'color',[0.5,0.5,0.5],'linewidth',4)
title('\textbf{Aggregate exports}','fontname','times','fontsize',20)
set(gca,'box','on','fontname','times','fontsize',font_l)
xlabel('Mass with contingent exports $\alpha$','fontname','times','fontsize',font_l)
axis tight


subplot(1,2,2)
plot(alpha_vec,squeeze(stats_dom(4,1,:))./(1-alpha_vec'),'k--','linewidth',4)
hold on
plot(alpha_vec,squeeze(stats_dom(5,1,:))./alpha_vec','k','linewidth',4)
title('\textbf{Exports by type}','fontname','times','fontsize',20)
set(gca,'box','on','fontname','times','fontsize',font_l)
legend('Non-contingent exports','Contingent exports') 
set(gca,'box','on','fontname','times','fontsize',font_l)
legend boxoff
xlabel('Mass with contingent exports $\alpha$','fontname','times','fontsize',font_l)
axis tight

% Saves figure
saveas(gcf,'fig8_financial.eps','epsc')

