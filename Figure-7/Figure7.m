%{
Can Global Uncertainty Promote International Trade?
Authors: Isaac Baley, Laura Veldkamp, and Mike Waugh 
Data: May 2020

---------------------------------------------------------------------------
---------------------------------------------------------------------------
            FIGURE 7: Uncertainty Decreases Trade Coordination 
                       and Increases Utility Correlation

Reads results  "results_risksharing.mat" and plots Figures 7.  
---------------------------------------------------------------------------
---------------------------------------------------------------------------
%}


load('results_risksharing.mat')


%% Figure
font_l = 22;
font_s = 20;
set(0,'DefaultTextinterpreter','latex')

figure;
subplot(1,2,1)
plot(eta_vec,real(corr_T),'k','linewidth',4);
title('\textbf{Export Correlation}','fontname','times','fontsize',20)
set(gca,'box','on','fontname','times','fontsize',font_l)
xlabel('Uncertainty','fontname','times','fontsize',font_l)
ylabel('corr ($T_x$, $T_y^*$)','fontname','times','fontsize',font_l)
title('\textbf{Trade correlation}')
xlim([0 3])
ylim([0 0.32])


subplot(1,2,2)
plot(eta_vec,real(corr_C),'k','linewidth',4);
title('\textbf{Utility Correlation}','fontname','times','fontsize',20)
set(gca,'box','on','fontname','times','fontsize',font_l)
xlabel('Uncertainty','fontname','times','fontsize',font_l)
ylabel('corr ($U$, $U^*$)','fontname','times','fontsize',font_l)
title('\textbf{Utility correlation}')
xlim([0 3])
ylim([0.76 0.87])

%saveas(gcf,'fig7_risksharing.eps','epsc')
