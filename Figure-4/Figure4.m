%{
Can Global Uncertainty Promote International Trade?
Authors: Isaac Baley, Laura Veldkamp, and Mike Waugh 
Data: May 2020

---------------------------------------------------------------------------
---------------------------------------------------------------------------
            FIGURE 4: Effect of Uncertainty in Trade 
             Depends on the Elasticity of Substitution

Reads results  "results_twoelasticities.mat" and plots Figures 4.  
---------------------------------------------------------------------------
---------------------------------------------------------------------------

%}


load('results_twoelasticities.mat')

%Picks average state
med  = ceil(N*0.5);

%Sets fonts
font_l = 22;
font_s = 18;
set(0,'DefaultTextinterpreter','latex')

figure;
subplot(1,2,1)
plot(eta_vec,reshape(Exports_eta_r(med,med,:),[1 N_eta]),'k','linewidth',4);
set(gca,'box','on','fontname','times','fontsize',font_l)
title('\textbf{High Substitution}','fontname','times','fontsize',20)
ylabel('Exports, $T_x$','fontname','times','fontsize',font_l)
xlabel('Uncertainty','fontname','times','fontsize',font_l)
xlim([0 3])
ylim([1.33 1.41])

subplot(1,2,2)
plot(eta_vec,reshape(Export_eta_r_1(med,med,:),[1 N_eta]),'k','linewidth',4);
set(gca,'box','on','fontname','times','fontsize',font_l)
title('\textbf{Low Substitution}','fontname','times','fontsize',20)
ylabel('Exports, $T_x$','fontname','times','fontsize',font_l)
xlabel('Uncertainty','fontname','times','fontsize',font_l)
xlim([0 3])
ylim([1.358 1.385])

%Save figure
saveas(gcf,'fig4_twoelasticities.eps','epsc')
