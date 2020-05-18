%{
Can Global Uncertainty Promote International Trade?
Authors: Isaac Baley, Laura Veldkamp, and Mike Waugh 
Data: May 2020

---------------------------------------------------------------------------
---------------------------------------------------------------------------
                         FIGURES 2 and 3 

Reads the simulation results  "results_simu.mat" and plots Figures 2 and 3. 
Computes trade correlation and terms of trade moments under the two info
scenarios.
---------------------------------------------------------------------------
---------------------------------------------------------------------------

%}


load('results_twoprecisions.mat')

%% FIGURE 2:  Export Coordination under Perfect and Imperfect Information


font_l = 20;
font_s = 18;
set(0,'DefaultTextinterpreter','latex')

%pick periods (for illustration, we pick 20 periods that deliver terms of trade spikes)
tmin = 25;
tmax = 44;

figure(1)
subplot(2,1,1)
plot(log(T_x_simu(tmin:tmax,1)),'k--','linewidth',3)
hold on
plot(log(T_y_simu(tmin:tmax,1)),'color',[0,0,0]+0.5,'linewidth',3)
title('\textbf{Perfect Information}','fontname','times','fontsize',font_l)
set(gca,'box','on','fontname','times','fontsize',font_s)
ylabel('Log Exports','fontname','times','fontsize',font_l)
xlabel('Periods','fontname','times','fontsize',font_s)
legend('Country X', 'Country Y')
legend boxoff
text(1,1.2,'$corr[T_x,T_y]$ = 0.9 ','fontname','times','Fontsize',font_s)
ylim([-1.5 1.5])

subplot(2,1,2)
plot(log(T_x_simu(tmin:tmax,2)),'k--','linewidth',3)
hold on
plot(log(T_y_simu(tmin:tmax,2)),'color',[0,0,0]+0.5,'linewidth',3)
title('\textbf{Imperfect Information}','fontname','times','fontsize',font_l)
set(gca,'box','on','fontname','times','fontsize',font_s)
ylabel('Log Exports','fontname','times','fontsize',font_l)
xlabel('Periods','fontname','times','fontsize',font_s)
legend('Country X', 'Country Y')
legend boxoff
text(1,1.2,'$corr[T_x,T_y] = 0$ ','fontname','times','Fontsize',font_s)
ylim([-1.5 1.5])

saveas(gcf,'fig2_simutradelevel.eps','epsc')

%% FIGURE 3: Uncertainty Increases Terms of Trade Volatility and Average


figure(2)
y1=plot(q_simu(tmin:tmax,1),'k','linewidth',3);
hold on
y2=plot(q_simu(tmin:tmax,2),'color',[0,0,0]+0.5,'linewidth',3);
hold on
plot(mean(q_simu(tmin:tmax,1))*ones(20,1),'k--','linewidth',3)
hold on
plot(mean(q_simu(tmin:tmax,2))*ones(20,1),'--','color',[0,0,0]+0.5,'linewidth',3)
set(gca,'box','on','fontname','times','fontsize',font_s)
ylabel('Terms of Trade','fontname','times','fontsize',font_l)
xlabel('Periods','fontname','times','fontsize',font_s)
h_leg=legend([y1 y2], 'Perfect Information','Imperfect Information');
set(h_leg,'box','on','location','Northeast','fontsize',font_l);
legend('Perfect Information', 'Imperfect Information')
legend boxoff
text(5,4,'Average Terms of Trade ','fontname','times','Fontsize',font_s)
ylim([0 8.5])
% Create arrow
annotation('arrow',[0.393767705382436 0.463172804532578],...
    [0.456831325301205 0.238955823293173]);

% Create arrow
annotation('arrow',[0.405099150141643 0.531161473087819],...
    [0.462855421686747 0.311244979919679]);

saveas(gcf,'fig3_simutermstrade.eps','epsc')

%% STATISTICS: %rade correlation and terms of trade moments under the two info scenarios

disp('trade correlation perfect info')
corr(log(T_x_simu(:,1)),log(T_y_simu(:,1)))
disp('trade correlation, imperfect info')
corr(log(T_x_simu(:,2)),log(T_y_simu(:,2)))
disp('terms of trade moments E, V perfect')
M_imp = [ mean(q_simu(:,1)), var(q_simu(:,1))]
disp('terms of trade moments E,V imperfect')
M_per = [ mean(q_simu(:,2)), var(q_simu(:,2))]

