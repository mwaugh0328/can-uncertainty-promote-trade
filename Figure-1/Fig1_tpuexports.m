
%{
Can Global Uncertainty Promote International Trade?
Authors: Isaac Baley, Laura Veldkamp, and Mike Waugh 
Data: May 2020

---------------------------------------------------------------------------
---------------------------------------------------------------------------
             FIGURE 1: TRADE POLICY UNCERTAINTY AND EXPORTS
---------------------------------------------------------------------------
---------------------------------------------------------------------------

This code reads the Excel data "TPUandExports.xlsx" and plots Figure 1.  
%}

clear all;
close all;
clc

%% Import the data
[~, ~, raw] = xlsread('TPUandExports.xlsx','Sheet1');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw);                       % Find non-numeric cells
raw(R) = {NaN};                                                             % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Read data
exports_dates = data(:,1);
exports       = data(:,2);
tpu_dates     = data(:,3);
tpu           = data(:,4);


%% Creates Figure 1
axes1 = axes(figure);
hold(axes1,'on');

% Activate the left side of the axes
yyaxis(axes1,'left');
% Create plot
plot(tpu_dates,tpu,'DisplayName','Trade Policy Uncertainty (left axis)',...
    'LineWidth',3,...
    'Color',[0.6 0.6 0.6]);

% Create ylabel
ylabel('Trade Policy Uncertainty','FontName','Times New Roman');
ylim(axes1,[0 2200]);

% Set the remaining axes properties
set(axes1,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(axes1,'right');
% Create plot
plot(exports_dates,exports,'DisplayName','Exports (right axis)','LineWidth',3,...
    'Color',[0 0 0]);

% Create ylabel
ylabel('Exports','FontName','Times New Roman');
ylim(axes1,[90 110]);

% Set the remaining axes properties
set(axes1,'YColor',[0 0 0]);
xlim(axes1,[2014 2019]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',18,'LineStyleOrderIndex',...
    2);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','northwest','FontSize',20);
legend boxoff

%Save figure
saveas(gcf,'fig1_tpuexports.eps','epsc')



