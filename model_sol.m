function [x,histout,jachist,rout,sc] = DriverBasic_fminsearch
close all
clear all

global ALLPARS INDMAP data

% Load data
data = load('ABEN.mat');
Adata = data.Adata;
Bdata = data.Bdata;
Edata = data.Edata; 
Ndata = data.Ndata;

% SPECIFY DATASET AS sub_id IN LINE 45 
% 1 for A, 2 for B, 3 for C, 4 for D, 5 for E
sub_id = 5;

%  select data for subject 
 Adata = Adata(sub_id,:);
 Bdata = Bdata(sub_id,:);
 Edata = Edata(sub_id,:);
 Ndata = Ndata(sub_id,:);

xdata = data.time_data;
% Indices for parameters estimated (sensitive parameters)
INDMAP = [10 11 12 13 15 16 17];

% Get nominal parameter values, upper and lower bound for optimization
[pars,Init,lo,hi] = load_global(Adata,Bdata,Ndata,Edata, INDMAP, sub_id);

% Create structure with data and initial conditions
data.Adata = Adata;
data.Bdata = Bdata;
data.Edata = Edata;
data.Ndata = Ndata;
data.xdata = xdata;
data.Init = Init;

%--------  Solution w/ Nominal parameter values -----------------------
% Integrate the system using ode45.m
sol = model_sol(pars,data);

Csol = sol(1,:)';
Asol = sol(2,:)';
Bsol = sol(3,:)';
Esol = sol(4,:)';
Nsol = sol(5,:)';


%----------------------- Nelder-Mead (fminsearch) -----------------------

ALLPARS = pars;

opts = optimset('MaxIter',500,'MaxFunEvals',7500,'TolX',1e-10);
[xopt, fval, exitflag, output] = fmincon(@model_fmin,pars(INDMAP),[],[],[],[],lo,hi,[],opts);

fval
exitflag
output

parsNM = ALLPARS;
% print final parameter set (fixed and optimized) 
parsNM(INDMAP) = xopt

% Integrate the system using ode45.m
sol = model_sol(parsNM,data);
Copt = sol(1,:)';
Aopt = sol(2,:)';
Bopt = sol(3,:)';
Eopt = sol(4,:)';
Nopt = sol(5,:)';

%% combined plots (used in manuscript)
figure(4); %hold on
set(gca, 'FontSize',12,'FontName','Arial');
t = tiledlayout(1,3);
set(gcf,'units','centimeters','position',[0,0,60,10])

% plot cortisol 
nexttile;
set(gca, 'FontSize',12,'FontName','Arial');
pp = plot(xdata-50, Bdata,'ob', xdata-50,Bopt,'-r');
hold on;
set(pp,'linewidth',2);
ylabel('Cortisol (\mug/dL)', 'FontSize',12)
xregion(0,37);

% plot epinephrine
nexttile;
set(gca, 'FontSize',12,'FontName','Arial');
qq = plot(xdata-50, Edata,'ob', xdata-50,Eopt,'-r');
set(qq,'linewidth',2);
ylabel('Epinephrine (\mug/dL)', 'FontSize',12)
xregion(0,37);

% plot norepinephrine
nexttile;
set(gca, 'FontSize',12,'FontName','Arial');
mm = plot(xdata-50,Ndata,'ob', xdata-50,Nopt,'-r' );
set(mm,'linewidth',2);
ylabel('Norepinephrine (\mug/dL)', 'FontSize',12)
xregion(0,37);

legend( 'Data', 'Optimized');
lgd = legend;
lgd.Layout.Tile = 'east';
% title plot by sub_id (subject 6 is missing from data set) 
if sub_id<6
    title(t,['Subject ', num2str(sub_id)],'FontSize',16);
else 
    title(t,['Subject ', num2str(sub_id+1)],'FontSize',16);  
end
xlabel(t,'Time (min)','FontSize',12)



