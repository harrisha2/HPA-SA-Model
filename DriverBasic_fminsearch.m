function [ss,sol,rout] = model_fmin(pars)

global data

% load data, ICs
Adata = data.Adata;
Bdata = data.Bdata;
Edata = data.Edata; 
Ndata = data.Ndata;
xdata = data.xdata; 
init  = data.Init; 

% get solution vectors
sol = model_wrap(pars);
A = sol(2,:)';
B = sol(3,:)';
E = sol(4,:)';
N = sol(5,:)';

% calculate residuals normalized by mean of observations
Arout = (A - Adata')/mean(Adata);
Brout = (B - Bdata')/mean(Bdata);
Erout = (E - Edata')/mean(Edata);
Nrout = (N - Ndata')/mean(Ndata);

% concatenate vectors of residuals for variables of interest
rout = [Brout; Erout; Nrout];

% calculate SSR 
ss = rout'*rout;

