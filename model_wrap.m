%--------------------------------------------------------------------------
%Solve the ODE using the given parameters (pars) and initial conditions
%(Init) at the given points (ydata)
%Set ODE tolerance and data as global
%td and xdata
%--------------------------------------------------------------------------

function sol = model_sol(pars,data)
global ODE_TOL 

% load data, ICs
Adata = data.Adata;
Bdata = data.Bdata;
Edata = data.Edata; 
Ndata = data.Ndata;
xdata = data.xdata; 
Init  = data.Init;

% run solver
options = odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL);
sol   = ode45(@modelBasic,xdata,Init,options,pars);
sol   = deval(sol,xdata);

% get solution vectors
C = sol(1,:)';
A = sol(2,:)';
B = sol(3,:)';
E = sol(4,:)';
N = sol(5,:)';


