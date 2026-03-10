function ydot = modelBasic(t,y,pars)

% parameters
a0 = pars(1);
a1 = pars(2); 
a2 = pars(3); 
mu = pars(4);
w1 = pars(5); 
a3 = pars(6); 
a4 = pars(7); 
w2 = pars(8);
a5 = pars(9); 
w3 = pars(10); 
d_e = pars(11);
d_n = pars(12);
a6 = pars(13);
k_e = pars(14);
alpha_e = pars(15);
k_n = pars(16);
alpha_n = pars(17);
alpha_s = pars(18);
k_s = pars(19);


% states
C = y(1);
A = y(2); 
B = y(3);
E = y(4);
N = y(5);

% DEs
dC = a0+circ(t-490)*a1/(1+a2*B^2)*(C/(mu+C))-w1*C;
dA = a3*C/(1+a4*B)-w2*A;
dB = a5*A-w3*B+a6*E;
dE = (alpha_e*B*N^2)/(k_e+B*N)-d_e*E; 
dN = (alpha_n*Stress(t))/(k_n+Stress(t))-d_n*N-(alpha_s*B*N)/(k_s+B*N);

% define circadian function
function C = circ(t)
Nc  = 0.5217;
eps = 0.01;
k = 5;
alpha = 300;
T = 1440;
l = 6;
delta = 83.8;
beta = 950;
C = Nc*((mod(t-delta,T))^k/((mod(t-delta,T))^k+alpha^k)*(T-(mod(t-delta,T)))^l/((T-(mod(t-delta,T)))^l+beta^l)+eps);
end 

% define stress function 
function S=Stress(t)
if 0<=t && t<83
    S = 31+90/(1+exp(-(t-60)/5));
else 
   S=121-90/(1+exp(-(t-92)/2)); 
end
end 

ydot = [dC; dA; dB; dE; dN];

end 