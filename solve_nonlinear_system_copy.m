syms T1 CA1 T2 CA2 Q1 CA0 Q2 CA3 

% Parameters
F0 = 4.998; %m3/h
F1 = 39.996; %m3/h
F3 = 30.0; %m3/h
Fr = 34.998; %m3/h
V1 = 1.0; %m3
V2 = 3.0; %m3
R = 8.314; %kJ/kmol K
T0 = 300.0; %K
T3 = 300.0; %K
CA0s = 4.0; %kmol/m3
CA1s = 1.77; %kmol/m3
CA2s = 1.75; %kmol/m3
CA3s = 2.0; %kmol/m3
delH1 = -5.0*10^4; %kJ/kmol
delH2 = -5.2*10^4; %kJ/kmol
delH3 = -5.4*10^4; %kJ/kmol 
k1 = 3.0*10^6; %h
k2 = 3.0*10^5; %h
k3 = 3.0*10^5; %h 
E1 = 5.0*10^4; %kJ/kmol
E2 = 7.53*10^4; %kJ/kmol
E3 = 7.53*10^4; %kJ/mol
rho = 1000.0; %kg/m3
Cp = 0.231; %kJ/kg K
T1s = 457.9; %K
T2s = 415.5; %K

% Unstable steady state point
T1_uss = 457.9; 
CA1_uss = 1.77; 
T2_uss = 415.5;
CA2_uss = 1.75;

f1 = (F0/V1)*(T0 - T1) + (Fr/V1)*(T2 - T1) - (CA1/rho/Cp)*( delH1*k1*exp(-E1/R/T1) + delH2*k2*exp(-E2/R/T1) + delH3*k3*exp(-E3/R/T1) ) + Q1/(rho*Cp*V1);
f2 = (F0/V1)*(CA0 - CA1) + (Fr/V1)*(CA2- CA1) - CA1*( k1*exp(-E1/R/T1) + k2*exp(-E2/R/T1) + k3*exp(-E3/R/T1) ) ;
f3 = (F1/V2)*(T1 - T2) + (F3/V2)*(T3 - T2) - (CA2/rho/Cp)*( delH1*k1*exp(-E1/R/T2) + delH2*k2*exp(-E2/R/T2) + delH3*k3*exp(-E3/R/T2) ) + Q2/(rho*Cp*V2);
f4 = (F1/V2)*(CA1 - CA2) + (F3/V2)*(CA3 - CA2) - CA2*( k1*exp(-E1/R/T2) + k2*exp(-E2/R/T2) + k3*exp(-E3/R/T2) );

F = [f1; f2; f3; f4];

% A matrix in symbolic form 
A_sym = jacobian(F, [T1, CA1, T2, CA2]);

% B matrix in symbolic form
B_sym = jacobian(F, [Q1, CA0, Q2, CA3]);

% Find roots (equilibrium points) of a multi-variable system 
% At steady state, Q1 = 0, Q2 =  Q1, CA0 = CA0s, CA3 = CA3s;
% Refer to root4d.m. To find one of the asymptoticall stable roots, use x0
% = [T1_uss, CA1_uss, T2_uss, CA2_uss]
% In the command window, type fun = @root4d THEN x0 THEN x=fsolve(fun,x0)
% To check value of f at roots x, do a = [fun(x)] and check each element
% s.t. a(1,1), a(1,2), a(1,3), a(1,4)

% Unstable s-s pt is (457.9428, 1.7702, 415.4585, 1.7522)
% One asymptotically stable eq pt is (300.3878, 2.4981, 300.3496, 2.2840) -
% Linearize about the stable eq pt (I had a sign error in my F initially
% which is why I got something else initially)

% A matrix at unstable steady state
A_J = subs(A_sym, {T1, CA1, T2, CA2}, {457.9, 1.77, 415.5, 1.75})

% B matrix at unstable steady state
% At steady state, Q1 = 0, CA0 = CA0s, Q2 =  Q1, CA3 = CA3s;
B_J = subs(B_sym, {Q1, CA0, Q2, CA3}, {0, CA0s, 0, CA3s})

% A matrix at asym stable eq pt 
% A_J = subs(A_sym, {T1, CA1, T2, CA2}, {300.3878, 2.4981, 300.3496, 2.2840})

% % B matrix at steady state
% % At steady state, Q1 = 0, CA0 = CA0s, Q2 =  Q1, CA3 = CA3s;
% B_J = subs(B_sym, {Q1, CA0, Q2, CA3}, {0, CA0s, 0, CA3s})


% double() gives numerical appx to A_J and B_J matrices!

