function find_jacobian
[r,F0,F3,Fr,F1,V1,V2,R,T03,T0,CA0s,CA03s,deltaH1,deltaH2,deltaH3,k10,k20,k30,E1,E2,E3,rho,cp,T1s,CA1s,T2s,CA2s] = systemParameters();
T1s = 300.38782;
CA1s = 2.49806;
T2s = 300.34962;
CA2s = 2.28400;

syms T1 CA1 T2 CA2
j = (jacobian([F0/V1*(T0 - T1) +  Fr/V1*(T2- T1) + (-deltaH1)/rho/cp*k10*exp(-E1/R/T1)*CA1 + (-deltaH2)/rho/cp*k20*exp(-E2/R/T1)*CA1  + (-deltaH3)/rho/cp*k30*exp(-E3/R/T1)*CA1;
    F0/V1*(CA0s - CA1) + Fr/V1*(CA2- CA1) - k10*exp(-E1/R/T1)*CA1- k20*exp(-E2/R/T1)*CA1- k30*exp(-E3/R/T1)*CA1;
    F1/V2*(T1- T2) + F3/V2*(T03-T2) + (-deltaH1)/rho/cp*k10*exp(-E1/R/T2)*CA2 + (-deltaH2)/rho/cp*k20*exp(-E2/R/T2)*CA2  + (-deltaH3)/rho/cp*k30*exp(-E3/R/T2)*CA2;
    F1/V2*(CA1 - CA2) + F3/V2*(CA03s- CA2) - k10*exp(-E1/R/T2)*CA2- k20*exp(-E2/R/T2)*CA2- k30*exp(-E3/R/T2)*CA2], [T1; CA1; T2;  CA2]));
js = subs(j, {T1, CA1, T2, CA2}, {T1s CA1s T1s CA2s})
end

