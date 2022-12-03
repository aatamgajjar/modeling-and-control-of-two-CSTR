function [g1] = nonlinearG1()
    [r,F0,F3,Fr,F1,V1,V2,R,T03,T0,CA0s,CA03s,deltaH1,deltaH2,deltaH3,k10,k20,k30,E1,E2,E3,rho,cp,T1s,CA1s,T2s,CA2s] = systemParameters();
    
    g1 = [1/rho/cp/V1/T1s 0;0 F0/V1/CA1s];
end