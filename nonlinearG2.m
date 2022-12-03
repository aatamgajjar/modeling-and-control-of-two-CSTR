function [g2] = nonlinearG2()
    [r,F0,F3,Fr,F1,V1,V2,R,T03,T0,CA0s,CA03s,deltaH1,deltaH2,deltaH3,k10,k20,k30,E1,E2,E3,rho,cp,T1s,CA1s,T2s,CA2s] = systemParameters();
    
    g2 = [1/rho/cp/V2/T2s 0;0 F3/V2/CA2s];
end