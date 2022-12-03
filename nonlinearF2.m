function [f2] = nonlinearF2(X1,X2)
    [r,F0,F3,Fr,F1,V1,V2,R,T03,T0,CA0s,CA03s,deltaH1,deltaH2,deltaH3,k10,k20,k30,E1,E2,E3,rho,cp,T1s,CA1s,T2s,CA2s] = systemParameters();
    %%%%%%%%%%%%%%%%%%%%%%% R's (reaction terms) %%%%%%%%%%%%%%%%%%%%%%%
    R1=k10*exp(-E1/R/T2s/(X2(1,1)+1))*CA2s*(X2(2,1)+1);
    R2=k20*exp(-E2/R/T2s/(X2(1,1)+1))*CA2s*(X2(2,1)+1);
    R3=k30*exp(-E3/R/T2s/(X2(1,1)+1))*CA2s*(X2(2,1)+1);
    %%%%%%%%%%%%%%%%%%%%%%% f1 %%%%%%%%%%%%%%%%%%%%%%%%
    f2(1,1)=F1/V2*((T1s*(1+X1(1,1))-T2s)/T2s-X2(1,1))+F3/V2*((T03-T2s)/T2s-X2(1,1))-deltaH1/rho/cp/T2s*R1-deltaH2/rho/cp/T2s*R2-deltaH3/rho/cp/T2s*R3;
    f2(2,1)=F1/V2*((CA1s*(1+X1(2,1))-CA2s)/CA2s-X2(2,1))+F3/V2*((CA03s-CA2s)/CA2s-X2(2,1))-(R1+R2+R3)/CA2s;
end