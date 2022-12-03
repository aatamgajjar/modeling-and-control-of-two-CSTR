function [f1] = nonlinearF1(X1,X2)
    [r,F0,F3,Fr,F1,V1,V2,R,T03,T0,CA0s,CA03s,deltaH1,deltaH2,deltaH3,k10,k20,k30,E1,E2,E3,rho,cp,T1s,CA1s,T2s,CA2s] = systemParameters();
    
    %%%%%%%%%%%%%%%%%%%%%%% R's (reaction terms) %%%%%%%%%%%%%%%%%%%%%%%
    R1=k10*exp(-E1/R/T1s/(X1(1,1)+1))*CA1s*(X1(2,1)+1);
    R2=k20*exp(-E2/R/T1s/(X1(1,1)+1))*CA1s*(X1(2,1)+1);
    R3=k30*exp(-E3/R/T1s/(X1(1,1)+1))*CA1s*(X1(2,1)+1);
    %%%%%%%%%%%%%%%%%%%%%%% f1 %%%%%%%%%%%%%%%%%%%%%%%%
    f1(1,1)=F0/V1*((T0-T1s)/T1s-X1(1,1))+Fr/V1*((T2s*(1+X2(1,1))-T1s)/T1s-X1(1,1))-deltaH1/rho/cp/T1s*R1-deltaH2/rho/cp/T1s*R2-deltaH3/rho/cp/T1s*R3;
    f1(2,1)=F0/V1*((CA0s-CA1s)/CA1s-X1(2,1))+Fr/V1*((CA2s*(1+X2(2,1))-CA1s)/CA1s-X1(2,1))-(R1+R2+R3)/CA1s;
end