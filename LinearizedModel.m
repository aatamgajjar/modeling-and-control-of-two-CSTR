%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% THIS IS USED TO GET THE LINEARIZED FORM OF    %%%%%%%%%%%%%%%
%%%%%%%%%%%%% THE MODELS OF THE TWO PROCESS                 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%This is doing something funky with the ODE system. . . I think it is
%%%%non dimentionalized but it may be done incorrectly. . . -James
function [A, B]=LinearizedModel(DELTA1,DELTA2,DELTA3)
%%%%%%%%%%%%%%%%%%%%%%% SYSTEM PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r,F0,F3,Fr,F1,V1,V2,R,T03,T0,CA0s,CA03s,deltaH1,deltaH2,deltaH3,k10,k20,k30,E1,E2,E3,rho,cp,T1s,CA1s,T2s,CA2s] = systemParameters;
%%% make the model imperfect %%%
% deltaH1_m=(1+DELTA1)*deltaH1;
% deltaH2_m=(1+DELTA2)*deltaH2;
% E1_m=(1+DELTA3)*E1;
% 
% deltaH1=deltaH1_m;
% deltaH2=deltaH2_m;
% % deltaH3=deltaH3_m;
% E1=E1_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% parameters of linearized model %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matrices describing models %
% first unit %
A1hat=zeros(2,2);
A12hat=zeros(2,2);
B1hat=zeros(2,2);
% second unit %
A2hat=zeros(2,2);
A21hat=zeros(2,2);
B2hat=zeros(2,2);

%%%%%%%%%%%%%%%% insert data to these matrices %%%%%%%%%%%%%%%%
% first unit %
%I think this is supposed to have a T1s term throughouth multiplying
%things. . . 
A1hat(1,1) = -(F0 + Fr)/V1 - deltaH1/rho/cp/T1s*k10*CA1s*exp(-E1/R/T1s)*E1/R/T1s - deltaH2/rho/cp/T1s*k20*CA1s*exp(-E2/R/T1s)*E2/R/T1s - deltaH3/rho/cp/T1s*k30*CA1s*exp(-E3/R/T1s)*E3/R/T1s;

%I believe this is not suppose to have a CA1s term but instead a T1s term?
%and I am not sure where this>|<  T1s terms comes from. . . 
A1hat(1,2) = -deltaH1/rho/cp/T1s*k10*exp(-E1/R/T1s)*CA1s - deltaH2/rho/cp/T1s*k20*exp(-E2/R/T1s)*CA1s - deltaH3/rho/cp/T1s*k30*exp(-E3/R/T1s)*CA1s;
A1hat(2,1) = -k10*exp(-E1/R/T1s)*E1/R/T1s-k20*exp(-E2/R/T1s)*E2/R/T1s-k30*exp(-E3/R/T1s)*E3/R/T1s;
A1hat(2,2) = -(F0+Fr)/V1 - k10*exp(-E1/R/T1s) - k20*exp(-E2/R/T1s) - k30*exp(-E3/R/T1s);

A12hat(1,1)=Fr*T2s/V1/T1s;
A12hat(2,2)=Fr*CA2s/V1/CA1s;

% B1hat(1,1)=1/rho/cp/V1/T1s; 
% B1hat(2,2)=F0/V1/CA1s;
% second unit %
A2hat(1,1)=-(F1+F3)/V2-deltaH1/rho/cp/T2s*k10*CA2s*exp(-E1/R/T2s)*E1/R/T2s-deltaH2/rho/cp/T2s*k20*CA2s*exp(-E2/R/T2s)*E2/R/T2s-deltaH3/rho/cp/T2s*k30*CA2s*exp(-E3/R/T2s)*E3/R/T2s;
A2hat(1,2)=-deltaH1/rho/cp/T2s*k10*exp(-E1/R/T2s)*CA2s-deltaH2/rho/cp/T2s*k20*exp(-E2/R/T2s)*CA2s-deltaH3/rho/cp/T2s*k30*exp(-E3/R/T2s)*CA2s;
A2hat(2,1)=-k10*exp(-E1/R/T2s)*E1/R/T2s-k20*exp(-E2/R/T2s)*E2/R/T2s-k30*exp(-E3/R/T2s)*E3/R/T2s;
A2hat(2,2)=-(F1+F3)/V2-k10*exp(-E1/R/T2s)-k20*exp(-E2/R/T2s)-k30*exp(-E3/R/T2s);

A21hat(1,1)=F1/V2*T1s/T2s;
A21hat(2,2)=F1/V2*CA1s/CA2s;

B1hat(1,1)=1/rho/cp/V1/T1s; 
B1hat(2,2)=F0/V1/CA1s;

B2hat(1,1)=1/rho/cp/V2/T2s;
B2hat(2,2)=F3/V2/CA2s;

A = [A1hat A12hat;
     A2hat A21hat];
B = [B1hat zeros(2,2);
     zeros(2,2),B2hat];
% B = zeros(4,4);
%     B(1,1)=1/rho/cp/V1/T1s; 
%     B(3,3)=1/rho/cp/V2/T2s;














