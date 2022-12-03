%To do
%Introduce disturbance Yulies work?
%XXX introduce ZOH
%Create line plots to figure out why H2 is better near unstable region
%Find a baseline h, tau, Phi (Stable and unstable), and alpha
%Introduce estimation? Check if possible?
%Create Contours for cases. . .
%Check H2 calc?



%Make plots similar to Sharou paper
%Make plots similar to Actuator delay paper
%Figure out the M^(k-1) in the performance calculation
%Add zero order hold
%Add no fault but disturbance
%Add phi estimation?
%%
%dt must be small here for things to work well
clear all
clc
close all

holdon = 0;
% Uncertainties
DELTA1=0; %deltaH1=(1+DELTA1)*deltaH1;
DELTA2=DELTA1; %deltaH2=(1+DELTA2)*deltaH2;
DELTA3=DELTA1; %deltaH3=(1+DELTA3)*deltaH3;
DELTA4=0; %E1=(1+DELTA4)*E1;

% Linearized models
%xdot=Ax + Bu + Ew
%y = Cx + Du

%XXXXXXX Need to fix symbolic issue XXXXXXX was unhappy with "syms T2s" So
%i change it do T2s=sym('T2s'); and that worked. . . once i repeated for
%all symbolic variables. . . 



%Initialization
[r,F0,F3,Fr,F1,...
    V1,V2,R,T03,T0,...
    CA0s,CA03s,...
    deltaH1,deltaH2,deltaH3,...
    k10,k20,k30,E1,E2,E3,...
    rho,cp,T1s,CA1s,T2s,CA2s] = systemParameters;

[Ahat,Bhat,Chat,Dhat,Kc] = YuleiSymsLinearization2CSTRImperfect;
stepdisturbance = 1;
X_offset = 5;
%introduce stable states (note that linearization give state space eqs in
%deviation form so X=zeros is at equilribum.
x1s = 1*T1s;
x2s = 1*CA1s;
x3s = 1*T2s;
x4s = 1*CA2s;
X_init = [413-T1s 1.94-CA1s 457-T2s 1.58-CA2s]';
% X_init = [0 0 0 0]';
tic  %timer starts
% A = 1000* [[0.0253, 1.2859, 0.035, 0]; [-0.0003, -0.0459, 0, 0.0350]; [0.0133, 0, -0.0028, 0.3357]; [0, 0.0133, -0.0001, -0.0249]];
% B = [[0.0043, 0, 0, 0]; [0, 4.998, 0, 0]; [0, 0, 0.0014, 0]; [0,0,0,10]];
A = [[25.3 4.97 31.75 0];[-78.03 -45.94 0 34.64];[-2.84 1.42 14.7 0];[-22.45 -24.88 0 13.47]];
B = [[0.945*10^-5 0 0 0];[0 2.82 0 0];[0 0 0.347*10^-5 0];[0 0 0 5.71]];
% A = [[-39.7777, 1.3111, 34.998, 0]; [-0.001, -40.0021, 0, 34.998]; [13.332, 0, -23.1329, 1.3078]; [0, 13.3320, -0.0009, -23.3380]];
% B = [[0.0043, 0, 0, 0]; [0, 4.998, 0, 0]; [0, 0, 0.0014, 0]; [0,0,0,10]];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulation Parameters  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time unit: min
t0 = 0;
dt = 0.001; %if this is 0.001 then the simulation does not return the expected results with regards to stability of the system
tn = 10;
t = t0:dt:tn;
Num = numel(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Accom plot
%0.884- This is on the limit of stability for h=0.1 tau=0.04 accom = 10
%0.999 This is on the limit of stability for h=0.4 tau=0.17+ accom = 10
%Very interesting simulation ^^^^^
%No accommodation plot:
%0.917- This is on the limit of stability for h=0.05 tau=0
%0.958- This is on the limit of stability for h=0.1 tau=0.05 accom = 10
%0.953- This is on the limit of stability for h=0.1 tau=0.04 accom = 10
%This simulation is giving some trouble ^^^ is it so close to the stability
%limit that it is taking forever or maybe try a smaller dt?


X = zeros(4,Num);
UD = zeros(4, Num);
X(:,1) = X_init;
linearsys = ss(A, B, eye(4), zeros(4));
p = [-0.5, -1, -2, -5];

for pole = 1:numel(p)
    lambda = [p(pole) p(pole) p(pole) p(pole)];
    K = place(A, B, lambda);
    
    if stepdisturbance == 0
        UD(:,:) = 0;
    else
        for i = 600:900
            UD(1,i) = 5*25000;
            UD(2,i) = 0.1;
            UD(3,i) = 5*20000;
            UD(4,i) = 0.1;
        end
        for i = 2100:2400
            UD(1,i) = 3*10000;
            UD(2,i) = 0.05;
            UD(3,i) = 3*10000;
            UD(4,i) = 0.05;
        end
        for i = 4200:4500
            UD(1,i) = 1*15000;
            UD(2,i) = 0.025;
            UD(3,i) = 1*15000;
            UD(4,i) = 0.025;
        end
    end
    for k = 1:Num-1
        X(:,k+1)        = X(:,k)     +dt*((A-B*K)*X(:,k) + B*UD(:,k));
    end

    T1 = zeros(1,Num);
    Ca1 = zeros(1,Num); 
    T2 = zeros(1,Num); 
    Ca2 = zeros(1,Num);
    T1(1,:) = X(1,:) + T1s; %*T1s + T1s;
    Ca1(1,:) = X(2,:) + CA1s; %*CA1s + CA1s;
    T2(1,:) = X(3,:) + T2s;% *T2s + T2s;
    Ca2(1,:) = X(4,:) + CA2s;% *CA2s + CA2s;
    clsys = ss(A-B*K, B, eye(4), zeros(4));
%     stepinfo(clsys(1,1))
%     stepinfo(clsys(2,2))
%     stepinfo(clsys(3,3))
%     stepinfo(clsys(4,4))
    % T2nl(1,:) = T2(1,:);
    % for i = 3000:5000
    %     T2nl(1,i) = T2(1,i) + 1;
    % end
    % for i = 7500:Num
    %     T2nl(1,i) = T2(1,i) + 0.5;
    % end
    toc  % timer ends
    
    NormalColor = [0/255 173/255 199/255]; %Blue
    DriftColor = [227/255 83/255 93/255]; %Red
    IdentificationColor = [71/255 171/255 108/255]; %Green
    UnitLineColor = [69/255 74/255 73/255]; %Dark Grey
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    Color = [[0/255 173/255 199/255]; [227/255 83/255 93/255]; [71/255 171/255 108/255]; [69/255 74/255 73/255]];
  
    NormalLinestyle = '-';
    
    
    DriftLinestyle = '-';
    IdentificationLinestyle = '-';
    LineWidth = 2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% Generate State Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Temperature 1
    figure(1)
    subplot(2,2,1);
    plot(t(1:Num),T1(1,1:Num),'Linestyle','-','color',Color(pole,:))
    
    hold on
    
    % plot(t(1:Num),T1nl(1,1:Num),'Linestyle',NormalLinestyle,'color',NormalColor)
    xlabel('Time [hr]','FontSize',20)
    ylabel('T_1 [K]','FontSize',20)
    set(gca,'FontSize',20);
    set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
    % legend('Linear','Nonlinear')
    subplot(2,2,2);
    % Concentration 1
    plot(t(1:Num),Ca1(1,1:Num),'Linestyle','-','color',Color(pole,:))
    hold on
    % plot(t(1:Num),Ca1nl(1,1:Num),'Linestyle',NormalLinestyle,'color',NormalColor)
    xlabel('Time [hr]','FontSize',20)
    ylabel('C_{A1}[kmol/m^{3}]','FontSize',20)
    set(gca,'FontSize',20);
    set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
    % ylim([2,4])
    % legend('Linear','Nonlinear')
    % Temperature 2
    subplot(2,2,3);
    plot(t(1:Num),T2(1,1:Num),'Linestyle','-','color',Color(pole,:))
    hold on
    % plot(t(1:Num),T2nl(1,1:Num),'Linestyle',NormalLinestyle,'color',NormalColor)
    xlabel('Time [hr]','FontSize',20)
    ylabel('T_2 [K]','FontSize',20)
    set(gca,'FontSize',20);
    set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
    % legend('Linear','Nonlinear')
    % Concentration 2
    subplot(2,2,4);
    plot(t(1:Num),Ca2(1,1:Num),'Linestyle','-','color',Color(pole,:))
    hold on
    % plot(t(1:Num),Ca2nl(1,1:Num),'Linestyle',NormalLinestyle,'color',NormalColor)
    xlabel('Time [hr]','FontSize',20)
    ylabel('C_{A2} [kmol/m^{3}]','FontSize',20)
    set(gca,'FontSize',20);
    set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
    % legend('Linear','Nonlinear')
end
