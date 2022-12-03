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

X_offset = 5;
%introduce stable states (note that linearization give state space eqs in
%deviation form so X=zeros is at equilribum.
x1s = 1*T1s;
x2s = 1*CA1s;
x3s = 1*T2s;
x4s = 1*CA2s;
X_init = [0 0 0 0]';

tic  %timer starts


A = [[-39.7777, 1.3111, 34.998, 0]; [-0.001, -40.0021, 0, 34.998]; [13.332, 0, -23.1329, 1.3078]; [0, 13.3320, -0.0009, -23.3380]];
B = [[0.0043, 0, 0, 0]; [0, 4.998, 0, 0]; [0, 0, 0.0014, 0]; [0,0,0,10]];

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

X(:,1) = X_init;

U = zeros(4,Num);

% 
% for i = 3000:5000
%     U(1,i) = 2*25000;
%     U(2,i) = 1;
%     U(3,i) = 2*20000;
%     U(4,i) = 1;
% end
% for i = 5000:7500
%     U(1,i) = 2*10000;
%     U(2,i) = 0.5;
%     U(3,i) = 2*10000;
%     U(4,i) = 0.5;
% end
% for i = 7500:Num
%     U(1,i) = 2*15000;
%     U(2,i) = 0.75;
%     U(3,i) = 2*15000;
%     U(4,i) = 0.75;
% end

for i = 3000:5000
    U(1,i) = 1*25000;
    U(2,i) = 10;
    U(3,i) = 1*20000;
    U(4,i) = 10;
end
for i = 5000:7500
    U(1,i) = 1*10000;
    U(2,i) = 15;
    U(3,i) = 1*10000;
    U(4,i) = 15;
end
for i = 7500:Num
    U(1,i) = 1*15000;
    U(2,i) = 7.5;
    U(3,i) = 1*15000;
    U(4,i) = 7.5;
end
% U(2,:) = CA0s;
% U(4,:) = CA03s;




for k = 1:Num-1

    X(:,k+1)        = X(:,k)     +dt*(A*X(:,k)       +B*U(:,k));

end

X1nl = zeros(2,Num);
X2nl = zeros(2,Num);

U1(1,:) = U(1,:);
U1(2,:) = U(2,:);
U2(1,:) = U(3,:);
U2(2,:) = U(4,:);

f1 = zeros(2,Num);
f2 = zeros(2,Num);

T1nl = zeros(1,Num);
Ca1nl = zeros(1,Num);
T2nl = zeros(1,Num);
Ca2nl = zeros(1,Num);

f1(:,1) = nonlinearF1(X1nl(:,1),X2nl(:,1));
f2(:,1) = nonlinearF2(X1nl(:,1),X2nl(:,1));
g1 = nonlinearG1();
g2 = nonlinearG2();

for  k = 2:Num
    X1nl(:,k) = X1nl(:,k-1) + dt*(f1(:,k-1) + g1*U1(:,k-1));
    X2nl(:,k) = X2nl(:,k-1) + dt*(f2(:,k-1) + g2*U2(:,k-1));
    f1(:,k) = nonlinearF1(X1nl(:,k),X2nl(:,k));
    f2(:,k) = nonlinearF2(X1nl(:,k),X2nl(:,k));
%     T1nl(1,k) = X1nl(1,k)*T1s + T1s;
%     Ca1nl(1,k) = X1nl(2,k) + CA1s;
%     T2nl(1,k) = X2nl(1,k) + T2s;
%     Ca2nl(1,k) = X2nl(2,k) + CA2s;
end
T1nl(1,:) = X1nl(1,:)*T1s + T1s;
Ca1nl(1,:) = X1nl(2,:) + CA1s;
T2nl(1,:) = X2nl(1,:) + T2s;
Ca2nl(1,:) = X2nl(2,:) + CA2s;


T1 = zeros(1,Num);
Ca1 = zeros(1,Num); 
T2 = zeros(1,Num); 
Ca2 = zeros(1,Num);
T1(1,:) = X(1,:) + T1s; %*T1s + T1s;
Ca1(1,:) = X(2,:)+ CA1s; %*CA1s + CA1s;
T2(1,:) = X(3,:) + T2s;% *T2s + T2s;
Ca2(1,:) = X(4,:) + CA2s;% *CA2s + CA2s;

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




NormalLinestyle = '-';


DriftLinestyle = '-';
IdentificationLinestyle = '-';
LineWidth = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Generate State Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature 1
figure(1)
subplot(2,1,2);
plot(t(1:Num),T1(1,1:Num),'Linestyle','--','color',DriftColor)

hold on

plot(t(1:Num),T1nl(1,1:Num),'Linestyle',NormalLinestyle,'color',NormalColor)
xlabel('Time [hr]','FontSize',20)
ylabel('T_1 [K]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
legend('Linear','Nonlinear')
subplot(2,1,1);
plot(t(1:Num),U(1,1:Num),'Linestyle',NormalLinestyle,'color',UnitLineColor)
xlabel('Time [hr]','FontSize',20)
ylabel('Q1 [J]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)

% Concentration 1
figure(2)
subplot(2,1,2);
plot(t(1:Num),Ca1(1,1:Num),'Linestyle','--','color',DriftColor)
hold on
plot(t(1:Num),Ca1nl(1,1:Num),'Linestyle',NormalLinestyle,'color',NormalColor)
xlabel('Time [hr]','FontSize',20)
ylabel('C_{A1}[kmol/m^{3}]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
ylim([2,4])
legend('Linear','Nonlinear')
subplot(2,1,1);
plot(t(1:Num),U(2,1:Num),'Linestyle',NormalLinestyle,'color',UnitLineColor)
xlabel('Time [hr]','FontSize',20)
ylabel('C_{A0} [kmol/m^{3}]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)


% Temperature 2
figure(3)
subplot(2,1,2);
plot(t(1:Num),T2(1,1:Num),'Linestyle','--','color',DriftColor)
hold on
plot(t(1:Num),T2nl(1,1:Num),'Linestyle',NormalLinestyle,'color',NormalColor)
xlabel('Time [hr]','FontSize',20)
ylabel('T_2 [K]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
legend('Linear','Nonlinear')
subplot(2,1,1);
plot(t(1:Num),U(3,1:Num),'Linestyle',NormalLinestyle,'color',UnitLineColor)
xlabel('Time [hr]','FontSize',20)
ylabel('Q2 [J]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)

% Concentration 2
figure(4)
subplot(2,1,2);
plot(t(1:Num),Ca2(1,1:Num),'Linestyle','--','color',DriftColor)
hold on
plot(t(1:Num),Ca2nl(1,1:Num),'Linestyle',NormalLinestyle,'color',NormalColor)
xlabel('Time [hr]','FontSize',20)
ylabel('C_{A2} [kmol/m^{3}]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
legend('Linear','Nonlinear')

ylim([2,4.5])
subplot(2,1,1);
plot(t(1:Num),U(4,1:Num),'Linestyle',NormalLinestyle,'color',UnitLineColor)
xlabel('Time [hr]','FontSize',20)
ylabel('C_{A3} [kmol/m^{3}]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
