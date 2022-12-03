clear all
clc
clf

[r,F0,F3,Fr,F1,...
    V1,V2,R,T03,T0,...
    CA0s,CA03s,...
    deltaH1,deltaH2,deltaH3,...
    k10,k20,k30,E1,E2,E3,...
    rho,cp,T1s,CA1s,T2s,CA2s] = systemParameters;


A = [[25.3 4.97 31.75 0];[-78.03 -45.94 0 34.64];[-2.84 1.42 14.7 0];[-22.45 -24.88 0 13.47]];
B = [[0.945*10^-5 0 0 0];[0 2.82 0 0];[0 0 0.347*10^-5 0];[0 0 0 5.71]];

C = eye(4);

D = zeros(4);

Q_states = [[5 0 0 0];
            [0 10 0 0];
            [0 0 5 0];
            [0 0 0 10]];

R_inputs = [[1 0 0 0];
            [0 2 0 0];
            [0 0 1 0];
            [0 0 0 2]];

[K,S,P] = lqr(A,B,Q_states,R_inputs)

p = [-2 -5 -2 -5];
Kpp = place(A, B, p);
% Closed-loop step response with generated gain matrix K
x1s = 1*T1s;
x2s = 1*CA1s;
x3s = 1*T2s;
x4s = 1*CA2s;
X_init = [413-T1s 1.94-CA1s 457-T2s 1.58-CA2s]';

sys = ss(A-B*K,B,C,D);
dt = 0.001;
t = 0:dt:10;
Num = numel(t);
UD = zeros(4, Num);
X = zeros(4,Num);
Xpp = zeros(4, Num);
for i = 600:900
    UD(1,i) = 1*25000;
    UD(2,i) = 0.1;
    UD(3,i) = 1*20000;
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
% Xpp(:,1) = X_init;

for k = 1:Num-1
    X(:,k+1)        = X(:,k)     +dt*((A - B*K)*X(:,k)  + B*UD(:,k));
    Xpp(:,k+1) = Xpp(:,k)     +dt*((A - B*Kpp)*Xpp(:,k)  + B*UD(:,k));
end
U(:,:) = -K*X(:,:) + UD(:,:);
Upp(:, :) = -Kpp*Xpp(:,:) + UD(:,:);
Q1(:,:) = U(1,:);
Q1pp(:,:) = Upp(1,:);
Q2(:,:) = U(3,:);
Q2pp(:,:) = Upp(3,:);
Ca0(:,:) = U(2,:) + CA0s;
Ca0pp(:,:) = Upp(2,:) + CA0s;
Ca3(:,:) = U(4,:) + CA03s;
Ca3pp(:,:) = Upp(4,:) + CA03s;
T1 = zeros(1,Num);
Ca1 = zeros(1,Num); 
T2 = zeros(1,Num); 
Ca2 = zeros(1,Num);
T1(1,:) = X(1,:) + T1s; %*T1s + T1s;
Ca1(1,:) = X(2,:) + CA1s; %*CA1s + CA1s;
T2(1,:) = X(3,:) + T2s;% *T2s + T2s;
Ca2(1,:) = X(4,:) + CA2s;% *CA2s + CA2s;
T1pp(1,:) = Xpp(1,:) + T1s; %*T1s + T1s;
Ca1pp(1,:) = Xpp(2,:) + CA1s; %*CA1s + CA1s;
T2pp(1,:) = Xpp(3,:) + T2s;% *T2s + T2s;
Ca2pp(1,:) = Xpp(4,:) + CA2s;% *CA2s + CA2s;

NormalColor = [0/255 173/255 199/255]; %Blue
DriftColor = [227/255 83/255 93/255]; %Red
IdentificationColor = [71/255 171/255 108/255]; %Green
UnitLineColor = [69/255 74/255 73/255]; %Dark Grey
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Color = [[0/255 173/255 199/255]; [227/255 83/255 93/255]; [71/255 171/255 108/255]; [69/255 74/255 73/255]];

NormalLinestyle = '-';
LineWidth = 2;

figure(1)
subplot(4,2,1);
plot(t(1:Num),Q1(1,1:Num),'Linestyle','-','color',UnitLineColor)
hold on
plot(t(1:Num),Q1pp(1,1:Num),'Linestyle','--','color',UnitLineColor)
xlabel('Time [hr]','FontSize',20)
ylabel('Q_1 [kJ]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
legend('LQR', 'Pole Placement')
subplot(4,2,2);
plot(t(1:Num),Ca0(1,1:Num),'Linestyle','-','color',UnitLineColor)
hold on
plot(t(1:Num),Ca0pp(1,1:Num),'Linestyle','--','color',UnitLineColor)
xlabel('Time [hr]','FontSize',20)
ylabel('C_{A0} [kJ]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
legend('LQR', 'Pole Placement')

subplot(4,2,3);
plot(t(1:Num),T1(1,1:Num),'Linestyle','-','color',NormalColor)
hold on
plot(t(1:Num),T1pp(1,1:Num),'Linestyle','-','color',DriftColor)
% plot(t(1:Num),T1nl(1,1:Num),'Linestyle',NormalLinestyle,'color',NormalColor)
xlabel('Time [hr]','FontSize',20)
ylabel('T_1 [K]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
legend('LQR', 'Pole Placement')
subplot(4,2,4);
% Concentration 1
plot(t(1:Num),Ca1(1,1:Num),'Linestyle','-','color',NormalColor)
hold on
plot(t(1:Num),Ca1pp(1,1:Num),'Linestyle','-','color',DriftColor)
xlabel('Time [hr]','FontSize',20)
ylabel('C_{A1}[kmol/m^{3}]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
legend('LQR', 'Pole Placement')
% legend('Linear','Nonlinear')
% Temperature 2

subplot(4,2,5);
plot(t(1:Num),Q2(1,1:Num),'Linestyle','-','color',UnitLineColor)
hold on
plot(t(1:Num),Q2pp(1,1:Num),'Linestyle','--','color',UnitLineColor)
xlabel('Time [hr]','FontSize',20)
ylabel('Q_2 [kJ]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
legend('LQR', 'Pole Placement')

subplot(4,2,6);
plot(t(1:Num),Ca3(1,1:Num),'Linestyle','-','color',UnitLineColor)
hold on
plot(t(1:Num),Ca3pp(1,1:Num),'Linestyle','--','color',UnitLineColor)
xlabel('Time [hr]','FontSize',20)
ylabel('C_{A3} [kJ]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
legend('LQR', 'Pole Placement')

subplot(4,2,7);
plot(t(1:Num),T2(1,1:Num),'Linestyle','-','color',NormalColor)
hold on
plot(t(1:Num),T2pp(1,1:Num),'Linestyle','-','color',DriftColor)
xlabel('Time [hr]','FontSize',20)
ylabel('T_2 [K]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
legend('LQR', 'Pole Placement')
% Concentration 2
subplot(4,2,8);
plot(t(1:Num),Ca2(1,1:Num),'Linestyle','-','color',NormalColor)
hold on
plot(t(1:Num),Ca2pp(1,1:Num),'Linestyle','-','color',DriftColor)
xlabel('Time [hr]','FontSize',20)
ylabel('C_{A2} [kmol/m^{3}]','FontSize',20)
set(gca,'FontSize',20);
set(findobj(gca, 'Type', 'Line', 'Linestyle', NormalLinestyle), 'LineWidth', LineWidth)
legend('LQR', 'Pole Placement')
