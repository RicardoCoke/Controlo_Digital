% Reads data from PID control output file, simulates FOPTD model response,
% and plots the real data and the simulation data for comparison.

close all
clear
clc

% Read TCLab test data from output file.
TCLab_data = readtable("PID_test_ITAE.txt");
TCLab_data.Properties.VariableNames = ["Time" "H1" "H2" "T1" "T2" "Setpoint"];

y_initial = TCLab_data.T1(1);

%% FOPTD model parameters (Two-points).
% K = 0.7087;
% T = 148.0700;
% L = 33.8800;

%% FOPTD model parameters (PSO).
K = 0.7243;
T = 162.4602;
L = 33.3980;

%% Cohen and Coon tuning rules for PID controller.
% Kp = (1/K)*(T/L)*(4/3 + L/(4*T));
% Ti = L*((32 + 6*(L/T))/(13 + 8*(L/T)));
% Td = L*(4/(11 + 2*(L/T)));

%% IAE.
% Kp = (1.435/K)*((T/L)^(0.921));
% Ti = 1/((0.878/T)*((L/T)^(-0.749)));
% Td = T*0.482*((L/T)^(1.137));

%% ITAE.
Kp = (1.357/K)*((T/L)^(0.947));
Ti = 1/((0.842/T)*((L/T)^(-0.738)));
Td = T*0.381*((L/T)^(0.995));

%% AMIGO. - mudar Simulink
% Kp = (1/K)*(0.2+0.45*(T/L));
% Ti = ((0.4*L+0.8*T)/(L+0.1*T))*L;
% Td = (0.5*L*T)/(0.3*L+T);

%% S-IMC. - mudar Simulink
% Kpp = 6.7954;
% Tii = 135.72;
% Tdd = 7.54;
% f = 1+(Tdd/Tii); % Fator de conversão para PID.
% Kp = Kpp*f;
% Ti = Tii*f;
% Td = Tdd/f;

%% Run Simulink model and get simulation data.
ObjFunc = 'PID_control_simulink';
sim_time = 600;
[t,x,y] = sim(ObjFunc,sim_time); % Read data.

%%

% Store simulation data.
r_P_sim = y(1:sim_time,1);
u_P_sim = y(1:sim_time,6);
y_P_sim = y(1:sim_time,3);
t_P_sim = y(1:sim_time,4);

% Plot overlayed results.
figure(1)

subplot(2,1,1)
plot(TCLab_data.Time,TCLab_data.Setpoint,'k--','LineWidth',1)
hold on
plot(TCLab_data.Time, TCLab_data.T1,'b-','LineWidth',2)
plot(y_P_sim,'r-','LineWidth',2)
hold off
title('System output')
ylabel('Temperature (°C)')
legend('Setpoint','TCLab','FOPTD model','Location','NorthWest')
axis([0 sim_time 10 80])
ax = gca;
ax.LineWidth = 2;

subplot(2,1,2)
plot(TCLab_data.Time,TCLab_data.H1,'b-','LineWidth',0.5)
hold on
plot(u_P_sim,'r-','LineWidth',1)
hold off
title('PID controller output')
xlabel('Time (s)')
ylabel('Heater (0-100%)')
legend('TCLab','FOPTD model','Location','NorthWest')
axis([0 sim_time -10 110])
ax = gca;
ax.LineWidth = 2;
