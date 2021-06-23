% Reads data from PI control output file, simulates FOPTD model response,
% and plots the real data and the simulation data for comparison.

close all
clear all
clc

% Read TCLab test data from output file.
TCLab_data = readtable("PI_test_1.txt");
TCLab_data.Properties.VariableNames = ["Time" "H1" "H2" "T1" "T2" "Setpoint"];

y_initial = TCLab_data.T1(1);

% FOPTD model parameters.
K = 0.7087;
T = 148.0700;
L = 33.8800;

% Cohen and Coon tuning rules for PI controller.
Kp = (1/K)*(T/L)*(0.9 + L/(12*T));
Ti = L*((30*T + 3*L)/(9*T + 20*L));

Ki = Kp/Ti;

% Run Simulink model and get simulation data.
ObjFunc = 'PI_control_simulink';
sim_time = 600;
[t,x,y] = sim(ObjFunc,sim_time); % Read data.

% Store simulation data.
r_P_sim = y(1:sim_time,1);
u_P_sim = y(1:sim_time,5);
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
title('PI controller output')
xlabel('Time (s)')
ylabel('Heater (0-100%)')
legend('TCLab','FOPTD model','Location','NorthWest')
axis([0 sim_time -10 110])
ax = gca;
ax.LineWidth = 2;
