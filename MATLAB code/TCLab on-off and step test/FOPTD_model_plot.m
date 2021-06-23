% Reads TCLab open-loop step response data from output file, simulates the
% two-point and PSO FOPTD models using Simulink, plots the three sets of
% data overlayed, and calculates the integral square error (ISE) of each
% model.

close all
clear
clc

SIMULINK_FOPTD_MODEL = 'FOPTD_model_simulink';

% Read TCLab open-loop step response data from output file.
TCLab_data = readtable("Step_test_1.txt");
TCLab_data.Properties.VariableNames = ["Time" "H1" "H2" "T1" "T2"];

% Time at which the step input was applied.
TIME_STEP = 10;

simulation_time = length(TCLab_data.Time);

% Initial temperature reading.
y_initial = TCLab_data.T1(TIME_STEP);
assignin('base','y_initial',y_initial)

%% Simulate two-point FOPTD model.

% K: DC_gain. T: time constant. L: time delay.
K = 0.7087;
T = 148.0700;
L = 33.8800;

% Assign parameters for Simulink.
assignin('base','K',K)
assignin('base','T',T)
assignin('base','L',L)

% Run Simulink FOPTD model open-loop step response.
[~,~,y1] = sim(SIMULINK_FOPTD_MODEL,simulation_time);

y_2pt = y1(1:simulation_time,2);
t_2pt = y1(1:simulation_time,3);

e1 = y_2pt - TCLab_data.T1(1:simulation_time);
ISE_2pt = sum(e1.*e1);

fprintf('Two-point FOPTD model ISE: %.4f\n',ISE_2pt)

%% Simulate PSO FOPTD model.

% K: DC_gain. T: time constant. L: time delay.
K = 0.7243;
T = 162.4602;
L = 33.3980;

% Assign parameters for Simulink.
assignin('base','K',K)
assignin('base','T',T)
assignin('base','L',L)

% Run Simulink FOPTD model open-loop step response.
[~,~,y1] = sim(SIMULINK_FOPTD_MODEL,simulation_time);

y_pso = y1(1:simulation_time,2);
t_pso = y1(1:simulation_time,3);

e1 = y_pso - TCLab_data.T1(1:simulation_time);
ISE_pso = sum(e1.*e1);

fprintf('PSO FOPTD model ISE: %.4f\n',ISE_pso)

%% Plot TCLab and simulation data.

figure(1)
plot(TCLab_data.Time,TCLab_data.T1(1:simulation_time),'b-','LineWidth',2)
hold on
plot(t_2pt,y_2pt(1:simulation_time),'r-','LineWidth',2)
plot(t_pso,y_pso(1:simulation_time),'g-','LineWidth',2)
hold off
title('Open-loop step response')
xlabel('Time (s)')
ylabel('Temperature (°C)')
legend('TCLab','Two-point FOPTD model','PSO FOPTD model','Location','SouthEast')
axis([0 simulation_time 10 80])
ax = gca;
ax.LineWidth = 2;
