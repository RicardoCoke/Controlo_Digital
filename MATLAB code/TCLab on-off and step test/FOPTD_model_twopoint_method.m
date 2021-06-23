% Reads TCLab open-loop step response data from output file and calculates
% the first order plus time delay (FOPTD) model parameters (K, T, L) using
% the two-points method.

close all
clear
clc

% Read TCLab open-loop step response data from output file.
TCLab_data = readtable("Step_test_1.txt");
TCLab_data.Properties.VariableNames = ["Time" "H1" "H2" "T1" "T2"];

% Time at which the step input was applied.
TIME_STEP = 10;

%% FOPTD model identification: two-point method.

% Initial input value.
u_initial = TCLab_data.H1(1);
% Initial output value.
y_initial = TCLab_data.T1(TIME_STEP);

% Final input value in steady state.
u_final = max(TCLab_data.H1);
% Final output value in steady state.
y_final = max(TCLab_data.T1);

% Find the two points (output at 35.2% and 85.3% of final value).
y_35 = y_initial + 0.352*(y_final - y_initial);
y_85 = y_initial + 0.853*(y_final - y_initial);

y_35_index = find(TCLab_data.T1 >= y_35,1,'first');
y_35_time = TCLab_data.Time(y_35_index);

y_85_index = find(TCLab_data.T1 >= y_85,1,'first');
y_85_time = TCLab_data.Time(y_85_index);

% Calculate FOPTD model parameters. K: DC gain. T: time constant. L: time
% delay.
K = (y_final - y_initial)/(u_final - u_initial);
T = 0.67*((y_85_time - TIME_STEP) - (y_35_time - TIME_STEP));
L = 1.3*(y_35_time - TIME_STEP) - 0.29*(y_85_time - TIME_STEP);

fprintf('FOPTD model parameters by the two-point method:\nK: %.4f\nT: %.4f\nL: %.4f\n',K,T,L)

%% Plot TCLab data with the two points highlighted.

figure(1)
plot(TCLab_data.Time,TCLab_data.T1,'b-','LineWidth',2)
hold on
plot(y_35_time,y_35,'ko','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','r')
plot(y_85_time,y_85,'ks','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','r')
hold off
title('FOPTD model identification: two-point method')
xlabel('Time (s)')
ylabel('Temperature (°C)')
legend('TCLab','1st point (35.2%)','2nd point (85.3%)','Location','SouthEast')
axis([0 length(TCLab_data.T1) 10 80])
ax = gca;
ax.LineWidth = 2;
