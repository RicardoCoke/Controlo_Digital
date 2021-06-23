% Reads and plots TCLab open-loop step response data from output file.

close all
clear
clc

% Read TCLab open-loop step response data from output file.
TCLab_data = readtable("Step_test_1.txt");
TCLab_data.Properties.VariableNames = ["Time" "H1" "H2" "T1" "T2"];

% Plot data.
figure(1)

subplot(2,1,1)
plot(TCLab_data.Time,TCLab_data.T1,'b-','LineWidth',2)
title('Open-loop step response')
ylabel('Temperature (°C)')
legend('TCLab','Location','SouthEast')
axis([0 length(TCLab_data.T1) 10 80])
ax = gca;
ax.LineWidth = 2;

subplot(2,1,2)
plot(TCLab_data.Time,TCLab_data.H1,'b-','LineWidth',0.5)
title('Heater input')
xlabel('Time (s)')
ylabel('Heater (0-100%)')
legend('TCLab','Location','SouthEast')
axis([0 length(TCLab_data.T1) -10 110])
ax = gca;
ax.LineWidth = 2;
