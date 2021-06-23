% Reads and plots TCLab on-off control test data from output file.

close all
clear
clc

%% On-off control without hysteresis band.

% Read TCLab on-off control test data from output file.
TCLab_data = readtable("OnOff_test_1.txt");
TCLab_data.Properties.VariableNames = ["Time" "H1" "H2" "T1" "T2" "Setpoint"];

% Plot data.
figure(1)

subplot(2,1,1)
plot(TCLab_data.Time,TCLab_data.Setpoint,'k--','LineWidth',1)
hold on
plot(TCLab_data.Time,TCLab_data.T1,'b-','LineWidth',2)
hold off
title('System output')
ylabel('Temperature (°C)')
legend('Setpoint','TCLab','Location','NorthWest')
axis([0 length(TCLab_data.T1) 10 80])
ax = gca;
ax.LineWidth = 2;

subplot(2,1,2)
plot(TCLab_data.Time,TCLab_data.H1,'b-','LineWidth',0.5)
title('Heater input')
xlabel('Time (s)')
ylabel('Heater (0-100%)')
legend('TCLab','Location','NorthWest')
axis([0 length(TCLab_data.T1) -10 110])
ax = gca;
ax.LineWidth = 2;

%% On-off control with hysteresis band.

% Read TCLab on-off control test data from output file.
TCLab_data = readtable("OnOff_test_2.txt");
TCLab_data.Properties.VariableNames = ["Time" "H1" "H2" "T1" "T2" "Setpoint"];

% Define the hysteresis band used.
HB_superior = TCLab_data.Setpoint + 0.7;
HB_inferior = TCLab_data.Setpoint - 0.7;

% Plot data.
figure(2)

subplot(2,1,1)
p1 = plot(TCLab_data.Time,HB_superior,'k--','LineWidth',1);
hold on
p2 = plot(TCLab_data.Time,HB_inferior,'k--','LineWidth',1);
plot(TCLab_data.Time,TCLab_data.T1,'b-','LineWidth',2)
hold off
title('System output')
ylabel('Temperature (°C)')
set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
legend('Hysteresis band (setpoint ± 0.7°C)','TCLab','Location','NorthWest')
axis([0 length(TCLab_data.T1) 10 80])
ax = gca;
ax.LineWidth = 2;

subplot(2,1,2)
plot(TCLab_data.Time,TCLab_data.H1,'b-','LineWidth',0.5)
title('Heater input')
xlabel('Time (s)')
ylabel('Heater (0-100%)')
legend('TCLab','Location','NorthWest')
axis([0 length(TCLab_data.T1) -10 110])
ax = gca;
ax.LineWidth = 2;
