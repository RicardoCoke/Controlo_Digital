% Orginal program by John Hendrengren
% Program modified by
% Paulo Moura Oliveira


close all; 
clear all; 
clc

% include tclab.m for initialization
tclab;

disp('Test Heater 1')
disp('LED Indicates Temperature')

figure(1)
t1s = [];   % PMO: For storing Temperature T1
t2s = [];   % PMO: For storing Temperature T2
h1s = [];   % PMO: For storing Heater 1: Q1
h2s = [];   % PMO: For storing Heater 2: Q2

% PMO-  For storing the elapsed time in each sample
t_plot=[]; 


% initial heater values
ht1 = 0;
ht2 = 0;
h1(ht1);
h2(ht2);

for i = 1:399
    tic  % PMO: Saves current time
    if i==5
        disp('Turn on heater 1 to 60%')
        ht1 = 60;    % PMO: Initial Step-value to Q1=60ºC
        h1(ht1);     % PMO: Sets heater 1 Q1- Temperature 
    end
    if i==105
        disp('Turn off heater 1')
        ht1 = 0;     % PMO: Set heater Q1 value to zero
        h1(ht1);
    end
    if i==150
        disp('Turn on heater 2 to 80%')
        ht2 = 80;   % PMO: Set Middle Step-value to Q1=80ºC
        h2(ht2);
    end
    if i==205
        disp('Turn off heaters')
        ht1 = 0;
        ht2 = 0;
        h1(ht1);
        h2(ht2);
    end
    % read temperatures
    t1 = T1C();
    t2 = T2C();

    % LED brightness
    brightness1 = (t1 - 30)/50.0;  % <30degC off, >100degC full brightness
    brightness2 = (t2 - 30)/50.0;  % <30degC off, >100degC full brightness
    brightness = max(brightness1,brightness2);
    brightness = max(0,min(1,brightness)); % limit 0-1
    led(brightness);
    
    % plot heater and temperature data
    h1s = [h1s,ht1];
    h2s = [h2s,ht2];
    t1s = [t1s,t1];
    t2s = [t2s,t2];
    
    % PMO: Determines size of vector
    n = length(t1s);
    time = linspace(0,n+1,n);
    
    clf
    subplot(2,1,1)
    plot(time,t1s,'r.','MarkerSize',10);
    hold on
    plot(time,t2s,'b.','MarkerSize',10);
    ylabel('Temperature (degC)')
    legend('Temperature 1','Temperature 2','Location','NorthWest')
    subplot(2,1,2)
    plot(time,h1s,'r-','LineWidth',2);
    hold on
    plot(time,h2s,'b--','LineWidth',2);
    ylabel('Heater (0-5.5 V)')
    xlabel('Time (sec)')
    legend('Heater 1','Heater 2','Location','NorthWest')
    drawnow;
    
    t = toc    % PMO Elapsed time
    t_plot(i)=t;
    pause_time= 1.0-t   % PMO 
    pause(max(0.01,1.0-t))
    
end

disp('Turn off heaters')
h1(0);
h2(0);

disp('Heater Test Complete')

% PMO
figure
plot(time,t_plot,'bo')
ylabel('Elapsed time in each sample')
xlabel('Time (sec)')




