% Runs TCLab open-loop step response test and outputs a file with the
% resulting data.

close all
clear
clc

% Include tclab.m for initialization.
tclab

disp('Step response test with TCLab.')
disp('LED indicates temperature.')

figure(1)
t1s = []; % Stores temperature 1 values.
t2s = []; % Stores temperature 2 values.
h1s = []; % Stores heater 1 values.
h2s = []; % Stores heater 2 values.

% Initial heater values.
ht1 = 0;
ht2 = 0;
h1(ht1);
h2(ht2);

% Stores the elapsed time in each sample.
t_plot = [];

sim_time_s = [];
real_time = 0;
real_time_s = [];

for i = 1:600
    tic;
    if i==10
        disp('Turn on heater 1.')
        ht1 = 80;
        h1(ht1);
    end
    
    % read temperatures
    t1 = T1C();
    t2 = T2C();
    
    % LED brightness
    %     brightness1 = (t1 - 30)/50.0;  % <30°C off, >100°C full
    %     brightness brightness2 = (t2 - 30)/50.0;  % <30°C off, >100°C
    %     full brightness brightness = max(brightness1,brightness2);
    %     brightness = max(0,min(1,brightness)); % limit 0-1
    %     led(brightness);
    
    % plot heater and temperature data
    h1s = [h1s,ht1];
    h2s = [h2s,ht2];
    t1s = [t1s,t1];
    t2s = [t2s,t2];
    n = length(t1s);
    time = linspace(0,n+1,n);
    
    clf
    subplot(2,1,1)
    plot(time,t1s,'r.','MarkerSize',10);
    hold on
    plot(time,t2s,'b.','MarkerSize',10);
    ylabel('Temperature (°C)')
    legend('Temperature 1','Temperature 2','Location','NorthWest')
    subplot(2,1,2)
    plot(time,h1s,'r-','LineWidth',1); % Original LineWidth = 2.
    hold on
    plot(time,h2s,'b--','LineWidth',1); % Original LineWidth = 2.
    ylabel('Heater (0-100%)')
    xlabel('Time (s)')
    legend('Heater 1','Heater 2','Location','NorthWest')
    drawnow;
    t = toc;
    t_plot(i)=t;
    real_time=real_time+t;
    real_time_s(i)=real_time;
    
    sim_time_s(i)=i;
    
    pause(max(0.01,1.0-t))
end

disp('Turn off heaters.')
h1(0)
h2(0)

led(0) % Turn off LED.

figure
plot(time,t_plot,'bo')
ylabel('Elapsed time in each sample')
xlabel('Time (s)')

disp('Heater test complete.')

% Save txt file with data.
data = [sim_time_s',h1s',h2s',t1s',t2s'];
csvwrite('Step_test_1.txt',data)
