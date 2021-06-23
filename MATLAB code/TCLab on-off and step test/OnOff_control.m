% Teste On-Off sem banda de histerese.

close all
clear
clc

% Include tclab.m for initialization.
tclab

disp('On Off Control')
disp('LED Indicates Temperature')

figure(1)
% Stores values
t1s = [];
t2s = [];
h1s = [];
h2s = [];
SP_T1s = [];    % Stores set-point 1 values
% initial heater values
ht1 = 0;
ht2 = 0;
h1(ht1);
h2(ht2);

% set point (°C)
%t1_sp = 35;
ns=600;                  % Number of samples
SP_T1(1:300)  = 30.0;    % Set-point step 1
SP_T1(300:ns) = 60.0;    % Set-point step 2

for i = 1:ns
    tic;

    % read temperatures
    t1 = T1C();
    t2 = T2C();

    % LED brightness
    brightness1 = (t1 - 30)/50.0;  % <30°C off, >100°C full brightness
    brightness2 = (t2 - 30)/50.0;  % <30°C off, >100°C full brightness
    brightness = max(brightness1,brightness2);
    brightness = max(0,min(1,brightness)); % limit 0-1
    led(brightness);
    
    % On Off Control
    %if t1>t1_sp
    if t1 > SP_T1(i)    
       ht1 = 0;
    else
       ht1 = 100;
    end
    h1(ht1);
    
    % plot heater and temperature data
    h1s = [h1s,ht1];    % Q1
    h2s = [h2s,ht2];    % Q2
    t1s = [t1s,t1];     % T1
    t2s = [t2s,t2];     % T2
    
    % Set-point 1
    SP_T1s = [SP_T1s,SP_T1(i)];  
    
    n = length(t1s);
    time = linspace(0,n+1,n);
    
    clf
    subplot(2,1,1)
    plot(time,t1s,'r.','MarkerSize',10);
    hold on
    plot([0,n+1],[SP_T1(i),SP_T1(i)],'b-','LineWidth',1); % Original LineWidth = 2.
    ylabel('Temperature (°C)')
    legend('Temperature 1','Temp 1 Set Point','Location','NorthWest')
    subplot(2,1,2)
    plot(time,h1s,'r-','LineWidth',1); % Original LineWidth = 2.
    hold on
    plot(time,h2s,'b--','LineWidth',1); % Original LineWidth = 2.
    ylabel('Heater (0-100%)')
    xlabel('Time (s)')
    legend('Heater 1','Heater 2','Location','NorthWest')
    drawnow;
    
    t = toc;
    
    sim_time_s(i)=i;
    
    pause(max(0.01,1.0-t))
end

disp('Turn off heaters and LED')
h1(0);
h2(0);

led(0); % Desligar o LED.

disp('Heater Test Complete')

% Save .txt file with data
data = [sim_time_s',h1s',h2s',t1s',t2s',SP_T1s'];
csvwrite('On_Off_test_1.txt',data);
