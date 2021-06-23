% Proportional derivative control with TCLab.

close all
clear all
clc

% Include tclab.m for initialization.
tclab

disp('Proportional derivative control with TCLab.')
disp('LED indicates temperature.')

figure(1)
T1s = []; % Stores temperature 1 values.
T2s = []; % Stores temperature 2 values.
Q1s = []; % Stores heater 1 values.
Q2s = []; % Stores heater 2 values.
SP_T1s = []; % Stores setpoint values for heater 1.

% Initial heater values.
Q1 = 0;
Q2 = 0;
h1(Q1);
h2(Q2);

% Stores the elapsed time in each sample.
t_plot = [];

sim_time_s = [];
real_time = 0;
real_time_s = [];

ns = 600; % Number of samples.
SP_T1(1:199) = 30.0; % Setpoint 1.
SP_T1(200:ns) = 60.0; % Setpoint 2.

Pk = []; % Proportional action.

Dks = [];

% FOPTD model parameters.
K = 0.7087;
T = 148.0700;
L = 33.8800;

% Cohen and Coon tuning rules for PD controller.
Kp = (1/K)*(T/L)*(1.25 + L/(6*T));
Td = L*((6*T - 2*L)/(22*T + 3*L));

Kd = Kp*Td;

uk=0;
yk = 0;
yk_1=0;
Tf=5;
Ts=1;
Dk = 0;
Dk_1=0;
T1_1=0;

for i = 1:ns
    tic;
    
    % Read temperatures.
    T1 = T1C();   % Temperature 1.
    T2 = T2C();   % Temperature 2.
    
    ek = SP_T1(i)-T1; % Error.
    dyk = (T1-T1_1)/Ts;
    Pk = Kp*ek; % Proportional action.
    
    % Derivative action.
    filter_aux = 1; % 0 = without filter, 1 = with filter.
    if(filter_aux == 1)
        Dk = (Tf/(Tf+Ts))*Dk_1 + (Ts/(Tf+Ts))*Dk; % With filter.
    else
        Dk = Kp*Td*dyk; % Without filter.
    end
    
    uk = Pk - Dk; % PD controller output.
    
    Q1 = uk;
    
    % Check the actuator limits.
    if Q1>=100
        Q1=100;
    end
    if Q1<=0
        Q1=0;
    end
    
    % Load disturbance.
    if i>=400
        Q1=Q1-40;
    end
    
    h1(Q1);      % Send the new heater value.
    yk_1=yk;
    Dk_1=Dk;
    T1_1=T1;
    
    % LED brightness
    %     brightness1 = (t1 - 30)/50.0;  % <30degC off, >100degC full brightness
    %     brightness2 = (t2 - 30)/50.0;  % <30degC off, >100degC full brightness
    %     brightness = max(brightness1,brightness2);
    %     brightness = max(0,min(1,brightness)); % limit 0-1
    %     led(brightness);
    
    % Plot heater and temperature data.
    Q1s = [Q1s,Q1];
    Q2s = [Q2s,Q2];
    T1s = [T1s,T1];
    T2s = [T2s,T2];
    % Setpoint T1
    SP_T1s = [SP_T1s,SP_T1(i)];
    
    
    n = length(T1s);
    time = linspace(0,n+1,n);
    
    time_up=time(n);
    
    clf
    subplot(2,1,1)
    plot(time,SP_T1s,'b--','LineWidth',1);
    hold on
    plot(time,T1s,'r-','LineWidth',1);
    %plot(time,T2s,'b-','LineWidth',1);
    ylabel('Temperature (°C)')
    %legend('Set-point T1','Temperature 1','Temperature 2','Location','NorthWest')
    legend('Set-point T1','Temperature 1','Location','NorthWest')
    axis([0 time_up 0 80])
    subplot(2,1,2)
    %plot(time,Q1s,'r.','MarkerSize',8);
    plot(time,Q1s,'r-','LineWidth',1);
    hold on
    %plot(time,Q2s,'b--','LineWidth',1); % Original LineWidth = 2.
    ylabel('Heater (0-100%)')
    xlabel('Time (s)')
    %legend('Heater 1','Heater 2','Location','NorthWest')
    legend('Heater 1','Location','NorthWest')
    axis([0 time_up -10 110])
    drawnow
    
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
data = [sim_time_s',Q1s',Q2s',T1s',T2s',SP_T1s'];
if(filter_aux == 1)
    csvwrite('PD_test_2.txt',data)
else
    csvwrite('PD_test_1.txt',data)
end
