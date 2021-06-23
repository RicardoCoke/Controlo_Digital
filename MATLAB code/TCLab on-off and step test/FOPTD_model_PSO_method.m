% Reads TCLab open-loop step response data from output file and calculates
% the first order plus time delay (FOPTD) model parameters (K, T, L) using
% a particle swarm optimization (PSO) algorithm.

close all
clear
clc

%% FOPTD model identification: PSO method.

fun = @PSO_min_fun;
lb = [0.1,120,16];
ub = [1.6,200,45];
options = optimoptions('particleswarm','SwarmSize',20,'HybridFcn',@fmincon,'MaxIterations',40);

rng default % For reproducibility.
nvars = 3;
x = particleswarm(fun,nvars,lb,ub,options); % Run PSO algorithm.

fprintf('FOPTD model parameters by PSO method:\nK: %.4f\nT: %.4f\nL: %.4f\n',x(1),x(2),x(3))

function ISE_pso = PSO_min_fun(x)
K = x(1);
T = x(2);
L = x(3);

SIMULINK_FOPTD_MODEL = 'FOPTD_model_simulink';

% Read TCLab open-loop step response data from output file.
TCLab_data = readtable("Step_test_1.txt");
TCLab_data.Properties.VariableNames = ["Time" "H1" "H2" "T1" "T2"];

% Time at which the step input was applied.
TIME_STEP = 10;

simulation_time = length(TCLab_data.Time);

% Initial temperature reading.
y_initial = TCLab_data.T1(TIME_STEP);

% Assign Simulink model parameters.
assignin('base','K',K)
assignin('base','T',T)
assignin('base','L',L)
assignin('base','y_initial',y_initial)

% Run Simulink FOPTD model open-loop step response.
[~,~,y1] = sim(SIMULINK_FOPTD_MODEL,simulation_time);

y_pso = y1(1:simulation_time,2);

e1 = y_pso - TCLab_data.T1(1:simulation_time);
ISE_pso = sum(e1.*e1);
end
