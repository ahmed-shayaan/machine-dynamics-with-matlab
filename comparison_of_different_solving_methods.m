%% DESCRIPTION
%
% This is a Script for solving the differential eqns of force excited
% single mass system
%
%% OUTPUT
%
% Formatted figure of displacement of a single mass system and its
% animation
%
%% VERSION
%
% Author: Shayaan Ahmed
% Creation Date: 09 July 2020
% MATLAB Version: R2020a
%
%% REVISION
%
% V1.0 | 9 July 2020 | shayaan.ahmed@vitap.ac.in
% Copyright (c) 2020 Shayaan Ahmed
%
%% Program
clear       % Clears the workspace
clc         % Clears the command window
close all   % Closes all figures

%% 1.) Definition
%% 1.)-Parametric Definition
mass = 750;            % mass of the body(kg)
stiffness = 50000;       % stiffness of the spring (N/m)
damping = 10;         % damping coeff of the damper (Ns/m)
time = 0:0.01:10;            % time (s)\

x_0 = 0.00;             % initial cond. for displacement (m)
x_dot_0 = 0.00;         % initial cond. for velocity (m/s)

force = 2000;           % force applied on the body (N)
omega = 1.88;           % Angular freq. of the  excitation (rad/s)

n = 10;
%% 2.) Computing
%%% Analytically%%%
for k = 1:n
    tic
    %% 2.)- Analytically | Parametric calculation
    dampingcoeff = damping/(2*mass); % calculation of damping coefficient
    ang_eig_freq = sqrt(stiffness/mass);    % calculation of angular eigen freq.
    
    %%2.)- Analytically | Calculation of Characteristic Polynomial
    lambda = roots([1, 2*damping, ang_eig_freq^2]); % calc of roots of eqn
    
    %% 2.)- Analytically | Calc. of Constants
    x_roof = sqrt(force^2/((damping*omega)^2 + (stiffness - mass*omega^2)^2)); % max amplitude of particular soln
    phi = atan((damping*omega)/(stiffness - mass*omega^2));
    
    k1 = x_0 -(lambda(1)*x_0 - x_dot_0 - x_roof*omega*sin(-phi) - x_roof*lambda(1)*cos(-phi))/(lambda(1) - lambda(2)); % integration constsnt
    k2 = (lambda(1)*x_0 - x_dot_0 - x_roof*omega*sin(-phi) - x_roof*lambda(1)*cos(-phi))/(lambda(1) - lambda(2)); % Integration Constant
    
    %% 2.)- Analytically | Calc. of Solution
    
    x_t_h = k1*exp(lambda(1)*time) + k2*exp(lambda(2)*time); % Homogeneous soln for the displacement
    v_t_h = k1*lambda(1)*exp(lambda(1)*time) + k2*lambda(2)*exp(lambda(2)*time); % Homogeneous soln for the velocity
    
    x_t_p = x_roof*cos(omega*time - phi); % Particular soln for the displacement
    v_t_p = -x_roof*omega*cos(omega*time - phi); % Particular soln for the velocity
    
    x_t = real(x_t_h) + x_t_p; % overall soln for displacement
    v_t = real(v_t_h) + v_t_p; % overall soln for velocity
    
    time_analytically(n) = toc;
end
%%% dsolve method %%%
%% 2.)- dsolve | Calculation of the Soln
for kk = 1:n
    tic
        syms x(t)           % defining dependent variable
        Dx = diff(x,1);     % Defining first derivation
        D2x = diff(x,2);    % defining second derivation
        x = dsolve(mass*D2x + damping*Dx + stiffness*x == force*cos(omega*t), x(0) == x_0,...
            Dx(0) == x_dot_0, 't');
    time_dsolve(n) = toc;
        
     %% 2.) dsolve | Evaluate the eqn
     x_fun = matlabFunction(x);         % creating function handle for x
     x_dot_fun = matlabFunction(diff(x));   % creating function handle for x_dot
     x_t = feval(x_fun,time);           % evaluate functions at points "time"
     v_t = feval(x_dot_fun,time);       % evaluate function at points "time"
end

%%% ode %%%
%% 2.)- ode | Set Initial Condition
for kkk = 1:n
    tic
w0 = [x_0 x_dot_0]; % creating vector with initial conditions

    %% 2.)- ode | Numerical Solution of motion
    [tsim wsim] = ode45(@state_space_eqn,time,w0,'options',mass,stiffness,damping,force,omega);
    time_ode(n) = toc;
    time_for_ode = tsim';
    x_t_ode = wsim(:,1)';
    v_t_ode = wsim(:,2)'; 
    
end

compare_analytical_dsolve = sum(time_dsolve)/sum(time_analytically);
compare_analytical_ode = sum(time_ode)/sum(time_analytically);

fprintf('Average calculation time for analytical solution %d seconds \n',sum(time_analytically)/n);
fprintf('Average calculation time for dsolve solution %d seconds \n',sum(time_dsolve)/n);
fprintf('Average calculation time for ode solution %d seconds \n',sum(time_ode)/n);

fprintf('Average dsolve solution takes %d longer than analytical solution \n',compare_analytical_dsolve);
fprintf('Average ode solution takes %d longer than analytical solution \n',compare_analytical_ode);