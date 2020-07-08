%% DESCRIPTION
%
% This is a Script for solving the differential eqn of a forced excited
% single mass system.
%
%% OUTPUT
%
% Formatted figure of the displacement of a single mass system and its
% animation
%
%% VERSION
%
% Author: Shayaan Ahmed
% Creation Date: 7 July 2020
% MATLAB Version R2020a

%% REVISION
% 
%V1.0 | 7 July 2020 | shayaan.ahmed@vitap.ac.in | 
% Copyright (c) 2020 Shayaan Ahmed
%
%% Progress
clear           % Clears all variables in the workspace
clc             % Clears the command window
close all       % Closes all the plots and diagrams if any
%
%% 1.) Definition
%% 1.)- Parametric Definition
mass = 2600;            % mass of the body in kg
stiffness = 75000;       % Stifness of the spring in N/m
damping = 25;         % damping coeff. of the damper in Ns/m
time_period = 0:0.1:15;     % time in seconds

x_0 = 0;             % initial displacement in m
x_dot_0 = 0;         % initial velocity in m/s

force = 5500;       % force excitation on th system in N
omega = 2;          % Angular freq of the system in rad/sec
%
%% 2.) Computing
%% 2.)- Parametric Calculation
damping_coeff = damping/(2*mass);
angular_eig_freq = sqrt(stiffness/mass);

lambda = roots([1,2*damping_coeff,angular_eig_freq^2]); % eigen values of the solution
x_roof = sqrt(force^2/((damping*omega)^2+(stiffness-mass*omega^2)^2)); % amplitude of the system
phi = atan((damping*omega)/(stiffness-(mass*omega^2)));

K1 = x_0-(lambda(1)*x_0 - x_dot_0 - x_roof*omega*sin(-phi) - x_roof*lambda(1)*cos(-phi))/(lambda(1) - lambda(2)); % Integration constant
K2 = (lambda(1)*x_0 - x_dot_0 - x_roof*omega*sin(-phi) - x_roof*lambda(1)*cos(-phi))/(lambda(1) - lambda(2)); % Integration Constsnnt

x_h = K1*exp(lambda(1)*time_period) + K2*exp(lambda(2)*time_period); % homogeneous eqn for displacement
v_h = K1*lambda(1)*exp(lambda(1)*time_period) + K2*lambda(2)*exp(lambda(2)*time_period); % homogeneous eqn for velocity

x_p = x_roof*cos(omega*time_period - phi); % particular soln for the displacement
v_p = -x_roof*omega*sin(omega*time_period - phi); % particular soln for the velocity

x_t = real(x_h + x_p); % overall soln for displacement
v_t = real(v_h + v_p); % overall soln for velocity

P1 = plot(x_t)
P1.LineStyle = '--'
P1.LineWidth = 2.0
P1.Color = 'r'
hold on
P2 = plot(v_t)
P2.LineWidth = 1.3
P2.Color = 'g'
P2.Marker = 'x'
P2.MarkerFaceColor = 'b'
hold off
title('Displacement and Velocity Plot for Forced Excitation')
legend('x_t','v_t')