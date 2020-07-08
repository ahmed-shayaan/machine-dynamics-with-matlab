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
% Creation Date: 8 July 2020
% MATLAB Version: R2020a
%
%% REVISION
%
% V1.0 | 8 July 2020 | shayaan.ahmed@vitap.ac.in
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

%% 2.) Computing
w0 = [x_0 x_dot_0];

[t_sim,w_sim] = ode45(@state_space_eqn,time,w0,'options',mass,stiffness,damping,force,omega);

time = t_sim';
x_t = w_sim(:,1);
v_t = w_sim(:,2);
plot(time,x_t,'r',time,v_t,'b','LineWidth',1)
xlabel('Time (s)')
ylabel('Frequency')
legend('x_t','v_t')
title('Displacement and Velociity vs Time')


