%% Description of the MATLAB File
%
% This is a Script to solve the differential equation of a single mass (1 DOF System)
% system.
%
%% 
%
% Formatted figure of the displacement of a single mass system and its
% animation.
%
%% Version
% Author: Shayaan Ahmed (ahmedshayaan7030@gmail.com)
% Creation date: 27-June-2020
% Matlab Version: R2020a
% Copyright (c) 2020 Shayaan Ahmed 
% 
%% Revision 
%
% Ver 1.0 | 27-June-2020 | Shayaan Ahmed | creation

%% Program
clear % Delete Workspace
close all % Closes all figures
clc % clears command window

%% 1.) Definition
%% 1.)- Parameter Definition
mass = 1150; % mass of the body (kg)
stiffness = 78000; % (N/m)
time = 0:0.01:1; % time(seconds)
x_0 = 0.01; % (meters)
x_dot_0 = 0.1; % (m/s)

%%- Symbolic Function definition
syms x(t) %gives symblic relation
Dx = diff(x,1);
D2x =  diff(x,2);

%% 2.) Computing
%% 2.) - Solve the equation
x = dsolve((mass*D2x)+(stiffness*x)==0, x(0)== x_0, Dx(0)== x_dot_0,'t');

%% 2.)- Evaluate the eqn
x_func = matlabFunction(x);
x_dot_func = matlabFunction(diff(x));

x_t = x_func(time);
v_t = x_dot_func(time);

%% 2.)- Calculate amplitude
x_roof = max(abs(x_t)); % amplitude of the systemclear
amp_display = sprintf('Amplitude of the system = %f',x_roof);
disp(amp_display)

%% 2.)- Calculate time period T
[maxima, max_location] = findpeaks(x_t,time);
time_period = diff(max_location);
t_p = sprintf('Time Period of the system = %f',time_period);
disp(t_p)
%% 2.)- Calculate eigenfrequency and angular eigenfrequency
eigen_freq = 1/time_period; % units in Hz
angular_eigen_freq = 2*pi*eigen_freq; % units in cycles/sec
e_f = sprintf('Eigen Frequency of the system = %f',eigen_freq);
disp(e_f)
a_e_f = sprintf('Angular Eigen Frequency of the system = %f',angular_eigen_freq);
disp(a_e_f)

%% 2.) - Calculate Phase angle
temp_variable = diff(sign(v_t));
index_up = find(temp_variable>0);
index_down = find(temp_variable<0);
index_up_time = time(index_up);
index_down_time = time(index_down);
first_zero_crossing = min([index_up_time index_down_time]);

% x_dot_t = -amplitude*ang_freq*sin((ang_freq*t(phase))-phase_angle) is the given
% equation of motion
% At first crossing x_dot_t=0 meaning
% -amplitude*ang_freq*sin((ang_freq*t(phase))-phase_angle) = 0
% But amplitude and ang_freq cannot be 0. Hence sin((ang_freq*t(phase))-phase_angle) = 0
% Hence, ang_freq*t(phase)-phase_angle=0 => phase_angle = ang_freq*t(phase)

phase_angle_radians = angular_eigen_freq*first_zero_crossing; % phase angle in radians
phase_angle_degrees = (180/pi)*phase_angle_radians; % phase angle in degrees
phase_ang = sprintf('Phase Angle of the system in degrees = %f',phase_angle_degrees);
disp(phase_ang)

%% 3.)Plots
%% 3.)-Plot Parameters
clr = [236/255 237/255 237/255]; %defines color of the plot curves
units = 'normalized';
lnwdth = 2; % defining the line width
fntsz = 22; % defines the font size 
pos_fig = [0.01 0.1 0.98 0.8]; % (1st val is distance from left of screen, 2nd value is dist. from bottom of screen)
                              % (3rd value is length b/w the rt and lt ends of screen and 4th val is the length b/w bottom & top)
title_graph = 'Displacement and velocity vs. Time';
xlabel_graph = 'Time t (s)';
ylabel_graph{1} = 'Displacement x (m)';
ylabel_graph{2} = 'velocity v (m/s)';

%% 3.)-Plot Graph
fig = figure('color',clr,'units',units,'position',pos_fig);
[axes_graph, Line1, Line2] = plotyy(time,x_t,time,v_t); 
set(Line1,'color','k','Linewidth',lnwdth);
set(Line2,'color','r','Linewidth',lnwdth);
set(axes_graph(1),'Ycolor','k','linewidth',lnwdth,'fontsize',fntsz);
set(axes_graph(2),'Ycolor','r','linewidth',lnwdth,'fontsize',fntsz);
xlabel(axes_graph(1),xlabel_graph,'fontsize',fntsz);
ylabel(axes_graph(1),ylabel_graph(1),'fontsize',fntsz);
ylabel(axes_graph(2),ylabel_graph(2),'fontsize',fntsz);

title(title_graph);

x_t_max_limit = max(abs(x_t))+0.08*max(abs(x_t));
ylim(axes_graph(1),[-x_t_max_limit,x_t_max_limit]);

v_t_max_limit = max(abs(v_t))+0.08*max(abs(v_t));
ylim(axes_graph(2),[-v_t_max_limit,v_t_max_limit]);

legend('displacement curve', 'velocity curve');