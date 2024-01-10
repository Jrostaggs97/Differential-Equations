clc;
clear;

%step size for space and time
delta_x = .1;
delta_t = .05;
cfl = delta_t/(delta_x^2);

%set up grid
x_grid = 0:delta_x/(1-delta_x):1;
t_grid = 0:delta_t/(1-delta_t):1;