
clc; clear; 

% Step sizes 
delta_x = .05; 
delta_t = .04; 
cfl = delta_t/(2*(delta_x));

% Set up grids
x_grid = 0:delta_x/(1-delta_x):5;
t_grid = 0:delta_t/(1-delta_t):5;
J = length(x_grid);
T = length(t_grid);

% Create solution matrix
U = zeros(T,J);

% Insert the initial condition
for j = 1:J
    if x_grid(j)<1
    U(1,j) = 1;
    else 
    U(1,j) = 0;
    end
end

% Insert the boundary condition
for n = 1:T
    U(n,1) = 1;
end


% Iterate in time?
for n = 1:T
    % Iterate in space
    for j = 2:J
        U(n+1,j) = U(n,j) -(delta_t/(2*delta_x))*(U(n,j)^2 - U(n,j-1)^2);
    end
end

for k = 1:20:T
    hold on;
    plot(x_grid,U(k,:))
end


title("Approximation of $$u_t + (\frac{u^2}{2})_x=0$$", "Interpreter","latex")
xlabel("x")
ylabel("Approximate Solution, u(x)")



lgd = legend("t = 1","t = 2","t = 3","t = 4","t = 5")
lgd.FontSize = 20
