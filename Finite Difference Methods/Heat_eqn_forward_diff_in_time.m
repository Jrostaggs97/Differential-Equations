clc;
clear;

%step size for space and time
delta_x = .02;
delta_t = .0001;
cfl = delta_t/(delta_x^2);

%set up grid
x_grid = 0:delta_x/(1-delta_x):1;
t_grid = 0:delta_t/(1-delta_t):1;

%spatial finite differencing matrix 
%(2nd order centered difference for u_xx)
J = length(x_grid);
T = length(t_grid);
v = ones(J,1);
A = cfl*spdiags([v,-2*v,v],-1:1,J-2,J-2); %solving inside domain and will deal with boundary conditions later
A= full(A); %make nonspare for debugging


%matrix that will hold soln at each time step
U = zeros(J,T); 
%notice that we are enforcing the zero boundary conditions
%so we will only be solving "inside" the domain and 
%artficially tack on the boundary conditions at the end

%initial condition
U(2:J-1,1) = (.5 -abs(x_grid(2:J-1) -.5));
%plot(U(:,1))

%matrix to store solutions at each time
%each time step is store in the column

for t =1:T-1
    
   U(2:J-1,t+1) = U(2:J-1,t) + A*U(2:J-1,t); 
   
   
end



xlim([0,1])
ylim([0,.5]) 
T_list = [1,50,100,350,2500,6000,T];
legendStrings = "T = " +string(T_list);

for t = T_list
    plot(x_grid',U(:,t))
    hold on;

end

xlabel("x")
ylabel("u(x,t)")
lgd = legend(legendStrings)
lgd.FontSize = 23
title("Numerical Approx. of $$u_t = u_{xx}, \Delta t = .0001,\Delta x = .02$$", "interpreter","latex")
