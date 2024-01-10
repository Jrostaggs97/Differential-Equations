clc; clear;

m = 100;
alpha = .7;
beta = .5;
T = 2*pi;
h=T/(m+1);
x=linspace(0,T,m+2);

theta =.7 + sin(x/2); %initial theta
%alpha*cos(x) + beta*sin(x)
%



% 
for k = 1:60
    
    G_vec = G(theta,h,m); %vector valued function containing non linear function data
    Jacob = pendulum_jacobian(theta,h,m); %Jacobian of G

    theta(2:m+1) = theta(2:m+1)' - Jacob\G_vec; %starting 2 and ending at m+1 enforces boundary conditions. sin(0)=cos(2pi) = 1 and sin(2pi)=cos(0)=0
    
    if max(abs(Jacob\G_vec)) <10^(-14) %infinity norm for stopping point
       break  
    end
    
    %plotting
%     plot(x,theta)
%     hold on
%     xlabel("theta (radians)")
%     ylabel("pendulum amplitude")
%     title("Pendulum amplitude vs theta, Newton iterates")
%     xlabel("Time")
%     ylabel("Angle (radians)")
%     title("Theta vs. T")

end
% lgd =legend("k=1","k=2","k=3","k=4","k=5","k=6")
% lgd.FontSize = 20

plot(x,theta)
title("first")


T2 = 8*pi;
h2=T2/(m+1);
x2=linspace(0,T2,m+2);

theta_u = theta;

for k = 1:100
    
    G_vec2 = G(theta_u,h2,m); %vector valued function containing non linear function data
    Jacob2 = pendulum_jacobian(theta_u,h2,m); %Jacobian of G

    theta_u(2:m+1) = theta_u(2:m+1)' - Jacob2\G_vec2; %starting 2 and ending at m+1 enforces boundary conditions. sin(0)=cos(2pi) = 1 and sin(2pi)=cos(0)=0
    
    if max(abs(Jacob2\G_vec2)) <10^(-14) %infinity norm for stopping point
       break  
    end
    
    %plotting
%     plot(x,theta)
%     hold on
%     xlabel("theta (radians)")
%     ylabel("pendulum amplitude")
%     title("Pendulum amplitude vs theta, Newton iterates")
%     xlabel("Time")
%     ylabel("Angle (radians)")
%     title("Theta vs. T")

end
max_theta = max(theta_u);
plot(x2,theta_u)
xlabel("Time")
ylabel("Angle (radians)")
title("Theta vs T")
