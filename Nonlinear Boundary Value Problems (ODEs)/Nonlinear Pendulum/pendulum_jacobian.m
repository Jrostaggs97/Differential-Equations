function [A] = pendulum_jacobian(theta,h,m)
B = zeros(m); %initialize matrix

for i =1:m
   B(i,i) = -2 + (h^2)*cos(theta(i+1)); %diagonal entries
   
   if i<m
   B(i,i+1) = 1; %upper diagonal
   B(i+1,i) =1; %lower diagonal  
   end 
   
end

A = (1/h^2)*B;

end