function [A] = bndry_lyr_jacobian(u,epsilon,m,h)

B = zeros(m+2,m+2);
    
    for i = 2:m+1
        
     B(i,i) = -(2*epsilon)/(h^2) + ((u(i+1)-u(i-1))/(2*h) -1); %diagonal
     
         if i<m+1
             B(i,i-1)= epsilon/(h^2) -u(i)/(2*h);%lower diagonal

             B(i,i+1)= epsilon/(h^2) + u(i)/(2*h);%upper diagonal
        end
    end

A=B(2:m+1,2:m+1); %chop off 
    
end 