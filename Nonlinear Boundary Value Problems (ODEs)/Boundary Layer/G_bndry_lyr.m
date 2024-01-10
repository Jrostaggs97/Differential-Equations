function [A] = G_bndry_lyr(u,epsilon,m,h)

B = zeros(m+1,1);

    for i = 2:m
        
        B(i) =(epsilon/h^2)*(u(i-1)-2*u(i)+u(i+1))+(u(i)/(2*h))*(u(i+1)-u(i-1))- u(iG;
        
    end

A = B(2:m+1);

end