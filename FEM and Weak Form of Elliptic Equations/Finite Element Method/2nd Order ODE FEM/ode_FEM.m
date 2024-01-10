clc;clear;
m=9; %for nicer looking grid spacing on uniform grid do 10^n - 1 

lbc = 0; %left boundary condition
rbc=0; %right boundary condition
grid = linspace(0,1,m+2); %set up gride
%uniform grid - comment out nonuniform for loop
x=grid(1:m+2);

%non uniform gird
for i = 1:length(x)
    
    x(i) = (i./(length(x))).^2;
    
end

h = x(2:m+2)-x(1:m+1); %grid spacing


%build FEM matrix
A = zeros(m,m);
    for i=1:m
    
    A(i,i) = (1/h(i)^2)*(cubic(x(i+1))-cubic(x(i)))+(1/h(i+1)^2)*(cubic(x(i+2))-cubic(x(i+1)));
    
        if i<m
            A(i,i+1) = -(1./h(i+1).^2).*(cubic(x(i+2))-cubic(x(i+1))); 
            A(i+1,i) = A(i,i+1); %matrix is symmetric
        end
    
    end 


h_mid = (x(2:m+2) - x(1:m+1))./2; %midpoint spacing
x_midpoint =x(2:m+2)-h_mid; %midpoint values

f_midpoint = rhs_func(x_midpoint); %right hand side function at midpoint
int_f = f_midpoint(1:m).*h_mid(1:m) + f_midpoint(2:m+1).*h_mid(2:m+1); %midpoint approximation for integral of f phi


u =A\int_f'; %solve system
u = [lbc; u; rbc]; 
soln = x.*(1-x); %real solution

maxh = max(h);
inf_norm = max(abs(u-soln')); %infinity norm
two_norm = sqrt(h(1))*norm(u-soln'); %two norm


