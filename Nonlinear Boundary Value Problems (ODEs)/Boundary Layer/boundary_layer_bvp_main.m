m = 159;
epsilon = .01;
x= linspace(0,1,m+2);
h = 1/(m+1);
a=0;
b=1;
alpha =-1;
beta=1.5;

x_bar = (a+b-alpha-beta)/2; 

w_bar = (a-b+beta-alpha)/2;

u = x-x_bar+w_bar*tanh(w_bar*(x-x_bar)/(2*epsilon)); %initial guess from book via asymptotics 


for k =1:50
   
    G_vec = G_bndry_lyr(u,epsilon,m,h);

    jacob = bndry_lyr_jacobian(u,epsilon,m,h);
    
    del = jacob\G_vec;
    
    u(2:m+1) = u(2:m+1)' - jacob\G_vec;
    
    if max(abs(jacob\G_vec))<10^(-14)
        break
    end

end 

plot(x,u)
xlabel("x")
ylabel("u(x), h=1/160")