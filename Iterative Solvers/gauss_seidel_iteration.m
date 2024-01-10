%Gauss Siedel iteration for 2d finite difference approx
function [u_approx,residuals,tol_iter] = gauss_seidel_iteration(u_init,b,A,max_iter,tol)
s = size(A);
N =s(1);
r = b-A*u_init;
L = tril(A);
M = L;
z = M\r;
u=u_init;
residuals = zeros(1,max_iter);
    for iter = 1:max_iter

        u = u + z;
        r = b-A*u;
        z=M\r;
        rel_resid = norm(r)/norm(b);
        residuals(iter) = rel_resid;
        if rel_resid <tol
            tol_iter = iter;
            break
        else
            
            tol_iter = max_iter;
            continue
        end
    end

    u_approx=u;

end