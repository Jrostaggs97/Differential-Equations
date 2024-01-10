%Jacobi iteration for 2d finite difference approx
function [u_approx,residuals,tol_iter] = SOR_iteration(u_init,b,A,omega,max_iter,tol)
s = size(A);
N =s(1);
I = eye(N,N);
d = diag(A);
L = tril(A,-1);
D = spdiags(d,0,I);
M = (1/omega)*D + L;
r = b-A*u_init;
z = M\r;
u=u_init;
residuals =zeros(1,max_iter);
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
            tol_iter= max_iter;
            continue
        end
    end

    u_approx=u;

end