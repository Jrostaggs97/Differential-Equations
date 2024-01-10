%Jacobi iteration for 2d finite difference approx
function [u_approx,rel_resid,tol_iter] = jacobi_iteration(u_init,b,A,max_iter,tol)
s = size(A);
N =s(1);
r = b-A*u_init;
I = eye(N,N);
d = diag(A);
M = full(spdiags(d,0,I));
z = M\r;
u=u_init;
resid = zeros(1,max_iter);
    for iter = 1:max_iter

        u = u + z;
        r = b-A*u;
        z=M\r;
        rel_resid = norm(r)/norm(b);
        resid(iter)= rel_resid;
        if rel_resid <tol
            tol_iter = iter;
            break
        else
            tol_iter = max_iter;
            continue

        end
    end

    u_approx=u;
    rel_resid = resid;
end