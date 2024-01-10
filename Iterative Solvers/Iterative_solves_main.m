n = 6*25;
h = 1/n;
N = (n-1)^2;

%  Form block tridiagonal finite difference matrix A and right-hand side 
%  vector b.

A = zeros(N,N);
b = ones(N,1);         % Use right-hand side vector of all 1's.

%  Loop over grid points in y direction.

for j=1:n-1,
  yj = j*h;
  yjph = yj+h/2;  yjmh = yj-h/2;

%    Loop over grid points in x direction.

  for i=1:n-1
    xi = i*h;
    xiph = xi+h/2;  ximh = xi-h/2;
    aiphj = 1 + xiph^2 + yj^2;
    aimhj = 1 + ximh^2 + yj^2;
    aijph = 1 + xi^2 + yjph^2;
    aijmh = 1 + xi^2 + yjmh^2;
    k = (j-1)*(n-1) + i;
    A(k,k) = aiphj+aimhj+aijph+aijmh;
    if i > 1, A(k,k-1) = -aimhj; end;
    if i < n-1, A(k,k+1) = -aiphj; end;
    if j > 1, A(k,k-(n-1)) = -aijmh; end;
    if j < n-1, A(k,k+(n-1)) = -aijph; end;
  end;
end;
A = (1/h^2)*A;

%Jacobi Iteration
u_init = zeros(N,1);
max_iter = 8000;
tol = 10^(-5);
% [jacobi_approx_u,jacobi_residuals,jacobi_iter_h] = jacobi_iteration(u_init,b,A,max_iter,tol);
% jacobi_iter_list(iter_h) = jacobi_iter_h

% Gauss Seidel Iteration
% [gs_approx_u,gs_residuals,gs_iter_h] = gauss_seidel_iteration(u_init,b,A,max_iter,tol);
% gs_iter_list(iter_h) = gs_iter_h
% SOR Iteration
% SOR_final_resid = zeros(1,10);
% 
% for omega = [1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
%     i = omega*10 -10;

omega= 1.9;
[SOR_approx_u,SOR_residuals,SOR_iter_h] = SOR_iteration(u_init,b,A,omega,max_iter,tol);
sor_iter_list(6) = SOR_iter_h
%     SOR_final_resid(i) = min(SOR_residuals(SOR_residuals>0));
% end

%Conjugate Gradient Method
%[cg_approx_u,flag,relres,cg_iter_h,cg_residuals] = pcg(A,b,10^(-5),4000);
%cg_iter_list(iter_h) = cg_iter_h;

% plot(iter_list,log10(residual_list))
% xlabel("$log_{10}$(iteration)")
% ylabel("$log_{10}$(residual error)")
% title("Residual error vs.iteration")


%Preconditioned Conjugate Gradient Method
%A = sparse(A);
%L = ichol(A);    
%[pcg_approx_u,flag,relres,pcg_iter_h,pcg_residuals] = pcg(A,b,10^(-5),500,L,L'); %this one smacks - can get really any desired tol within 150 iterations
%pcg_iter_list(iter_h) = pcg_iter_h;
%end