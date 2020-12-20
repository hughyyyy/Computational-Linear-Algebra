function [x, iter] = Jacobi(A, b, x_initial, tol, maxiter)
D = diag(diag(A));
R = A - D;
iter = 1;
x0 = x_initial;
r = b-A*x;
while (norm(r)>tol*norm(b) || iter <= maxiter)
    x = x0;    
    x0 = inv(D)*(b - R*x);  
    iter = iter+1
end
