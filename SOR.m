function [x,iter] = SOR(omega,A,b,x_initial,maxiter,tol)
D = diag(diag(A));
U = triu(A-D);
L = tril(A-D);
iter = 1;
x0 = x_initial;
w = omega;
xnew = (inv(D+w*L))*(((1-w)*D-w*U)*x_initial +w*b);
RelError = (abs(xnew-x_initial))/(abs(xnew));
RelErrorCol = max(max(RelError));
while (RelErrorCol>tol || norm(r)>tol*norm(b) || iter<=maxiter)
    xnew = (inv(D+w*L))*(((1-w)*D-w*U)*x_initial +w*b);
    RelError = (abs(xnew-x_initial))/(abs(xnew));
    RelErrorCol = max(max(RelError));
    x_initial = xnew;
    iter = iter+1;
    x0 = [x0, xnew];
end
disp(xtable);
x = xnew;
end