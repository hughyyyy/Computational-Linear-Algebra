function [x,iter] = GS(A,b,x_initial,maxiter,tol)
x0 = x_initial;
iter = 0;
D = diag(diag(A));
U = triu(A-D);
Inv = inv(D+D);
r=b-A*x0;
      while (norm(r)>tol*norm(b) || iter<=maxiter)
          x1 = x0;
          x0 = Inv*(b-(U*x1));
          i = i + 1;
      end
x = x0;
end