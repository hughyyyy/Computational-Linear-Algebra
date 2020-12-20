function [x,iter] = CG(A,b,x_initial,maxiter,tol) 
    x0 = x_initial;
    r = b - A*x;
    y = -r;
    z = A*y;
    s = y'*z;
    t = (r'*y)/s;
    x = x + t*y;
    iter=1;
    for k = 1:maxiter
       r = r - t*z;
       if(norm(r) <= tol*norm(b))
            return;
       end
       B = (r'*z)/s;
       y = -r + B*y;
       z = A*y;
       s = y'*z;
       t = (r'*y)/s;
       x = x + t*y;
       iter=iter+1;
    end
 end