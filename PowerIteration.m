function [v,lambda,iter] = PowerIteration(A,v0,maxiter,tol)
v = v0;
for iter = 1:maxiter
    w = A * v;
    v = w/norm(w);
    lambda = transpose(v)* A * v;
    if norm(A*v-lambda*v) < tol;
        break
    end
end
end
