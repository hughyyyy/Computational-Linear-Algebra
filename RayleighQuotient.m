function [v,lambda,iter] = RayleighQuotient(A,v0,maxiter,tol)
[m,n] = size(A);
v = v0;
lambda = transpose(v)*A*v;
for iter = 1:maxiter
    w = (A-lambda*eye(n))\v;
    v = w/norm(w);
    lambda = transpose(v)*A*v;
    if norm(A*v-lambda*v) < tol;
        break
    end
end
end