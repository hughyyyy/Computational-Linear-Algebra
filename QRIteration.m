function [V,Lambda,iter] = QRIteration(A,maxiter,tol)
Ai = A;
[m,n] = size(A);
V = eye(n);
for iter = 1:maxiter
    tolcount=0;
    [Q,R] = qr(Ai);
    V = V * Q;
    Ai = R*Q;
    Lambda = diag(Ai);
    for i = 1:n
        if norm(A*V(:,i)-Lambda(i,1)*V(:,i)) < tol
            tolcount = tolcount + 1;
        end
    end
    if tolcount == n
        return
    end
end
end
