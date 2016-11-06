function [Q, W] = calc_Q(X, A, epsilon,C )

[p,T] = size(X);
if nargin<4
    n = p;
else
    n = size(C,1);
end

if ~isequal(n,p)
    X = pinv(C)*X;
end
Xdiff = zeros(size(X));

for i = 1:p
    % filtering dimension does not matter
% Xdiff(i,:) = sqrt( (filter([1 -full(A(i,:))],1,X(i,:),[],2)).^2+epsilon^2);
Xdiff(i,:) = sqrt(filter([1 -full(A(i,:))],1,X(i,:)).^2+epsilon^2);
end
S = sparse(repmat(1:p,1,T),1:p*T,Xdiff(1:end),p,p*T );
Q = reshape(full(S),p,p,T);

XdiffInv = 1./Xdiff;

WS = sparse(repmat(1:p,1,T),1:p*T,XdiffInv(1:end),p,p*T );
W = reshape(full(WS),p,p,T);

% 
% WS0 = vec(XdiffInv)';
% WS1 = repmat(WS0,p,1);
% W = reshape(WS1,p,p,T);

end

