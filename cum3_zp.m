function [c3] = cum3_zp(X1,X2,X3,prewhiten)
%CUM4 Fourth-order cumulant tensor.¸Ä±à×Ôcum4()

if nargin < 2, prewhiten = false; end
if ischar(prewhiten), prewhiten = strcmpi(prewhiten,'prewhiten'); end

% Center the variables.
X1 = bsxfun(@minus,X1,mean(X1,1));
X2 = bsxfun(@minus,X2,mean(X2,1));
X3 = bsxfun(@minus,X3,mean(X3,1));
% Apply a prewhitening to X if requested.
n = size(X1,1);

if prewhiten
    [U,S,~] = svd(X1,'econ');
    X1 = U*(S*pinv(S))*sqrt(n);
    [U,S,~] = svd(X2,'econ');
    X2 = U*(S*pinv(S))*sqrt(n);
    [U,S,~] = svd(X3,'econ');
    X3 = U*(S*pinv(S))*sqrt(n);
end
% Compute c3 = E[xi*conj(xj)*conj(xk)].
c3 = mean(bsxfun(@times,bsxfun(@times,permute(X1,[2 3 4 1]), ...
    permute(conj(X2),[3 2 4 1])),permute(conj(X3),[3 4 2 1])),4);
end
