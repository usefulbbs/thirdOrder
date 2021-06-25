function [ Y, W ] = center_whitening( X )
%center_whitening is used for center and whitening the signal
%   X is the Original signal
%   Y is the preprocessed signal and the W is the whitening matrix
% centering
X = X';
[n,v] = size(X);
X = bsxfun(@minus,X,mean(X));
% whitening
[U,S,~] = svd(X,'econ');
Y = U*(S*pinv(S))*sqrt(n-1);
W = Y'*pinv(X');
Y = Y';
end