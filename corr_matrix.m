function [ corr_val ] = corr_matrix( X1,X2 )
% this function is used to culculate the correlation matrix between two
% matrices

if(size(X1,1)>=size(X1,2))
    X1 = X1';
end
if(size(X2,1)>=size(X2,2))
    X2 = X2';
end

for i = 1:size(X1,1)
    for j = 1:size(X2,1)
        temp = corrcoef(X1(i,:),X2(j,:));
        corr_val(i,j) = temp(1,2);
    end
end

end

