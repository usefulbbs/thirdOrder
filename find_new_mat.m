function [A2_new,seq_perm] = find_new_mat(A1,A2)
%%% this function is used to solve the permutation problem
%%% A1 is the real mixing matrix and the A2 is the estimation, the jjj is
%%% a number to show whether the estimate is right or not (the sign).
num_com = size(A1,2);
all_perm = perms([1:num_com]);
prod_all = zeros(1,size(all_perm,1));

for k = 1:size(all_perm,1)
    each_perm = all_perm(k,:);
    temp = abs(A1./(A2(:,each_perm)));
    temp = mat2gre1(temp);
    prod_all(k) = prod(prod(temp,2));
end
[~,min_ind] = min(prod_all);

seq_perm = all_perm(min_ind,:); % select the one which has the minimum A1./A2
A2_new = A2(:,seq_perm);
ratio_temp = A1./A2_new;
for k = 1:size(A2_new,2)
    each_col = sum(sign(ratio_temp(:,k)));
    diff_each_1 = abs(sum(A1(:,k)- A2_new(:,k)));
    diff_each_2 = abs(sum(A1(:,k) + A2_new(:,k)));
    if(each_col == size(A2_new,1)) 
        continue; % this condition is added by Liang on Aug.10, 2016
    elseif(each_col == -size(A2_new,1))
        A2_new(:,k) = -A2_new(:,k);
    elseif(diff_each_1>diff_each_2)
        A2_new(:,k) = -A2_new(:,k);
%         disp(j)
    end
end
end

function A = mat2gre1(A) %% to let each element in a positive matrix greater than 1
for i = 1:size(A,1)
    for j = 1:size(A,2)
        if(A(i,j)<1)
            A(i,j) = 1/A(i,j);
        end
    end
end
end










