function [A_new] = single_three_datas(Xa,Aa,K,lag)
[m,n] = size(Aa);
[Xa Wa] = center_whitening(Xa); %Whitening the signal
 %Generate delay matrix
 H=2; 
 M=size(Xa,1);       %H=2£¬Third-order statistics
 tau=zeros(K,H);
 for i=1:K
 tau(i,:)=(lag*(i-1)+1);
 end 
 C22=zeros(M,M,M,K);
 
 for k=1:K        
 X1=circshift(Xa,0,2);
 X2=circshift(Xa,tau(k,1),2);%Continuously delay the observation signal
 X3=circshift(Xa,tau(k,2),2);
 C22(:,:,:,k) = cum3_zp(X1',X2',X3',0);%Calculate the third-order tensor, the size of the tensor is M*M*M*K
 end
 
[U,output] = cpd(C22,n);   %Do tensor decomposition
b1 = U{1};b2 = U{2};b3 = U{3};
%---------------------------------------------------------------------------
%Perform SVD decomposition on the result of tensor decomposition to obtain the mixed matrix A
for i = 1:size(Aa,2)
    b = [b1(:,i),b2(:,i),b3(:,i)];
    [UU S V] = svd(b);
    A_new(:,i) = UU(:,1);
end
 A_new = inv(Wa)*A_new;
 %Normalize the newly updated mixing matrix
 for i = 1:size(A_new,2)
    A_new(:,i) = A_new(:,i)/sqrt(sum(A_new(:,i).^2));
 end
end

