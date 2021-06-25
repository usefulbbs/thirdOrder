clc;
clear all;
%————————————————————————————————————————
%Initial parameter setting
%The algorithm for signal recovery here is literature：
%Kim S G, Yoo C D. Underdetermined blind source separation based on subspace representation[J]. 
%IEEE Transactions on Signal processing, 2009, 57(7): 2604-2614.
delay_num=10;        %The number of signal delays
step_size=2;         %Lengthen the step
addpath(genpath('三阶统计量代码上传'))
addpath(genpath('tensorlab_2014-05-07'))
%---------------------------------------------------------------------------
%Signal reading
[s1, fs1] = audioread('speech1.wav');
[s2, fs2] = audioread('speech2.wav');
[s3, fs3] = audioread('speech3.wav');
[s4, fs4] = audioread('speech4.wav');
Sa = [s1 s2 s3 s4]';
%---------------------------------------------------------------------------
load ori_mixing3x4
Xa = A*Sa;
%Use third-order statistics to calculate the mixed matrix
Aa =single_three_datas(Xa,A,delay_num,step_size);
Aa = find_new_mat(A,Aa);
error = 10*log10(sum((Aa-A).^2)/sum(A.^2)); %Evaluate the mixing matrix A
%---------------------------------------------------------------------------
%After evaluating the mixing matrix, restore the source signal
Sb = estimate_s_3(Xa,Aa);
%---------------------------------------------------------------------------
%Evaluate the recovery signal
corr_mat = corr_matrix(Sa(:,1:size(Sb,2)),Sb);
PCC = 0;
for k = 1:size(corr_mat,1);
PCC = PCC + abs(corr_mat(k,k));   
end
PCC = PCC/(size(corr_mat,1));
%Draw the source signal
figure(1)
subplot(4,1,1)
plot(s1)
ylim([-0.6 0.6])
set(gca,'FontSize',9)
ylabel('源信号1','FontSize',13)
subplot(4,1,2)
plot(s2)
set(gca,'FontSize',9)
ylim([-0.6 0.6])
ylabel('源信号2','FontSize',13)
subplot(4,1,3)
plot(s3)
set(gca,'FontSize',9)
ylabel('源信号3','FontSize',13)
ylim([-0.6 0.6])
subplot(4,1,4)
plot(s4)
set(gca,'FontSize',9)
ylabel('源信号4','FontSize',13)
ylim([-0.6 0.6])
sgtitle('源信号','FontSize',14)

%Draw the mixed signal
figure(2)
subplot(3,1,1)
plot(Xa(1,:))
set(gca,'FontSize',9)
ylabel('观测信号1','FontSize',13)
subplot(3,1,2)
plot(Xa(2,:))
set(gca,'FontSize',9)
ylabel('观测信号2','FontSize',13)
subplot(3,1,3)
plot(Xa(3,:))
set(gca,'FontSize',9)
ylabel('观测信号3','FontSize',13)
sgtitle('混合信号','FontSize',14)

%Draw the recovered source signal
figure(3)
subplot(4,1,1)
plot(Sb(1,:))
ylim([-0.6 0.6])
set(gca,'FontSize',9)
ylabel('恢复信号1','FontSize',13)
subplot(4,1,2)
plot(Sb(2,:))
set(gca,'FontSize',9)
ylim([-0.6 0.6])
ylabel('恢复信号2','FontSize',13)
subplot(4,1,3)
plot(Sb(3,:))
set(gca,'FontSize',9)
ylim([-0.6 0.6])
ylabel('恢复信号3','FontSize',13)
subplot(4,1,4)
plot(Sb(4,:))
set(gca,'FontSize',9)
ylim([-0.6 0.6])
ylabel('恢复信号4','FontSize',13)
sgtitle('恢复的源信号','FontSize',14)

