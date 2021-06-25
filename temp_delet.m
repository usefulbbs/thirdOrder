
% for snr = 50:-5:5          % Simulations under various SNRs
%     snr
%     SNR = num2str(snr);
clc
ep = 1/30;            % Epsilon (SSD)
n = 4;                 % The number of sources
m = 3;                 % The number of mixtures
L = 10000;             % The length of sources

numch=4;
% Load mixing matrix
err = zeros(500,3);
fs=250;

%     s1 = wavread('speech1.wav')';
%     s2 = wavread('speech2.wav')';
%     s3 = wavread('speech3.wav')';
%     s4 = wavread('speech4.wav')';
% 
%     % Energy normalization
%     s1 = s1(1:L);
%     s2 = s2(1:L);
%     s3 = s3(1:L);
%     s4 = s4(1:L);
%     % Make the variance one
%     E_s1 = sum(s1.^2);
%     E_s2 = sum(s2.^2);
%     E_s3 = sum(s3.^2);
%     E_s4 = sum(s4.^2);
%     s1 = sqrt(L/E_s1)*s1;
%     s2 = sqrt(L/E_s2)*s2;
%     s3 = sqrt(L/E_s3)*s3;
%     s4 = sqrt(L/E_s4)*s4;           
% 
%     % Source matrix
%     S = [s1; s2; s3; s4];

    
for i = 1:300
    d = fdesign.bandpass('N,F3dB1,F3dB2',30,20,60,fs);
    Hd = design(d,'butter');
    EEGsim=[];
    [EEGsim] = eeg_simu_para(numch,fs,L);
    % EMG generation
    temp = randn(1,L+200);  % can assume any fs cause white noise's PSD is plain; fs depend on the sampling rate of bp
%     wn = zeros(1,eeg_length);
    tt = filter(Hd,temp);   % band pass filtering 20-60 Hz
    wn = tt(100:(L+99));
    EMGsim = wn(1,:);
    S = [EEGsim;EMGsim];

    
    disp(i);
    A = rand(4,5)*2-1;
    A = scale_mixingmatrix_para(A);
    
    % Mixtures
    X = A*S;

    EA1 = ubss_mix4_s5(X);
    EA1=scale_mixingmatrix_para(EA1);
    [EA1] = find_new_mat(A,EA1);
    err(i,1) = norm(A-EA1,'fro')/norm(A,'fro');
    
%     N =1024;
%     
%     % DFT of mixture
%     for k = 1:m
%         Xf{k} = analyze(X(k,:),N);
%     end
%     
%     % Size of Xf
%     [F,T] = size(Xf{1});
%     
%     % First, Autosource point selection
%     for t = 1:T
%         for f = 2:F/2+1
%             if norm([real(Xf{1}(f,t)) real(Xf{2}(f,t)) real(Xf{3}(f,t))]) < 0.001 | norm([imag(Xf{1}(f,t)) imag(Xf{2}(f,t)) imag(Xf{3}(f,t))]) < 0.001
%                 Xf{1}(f,t) = 0;
%                 Xf{2}(f,t) = 0;
%                 Xf{3}(f,t) = 0;
%             end
%         end
%     end
%     
%     % Second, SSS
%     SS1 = [];
%     SS2 = [];
%     SS3 = [];
%     for t = 1:T
%         for f = 2:F/2+1
%             if Xf{1}(f,t) ~= 0
%                 %                if (abs(imag(Xf{2}(f,t)/Xf{1}(f,t))) < ep) & (abs(imag(Xf{3}(f,t)/Xf{1}(f,t))) < ep) & (abs(imag(Xf{3}(f,t)/Xf{2}(f,t))) < ep)
%                 if (abs(angle(Xf{2}(f,t)/Xf{1}(f,t))) < ep || abs(angle(Xf{2}(f,t)/Xf{1}(f,t))) > pi-ep) & (abs(angle(Xf{3}(f,t)/Xf{1}(f,t))) < ep || abs(angle(Xf{3}(f,t)/Xf{1}(f,t))) > pi-ep) & (abs(angle(Xf{3}(f,t)/Xf{2}(f,t))) < ep || abs(angle(Xf{3}(f,t)/Xf{2}(f,t))) > pi-ep)
%                     SS1 = [SS1 Xf{1}(f,t)];
%                     SS2 = [SS2 Xf{2}(f,t)];
%                     SS3 = [SS3 Xf{3}(f,t)];
%                 end
%             end
%         end
%     end
%     
%     SS = [SS1; SS2; SS3];
%     
%     % ==================
%     % Estimation of error
%     % ==================
%     SSa = [real(SS) imag(SS)];
%     
%     for nn = 1:length(SSa(1,:))
%         SSa1(:,nn) = SSa(:,nn)/SSa(1,nn);
%         SSa1(:,nn) = SSa1(:,nn)/norm(SSa1(:,nn),2);
%     end
%     
%     
%     % Kmeans clustering
%     [IDX,C] = kmeans(SSa1',n,'emptyaction','singleton');
%     %     A(:,4) = -A(:,4);
%     %     [IDX,C] = kmeans(SSa1',4,'start',A');         % Performance of Kmeans depends on initialization.
%     % You can use another
%     % more efficient method
%     % for clustering.
%     % I recommend it.
%     
%     % Find the points in each cluster
%     C1 = [];
%     C2 = [];
%     C3 = [];
%     C4 = [];
%     for ii = 1:length(IDX)
%         if IDX(ii) == 1
%             C1 = [C1 SSa1(:,ii)];
%         elseif IDX(ii) == 2;
%             C2 = [C2 SSa1(:,ii)];
%         elseif IDX(ii) == 3
%             C3 = [C3 SSa1(:,ii)];
%         else
%             C4 = [C4 SSa1(:,ii)];
%         end
%     end
%     
%     % Centroid of each cluster (normalized)
%     cent1 = mean(C1'); cent1 = cent1/norm(cent1,2);
%     cent2 = mean(C2'); cent2 = cent2/norm(cent2,2);
%     cent3 = mean(C3'); cent3 = cent3/norm(cent3,2);
%     cent4 = mean(C4'); cent4 = cent4/norm(cent4,2);
%     
%     W = [cent1' cent2' cent3' cent4'];
%     W=scale_mixingmatrix_para(W);
%     [W]=find_new_mat(A,W);
%     err(i,2) = norm(A-W,'fro')/norm(A,'fro');
% 
% % %     
% %     
% %     %%
% %     [Aa_2,~] = est_one_matrix(X,A,10,1); % based on single dataset,t_1{k,m}(n,j) corresponding to single set
% %     [Aa_2] = find_new_mat(A,Aa_2);
% %     err(i,3) = norm(A-Aa_2,'fro')/norm(A,'fro');
% %     err(i,1:3)
%      Aa_2= rand(3,4)*2-1;
%      Aa_2 = scale_mixingmatrix_para(Aa_2);
%     [Aa_2] = find_new_mat(A,Aa_2);
%     err(i,3) = norm(A-Aa_2,'fro')/norm(A,'fro');
%     err(i,1:3)
end