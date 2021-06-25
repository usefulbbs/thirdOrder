
clear all;close all;

[s1,fs]= audioread('S013330161.wav');
s2 = audioread('S213031143.wav');
s3 = audioread('S013330147.wav');
s4 = audioread('S213030267.wav');

sample_num = 10000;
S = zeros(4, sample_num);
S(1,:) = s1(14001:14000+sample_num)';
S(2,:) = s2(14001:14000+sample_num)';
S(3,:) = s3(14001:14000+sample_num)';
S(4,:) = s4(14001:14000+sample_num)';


figure;
axis off;
subplot(4,1,1);
plot(S(1,:));
set(gca,'xtick',-inf:inf:inf);
set(gca,'ytick',-inf:inf:inf);
subplot(4,1,2);
plot(S(2,:));
set(gca,'xtick',-inf:inf:inf);
set(gca,'ytick',-inf:inf:inf);
subplot(4,1,3);
plot(S(3,:));
set(gca,'xtick',-inf:inf:inf);
set(gca,'ytick',-inf:inf:inf);
subplot(4,1,4);
plot(S(4,:));
set(gca,'xtick',-inf:inf:inf);
set(gca,'ytick',-inf:inf:inf);

A = [0.2828 0.1701 0.5517 0.3245
     0.5896 0.6763 0.9674 0.5210
     0.3720 0.5696 0.2284 0.8649];

X = A * S(1:4,1:sample_num);

source_num = 4;


figure;
axis off;
subplot(3,1,1);
plot(X(1,:));
set(gca,'xtick',-inf:inf:inf);
set(gca,'ytick',-inf:inf:inf);
subplot(3,1,2);
plot(X(2,:));
set(gca,'xtick',-inf:inf:inf);
set(gca,'ytick',-inf:inf:inf);
subplot(3,1,3);
plot(X(3,:));
set(gca,'xtick',-inf:inf:inf);
set(gca,'ytick',-inf:inf:inf);
t0=clock;
wlen=256;
timestep=128;
numfreq=256;

awin=hamming(wlen);%analysis window is a Hammingwindow
swin =awin;
XUSED = X(:,1:sample_num);
tfmat1 = tfanalysis(XUSED(1,:)',awin,timestep, numfreq);
tfmat2 = tfanalysis(XUSED(2,:)',awin,timestep, numfreq);
tfmat3 = tfanalysis(XUSED(3,:)',awin,timestep, numfreq);


rtfmat1 = reshape(tfmat1,1,[]);
rtfmat2 = reshape(tfmat2,1,[]);
rtfmat3 = reshape(tfmat3,1,[]);

Xtfmat = [rtfmat1; rtfmat2; rtfmat3];


tfmat_norm = sum(abs(Xtfmat).^2, 1);

norm_threshold = mean(tfmat_norm)/100;

tfmat = Xtfmat(:,find(tfmat_norm> norm_threshold));
tfmat_norm = tfmat_norm(1,find(tfmat_norm> norm_threshold));

tfmat = abs(tfmat).*sign(real(tfmat));

XW = tfmat./repmat(sum(tfmat.^2).^(1/2),3,1).*repmat(sign(tfmat(1,:)),3,1);

dist = abs(XW'*XW);
threshold = 0.9999;
groups_num = size(XW,1) +1;
count = sum(dist > threshold, 2);


eta = 60; % It is determined by viewing the results of the plot of OMEGA, and is set as nuber which leads to source_num groups.
OMEGA = XW(:, count> eta);


MAXiter = 100; % Maximum iteration for KMeans Algorithm
REPlic = 100; % Replication for KMeans Algorithm
[IDX,C] = kmeans(OMEGA',source_num,'start','sample','maxiter',MAXiter,'replicates',REPlic);
C = C';

disp(['The time cost for the mixing matrix estimation = ',num2str(etime(clock,t0)) 'seconds']);

AA = A./repmat(sum(A.^2).^(1/2),size(A,1),1);

EA = zeros(3,4);

for i = 1:source_num
    [u1,u2]=max(abs(C(:,i)'*AA));
    EA(:,u2) = C(:,i);
end

error = norm(AA - EA)

t0=clock;
% Recover the time-frequency representations of the source signals

ES = source_recovery2(Xtfmat,EA);
% Reshape the estimated signals as the time-frequency matrix form for each
% source.
estfmat1 = reshape(ES(1,:),numfreq,[]);
estfmat2 = reshape(ES(2,:),numfreq,[]);
estfmat3 = reshape(ES(3,:),numfreq,[]);
estfmat4 = reshape(ES(4,:),numfreq,[]);


% Apply the inverse sort time-frequency transform to the estimated time-frequency
%representations of the sources.
est1 = tfsynthesis(estfmat1,sqrt(2)*awin/256, timestep, numfreq);
est2 = tfsynthesis(estfmat2,sqrt(2)*awin/256, timestep, numfreq);
est3 = tfsynthesis(estfmat3,sqrt(2)*awin/256, timestep, numfreq);
est4 = tfsynthesis(estfmat4,sqrt(2)*awin/256, timestep, numfreq);

% Rescale the output recovered signals to between 0 and 1;
E_source = zeros(4, sample_num);
E_source(1,:) = est1(1:sample_num)./max(abs(est1(1:sample_num)));
E_source(2,:) = est2(1:sample_num)./max(abs(est2(1:sample_num)));
E_source(3,:) = est3(1:sample_num)./max(abs(est3(1:sample_num)));
E_source(4,:) = est4(1:sample_num)./max(abs(est4(1:sample_num)));
ES = E_source;

disp(['The time cost for source recovery = ',num2str(etime(clock,t0)) 'seconds']);

delta2 = zeros(size(S,1),1);
for i = 1:4
    delta2(i,1) = (S(i,:)*E_source(i,:)')/(E_source(i,:)*E_source(i,:)');
    E_source(i,:) = E_source(i,:).*delta2(i,1);
end

figure;
subplot('position',[0.025 0.1 0.2 0.8]) ;
% subplot('position',[left bottom width height]) 
plot(S(1,:));
set(gca,'fontsize',16, 'GridLineStyle', '-.','LineWidth', 1.5);
subplot('position',[0.275 0.1 0.2 0.8]) ;
plot(S(2,:));
set(gca,'fontsize',16, 'GridLineStyle', '-.','LineWidth', 1.5);
subplot('position',[0.525 0.1 0.2 0.8]) ;
plot(S(3,:));set(gca,'fontsize',16, 'GridLineStyle', '-.','LineWidth', 1.5);
subplot('position',[0.775 0.1 0.2 0.8]) ;
plot(S(4,:));set(gca,'fontsize',16, 'GridLineStyle', '-.','LineWidth', 1.5);

% Plot mixed signals
figure;
subplot('position',[0.15 0.1 0.2 0.8]) ;
plot(X(1,:));
set(gca,'fontsize',16, 'GridLineStyle', '-.','LineWidth', 1.5);
subplot('position',[0.4 0.1 0.2 0.8]) ;
plot(X(2,:));
set(gca,'fontsize',16, 'GridLineStyle', '-.','LineWidth', 1.5);
subplot('position',[0.65 0.1 0.2 0.8]) ;
plot(X(3,:));
set(gca,'fontsize',16, 'GridLineStyle', '-.','LineWidth', 1.5);

figure;
subplot(4,1,1);
plot(E_source(1,:));

set(gca,'fontsize',16, 'GridLineStyle', '-.','LineWidth', 1.5);
subplot(4,1,2);
plot(E_source(2,:));

set(gca,'fontsize',16, 'GridLineStyle', '-.','LineWidth', 1.5);
subplot(4,1,3);
plot(E_source(3,:));

set(gca,'fontsize',16, 'GridLineStyle', '-.','LineWidth', 1.5);
subplot(4,1,4);
plot(E_source(4,:));
set(gca,'fontsize',16, 'GridLineStyle', '-.','LineWidth', 1.5);


figure;
axis off;
subplot(4,1,1);
plot(E_source(1,:));
set(gca,'xtick',-inf:inf:inf);
set(gca,'ytick',-inf:inf:inf);
subplot(4,1,2);
plot(E_source(2,:));
set(gca,'xtick',-inf:inf:inf);
set(gca,'ytick',-inf:inf:inf);
subplot(4,1,3);
plot(E_source(3,:));
set(gca,'xtick',-inf:inf:inf);
set(gca,'ytick',-inf:inf:inf);
subplot(4,1,4);
plot(E_source(4,:));
set(gca,'xtick',-inf:inf:inf);
set(gca,'ytick',-inf:inf:inf);

error2 = norm(E_source - S)
