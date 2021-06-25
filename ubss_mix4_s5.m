function [A,fail,t] = ubss_mix4_s5(X)
%% The value of ep is quite important, for the simulation data, we select it as 0.45;
% Underdetermined Blind Source Separation Based on Subspace Representation
% Author: SangGyun Kim
%         MMP, EECS, KAIST
% Please, DO NOT share this code with anyone outside of the lab.


tic;
% ep = 1/30;            % Epsilon (SSD)
ep = 0.45;            % Epsilon (SSD)
n = 5;                 % The number of sources
m = 4;                 % The number of mixtures
L = size(X,2);             % The length of sources

N =256;
fail=0;
% DFT of mixture
for k = 1:m
    Xf{k} = analyze(X(k,:),N);
end

% Size of Xf
[F,T] = size(Xf{1});

% First, Autosource point selection
for t = 1:T
    for f = 2:F/2+1
        if (norm([real(Xf{1}(f,t)) real(Xf{2}(f,t)) real(Xf{3}(f,t)) real(Xf{4}(f,t))]) < 0.001 || norm([imag(Xf{1}(f,t)) imag(Xf{2}(f,t)) imag(Xf{3}(f,t)) imag(Xf{4}(f,t))]) < 0.001)
            Xf{1}(f,t) = 0;
            Xf{2}(f,t) = 0;
            Xf{3}(f,t) = 0;
            Xf{4}(f,t) = 0;
        end
    end
end

% Second, SSS
SS1 = [];
SS2 = [];
SS3 = [];
SS4 = [];
for t = 1:T
    for f = 2:F/2+1
        if Xf{1}(f,t) ~= 0
            %                if (abs(imag(Xf{2}(f,t)/Xf{1}(f,t))) < ep) & (abs(imag(Xf{3}(f,t)/Xf{1}(f,t))) < ep) & (abs(imag(Xf{3}(f,t)/Xf{2}(f,t))) < ep)
            if ((abs(angle(Xf{2}(f,t)/Xf{1}(f,t))) < ep || abs(angle(Xf{2}(f,t)/Xf{1}(f,t))) > pi-ep) && (abs(angle(Xf{3}(f,t)/Xf{1}(f,t))) < ep || abs(angle(Xf{3}(f,t)/Xf{1}(f,t))) > pi-ep) && (abs(angle(Xf{3}(f,t)/Xf{2}(f,t))) < ep || abs(angle(Xf{3}(f,t)/Xf{2}(f,t))) > pi-ep)...
                    &&(abs(angle(Xf{4}(f,t)/Xf{1}(f,t))) < ep || abs(angle(Xf{4}(f,t)/Xf{1}(f,t))) > pi-ep) && (abs(angle(Xf{4}(f,t)/Xf{2}(f,t))) < ep || abs(angle(Xf{4}(f,t)/Xf{2}(f,t))) > pi-ep) && (abs(angle(Xf{4}(f,t)/Xf{3}(f,t))) < ep || abs(angle(Xf{4}(f,t)/Xf{3}(f,t))) > pi-ep))
                SS1 = [SS1 Xf{1}(f,t)];
                SS2 = [SS2 Xf{2}(f,t)];
                SS3 = [SS3 Xf{3}(f,t)];
                SS4 = [SS4 Xf{4}(f,t)];
            end
        end
    end
end

SS = [SS1; SS2; SS3; SS4];
if(size(SS,1)<4)
    fail =1;
    A =[];
    return;
end
% ==================
% Estimation of error
% ==================
SSa = [real(SS) imag(SS)];

for nn = 1:length(SSa(1,:))
    SSa1(:,nn) = SSa(:,nn)/SSa(1,nn);
    SSa1(:,nn) = SSa1(:,nn)/norm(SSa1(:,nn),2);
end
if(size(SSa1,2)<n)
    SSa1,
    fail=1;
    A=[];
    return;
end

% Kmeans clustering
[IDX,C] = kmeans(SSa1',n,'emptyaction','singleton');
%     A(:,4) = -A(:,4);
%     [IDX,C] = kmeans(SSa1',4,'start',A');         % Performance of Kmeans depends on initialization.
% You can use another
% more efficient method
% for clustering.
% I recommend it.

% Find the points in each cluster
C1 = [];
C2 = [];
C3 = [];
C4 = [];
C5 = [];
for ii = 1:length(IDX)
    if IDX(ii) == 1
        C1 = [C1 SSa1(:,ii)];
    elseif IDX(ii) == 2;
        C2 = [C2 SSa1(:,ii)];
    elseif IDX(ii) == 3;
        C3 = [C3 SSa1(:,ii)];
    elseif IDX(ii) == 4;
        C4 = [C4 SSa1(:,ii)];
    else
        C5 = [C5 SSa1(:,ii)];
    end
end

% Centroid of each cluster (normalized)
cent1 = mean(C1'); cent1 = cent1/norm(cent1,2);
cent2 = mean(C2'); cent2 = cent2/norm(cent2,2);
cent3 = mean(C3'); cent3 = cent3/norm(cent3,2);
cent4 = mean(C4'); cent4 = cent4/norm(cent4,2);
cent5 = mean(C5'); cent5 = cent5/norm(cent5,2);
if((length(cent1)+length(cent2)+length(cent3)+length(cent4)+length(cent5))~=5*length(cent1))
    cent1,cent2,cent3,cent4,cent5,
    A=[];
    fail=1;
    return;
end
W = [cent1' cent2' cent3' cent4' cent5'];
A=W;
t=toc;
end

% end