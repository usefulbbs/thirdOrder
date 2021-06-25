% Underdetermined Blind Source Separation Based on Subspace Representation
% Author: SangGyun Kim
%         MMP, EECS, KAIST
%Modified by Liang Zou


function U = estimate_s_3(X,W)
U =[];

% for snr = 50:-5:5          % Simulations under various SNRs
%     snr
%     SNR = num2str(snr);
%     ep = 1/30;            % Epsilon (SSD),defalut is 1/30
L = size(X,2);             % The length of sources
[m,n] = size(W);

% DFT points
N =1024;
% DFT of mixture
for k = 1:m
    Xf{k} = analyze(X(k,:),N);
end
% Size of Xf
[F,T] = size(Xf{1});
% First, Autosource point selection
for t = 1:T
    for f = 2:F/2+1
        temp_real = [];
        temp_img =[];
        for k = 1:m
            temp_real = [temp_real real(Xf{k}(f,t))];
            temp_img = [temp_img imag(Xf{k}(f,t))];
        end
        if (norm(temp_real) < 0.001 || norm(temp_img) < 0.001)
            for k = 1:m
                Xf{k}(f,t) = 0;
            end
        end
    end
end

%pause;
% Lengthen the mixtures
Leng = ceil(L/N)*N;
Dleng = Leng - L;
X = [X zeros(m,Dleng)];

[T,F] = size(Xf{1});
% Parameter Updates
% Real part
for ii = 1:n
    for f = 1:F/2+1
        betar{ii}(f) = 1;
        sigr{ii}(f) = 1;
        mur{ii}(f)  = 0;
    end
end
% Imaginary part
for ii = 1:n
    for f = 1:F/2+1
        betai{ii}(f) = 1;
        sigi{ii}(f) = 1;
        mui{ii}(f)  = 0;
    end
end

for k = 1:m
    Xf{k} = analyze(X(k,:),N);
end
%=================================================
%                Main Algorithm
%=================================================

% SVD of W
[u,d,v] = svd(W);
inverA = pinv(W);
n_C = nchoosek(n,n-m);
N_C = nchoosek(1:n,n-m);

for iter = 1:1 % default is 3
%     iter
%     ITER = num2str(iter);
    
    % Iterations
    % Source Estimation
    
    for f = 1:1
        % Processing sample by sample
        for t = 1:T
            
            % MMSE based estimation
%             xr = [real(Xf{1}(t,f)); real(Xf{2}(t,f)); real(Xf{3}(t,f))];
            xr=[];
            for k = 1:m
                xr = [xr;real(Xf{k}(t,f))];
            end
            O = zeros(n,1);
            O = inverA*xr;
            
            % Minimum mean square error (MMSE) estimation
            mag_xr = sqrt(sum(xr.^2));
            if mag_xr == 0
                zr = [0];
                shatr = v(:,1:m)*inv(d(:,1:m))*u'*xr + v(:,m+1:n)*zr;           % Subspace representation
                Zr(t,f) = zr;
            else
                zr_cand = zeros(n-m,n_C);
                for ii=1:n_C
                    zr_cand(:,ii) = -(inv(v(N_C(ii,:),m+1:n))*O(N_C(ii,:)));
                end                         % This procedure reduces the computational burden of sampling
                for ii = 1:length(zr_cand)
                    all_betar =zeros(n,1);
                    all_sigr = zeros(n,1);
                    all_mur = zeros(n,1);
                    for k = 1:n
                        all_betar(k) = betar{k}(f);
                        all_sigr(k) = sigr{k}(f);
                        all_mur(k) = mur{k}(f);
                    end
                    Pzr_cand(ii) = pz(zr_cand(ii),xr,W,all_betar,all_sigr,all_mur);
                end
                if sum(Pzr_cand) == 0
                    f;
                    t;
                    shatr = bsr_linprog(xr,W);
                else
                    Pzr_cand = Pzr_cand/sum(Pzr_cand);
                    zr = sum(zr_cand.*Pzr_cand);
                    Zr(t,f) = zr;
                    
                    % Estimated sources
                    shatr = v(:,1:m)*inv(d(:,1:m))*u'*xr + v(:,m+1:n)*zr;
                    clc;
                end
            end
            
            % Estimted sources
            for k =1:n
                Shat{k}(t,f) = shatr(k);
            end
        end
        
    end
    
    for f = 2:F/2+1
        % Processing sample by sample
        for t = 1:T
            
            % Data augmentation (DA) algorithm
            xr = zeros(m,1);
            xi = zeros(m,1);
            for k =1:m
                xr(k) = real(Xf{k}(t,f));
                xi(k) = imag(Xf{k}(t,f));
            end
            Or = zeros(n,1);
            Or = inverA*xr;
            Oi = zeros(n,1);
            Oi = inverA*xi;
            
            % Real part
            mag_xr = sqrt(sum(xr.^2));
            if mag_xr == 0
                zr = [0];
                shatr = v(:,1:m)*inv(d(:,1:m))*u'*xr + v(:,m+1:n)*zr;
                Zr(t,f) = zr;
            else
                zr_cand = zeros(n-m,n_C);
                for ii=1:n_C
                    zr_cand(:,ii) = -(inv(v(N_C(ii,:),m+1:n))*Or(N_C(ii,:)));
                end                         % This procedure reduces the computational burden of sampling
                all_betar =zeros(n,1);
                all_sigr = zeros(n,1);
                all_mur = zeros(n,1);
                for k = 1:n
                    all_betar(k) = betar{k}(f);
                    all_sigr(k) = sigr{k}(f);
                    all_mur(k) = mur{k}(f);
                end
                for ii = 1:length(zr_cand)
                    Pzr_cand(ii) = pz(zr_cand(ii),xr,W,all_betar,all_sigr,all_mur);
                end
                if sum(Pzr_cand) == 0
                    f
                    t
                    shatr = bsr_linprog(xr,W);
                else
                    Pzr_cand = Pzr_cand/sum(Pzr_cand);
                    zr = sum(zr_cand.*Pzr_cand);
                    Zr(t,f) = zr;
                    
                    % Estimated sources
                    shatr = v(:,1:m)*inv(d(:,1:m))*u'*xr + v(:,m+1:n)*zr;
                end
            end
            % Imaginary part
            mag_xi = sqrt(sum(xi.^2));
            if mag_xi == 0
                zi = [0];
                shati = v(:,1:m)*inv(d(:,1:m))*u'*xi + v(:,m+1:n)*zi;
                Zi(t,f) = zi;
            else
                zi_cand = zeros(n-m,n_C);
                for ii=1:n_C
                    zi_cand(:,ii) = -(inv(v(N_C(ii,:),m+1:n))*Oi(N_C(ii,:)));
                end                         % This procedure reduces the computational burden of sampling
                all_betai =zeros(n,1);
                all_sigi = zeros(n,1);
                all_mui = zeros(n,1);
                for k = 1:n
                    all_betai(k) = betai{k}(f);
                    all_sigi(k) = sigi{k}(f);
                    all_mui(k) = mui{k}(f);
                end
                for ii = 1:length(zi_cand)
                    Pzi_cand(ii) = pz(zi_cand(ii),xi,W,all_betai,all_sigi,all_mui);
                end
                if sum(Pzi_cand) == 0
                    f
                    t
                    shati = bsr_linprog(xi,W);
                else
                    Pzi_cand = Pzi_cand/sum(Pzi_cand);
                    zi = sum(zi_cand.*Pzi_cand);
                    Zi(t,f) = zi;
                    
                    % Estimated sources
                    shati = v(:,1:m)*inv(d(:,1:m))*u'*xi + v(:,m+1:n)*zi;
                end
            end
            
            % Estimted sources
            for k =1:n
                Shat{k}(t,f) = shatr(k)+sqrt(-1)*shati(k);
            end
        end
        
    end
    
    for w = 1:F/2-1
        for k = 1:n
            Shat{k}(:,F/2+1+w) = conj(Shat{k}(:,F/2+1-w));
        end
    end
    
    % Inverse DFT
    u1 = resynth(Shat{1});
    U = zeros(k,length(u1));
    for k =1:n
        U(k,:) = resynth(Shat{k});
    end
    % Save the results
    %        save(['ubss_speech_sir_vs_snr' '_iter' ITER],'W','Shat','betar','betai','S','A','U','N','sigr','sigi','mur','mui','Zr','Zi');
    
end