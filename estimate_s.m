% Underdetermined Blind Source Separation Based on Subspace Representation
% Author: SangGyun Kim
%         MMP, EECS, KAIST
% Please, DO NOT share this code with anyone outside of the lab.


function U = estimate_s(X,W,n,m,ep)

% for snr = 50:-5:5          % Simulations under various SNRs
%     snr
%     SNR = num2str(snr);
%     ep = 1/30;            % Epsilon (SSD),defalut is 1/30
    L = size(X,2);             % The length of sources


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
            if norm([real(Xf{1}(f,t)) real(Xf{2}(f,t)) real(Xf{3}(f,t))]) < 0.001 | norm([imag(Xf{1}(f,t)) imag(Xf{2}(f,t)) imag(Xf{3}(f,t))]) < 0.001
                Xf{1}(f,t) = 0;
                Xf{2}(f,t) = 0;
                Xf{3}(f,t) = 0;
            end
        end
    end

    % Second, SSS
    SS1 = [];
    SS2 = [];
    SS3 = [];
    for t = 1:T
        for f = 2:F/2+1
            if Xf{1}(f,t) ~= 0
%                if (abs(imag(Xf{2}(f,t)/Xf{1}(f,t))) < ep) & (abs(imag(Xf{3}(f,t)/Xf{1}(f,t))) < ep) & (abs(imag(Xf{3}(f,t)/Xf{2}(f,t))) < ep)
                if (abs(angle(Xf{2}(f,t)/Xf{1}(f,t))) < ep || abs(angle(Xf{2}(f,t)/Xf{1}(f,t))) > pi-ep) & (abs(angle(Xf{3}(f,t)/Xf{1}(f,t))) < ep || abs(angle(Xf{3}(f,t)/Xf{1}(f,t))) > pi-ep) & (abs(angle(Xf{3}(f,t)/Xf{2}(f,t))) < ep || abs(angle(Xf{3}(f,t)/Xf{2}(f,t))) > pi-ep)
                    SS1 = [SS1 Xf{1}(f,t)];
                    SS2 = [SS2 Xf{2}(f,t)];
                    SS3 = [SS3 Xf{3}(f,t)];
                end
            end
        end
    end

    SS = [SS1; SS2; SS3];

    % ==================
    % Estimation of error
    % ==================
    SSa = [real(SS) imag(SS)];

    for nn = 1:length(SSa(1,:))
        SSa1(:,nn) = SSa(:,nn)/SSa(1,nn);
        SSa1(:,nn) = SSa1(:,nn)/norm(SSa1(:,nn),2);
    end

    
    % Kmeans clustering
    [IDX,C] = kmeans(SSa1',n,'emptyaction','singleton');          
%     A(:,4) = -A(:,4);
%     [IDX,C] = kmeans(SSa1',4,'start',A');         % Performance of Kmeans depends on initialization. 
                                                    % You can use another
                                                    % more efficient method
                                                    % for clustering.
                                                    % I recommend it.



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
        iter
        ITER = num2str(iter);

        % Iterations
        % Source Estimation

        for f = 1:1
            % Processing sample by sample
            for t = 1:T

                % MMSE based estimation
                xr = [real(Xf{1}(t,f)); real(Xf{2}(t,f)); real(Xf{3}(t,f))];
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
                        Pzr_cand(ii) = pz(zr_cand(ii),xr,W,[betar{1}(f); betar{2}(f); betar{3}(f); betar{4}(f)],[sigr{1}(f); sigr{2}(f); sigr{3}(f); sigr{4}(f)],[mur{1}(f); mur{2}(f); mur{3}(f); mur{4}(f)]);
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

                % Estimted sources
                Shat{1}(t,f) = shatr(1);
                Shat{2}(t,f) = shatr(2);
                Shat{3}(t,f) = shatr(3);
                Shat{4}(t,f) = shatr(4);
            end

        end

        for f = 2:F/2+1
            % Processing sample by sample
            for t = 1:T

                % Data augmentation (DA) algorithm
                xr = [real(Xf{1}(t,f)); real(Xf{2}(t,f)); real(Xf{3}(t,f))];
                xi = [imag(Xf{1}(t,f)); imag(Xf{2}(t,f)); imag(Xf{3}(t,f))];
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
                    for ii = 1:length(zr_cand)
                        Pzr_cand(ii) = pz(zr_cand(ii),xr,W,[betar{1}(f); betar{2}(f); betar{3}(f); betar{4}(f)],[sigr{1}(f); sigr{2}(f); sigr{3}(f); sigr{4}(f)],[mur{1}(f); mur{2}(f); mur{3}(f); mur{4}(f)]);
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
                    for ii = 1:length(zi_cand)
                        Pzi_cand(ii) = pz(zi_cand(ii),xi,W,[betai{1}(f); betai{2}(f); betai{3}(f); betai{4}(f)],[sigi{1}(f); sigi{2}(f); sigi{3}(f); sigi{4}(f)],[mui{1}(f); mui{2}(f); mui{3}(f); mui{4}(f)]);
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
                Shat{1}(t,f) = shatr(1) + sqrt(-1)*shati(1);
                Shat{2}(t,f) = shatr(2) + sqrt(-1)*shati(2);
                Shat{3}(t,f) = shatr(3) + sqrt(-1)*shati(3);
                Shat{4}(t,f) = shatr(4) + sqrt(-1)*shati(4);
            end

        end

        for w = 1:F/2-1
            Shat{1}(:,F/2+1+w) = conj(Shat{1}(:,F/2+1-w));
            Shat{2}(:,F/2+1+w) = conj(Shat{2}(:,F/2+1-w));
            Shat{3}(:,F/2+1+w) = conj(Shat{3}(:,F/2+1-w));
            Shat{4}(:,F/2+1+w) = conj(Shat{4}(:,F/2+1-w));
        end

        % Inverse DFT
        u1 = resynth(Shat{1});
        u2 = resynth(Shat{2});
        u3 = resynth(Shat{3});
        u4 = resynth(Shat{4});
        U = [u1; u2; u3; u4];
        % Save the results
%        save(['ubss_speech_sir_vs_snr' '_iter' ITER],'W','Shat','betar','betai','S','A','U','N','sigr','sigi','mur','mui','Zr','Zi');

    end

