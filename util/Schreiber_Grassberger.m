function x = Schreiber_Grassberger(signal,iterations,dE,tau,noise_lvl)
% This function implements the Schreiber-Grassberger de-noising technique
% for the chaotic Lorenz system. 

    D = 2.06;
    Q = 5;
    
    %% Time-delay embedding
    emb{1} = signal((dE-1)*tau+1:end,1);
    T = length(emb{1});
    for i = 2:dE
        emb{i} = signal((dE-i)*tau+1:end-(i-1)*tau,1);
    end
    x = zeros(T,dE);
    for i = 1:dE
        x(:,i) = emb{i};
    end
    
    % Setup parameters for de-noising
    K = ceil((noise_lvl*T^(2/D))^(2*D/(4+D)));
    
    R = ones(1,dE);
    R(1) = 1/sqrt(0.001);
    R(end) = 1/sqrt(0.001);
    
    disp(newline);
    disp('Starting Schreiber-Grassberger de-noising');

    for Z=1:iterations
        for I=1:T
            %% K-nearest neighbours mass centers and weighted covariance matrix
            % K-nearest neighbours center of mass
            Idx = knnsearch(x,x(I,:),'K',K);
            xi = mean(x(Idx,:));
        
            % Weighted covariance matrix
            Gamma = zeros(dE,dE);
            for i=1:dE
                for j=1:dE
                    for z=1:K
                        Gamma(i,j) = Gamma(i,j) + R(i)*1/K*(x(Idx(z),i)*x(Idx(z),j) - xi(i)*xi(j))*R(j);
                    end
                end
            end
            
            %% Correct the embedding vectors
            % Compute eigenvector corresponding to smallest eigenvalue
            [e_n,~] = eigs(Gamma,Q,'SM');
        
            % Compute correction
            theta = zeros(1,dE);
            for i=1:dE
                for j=1:dE
    
                    D_ij = 0;
                    for q=1:Q
                        D_ij = D_ij + e_n(i,q)*e_n(j,q);
                    end
    
                    theta(i) = theta(i) + 1/R(i)*D_ij*R(j)*(xi(j)-x(I,j));
                end
            end
            
            % Corrected point in embedded space
            x(I,:) = x(I,:) + theta;
    
            % To see progress
            if rem(I,100) == 0
                disp(newline);
                disp(['Progress: ',num2str(I),'/',num2str(T)]);
                disp(['Iteration: ',num2str(Z),'/',num2str(iterations)])
            end
        end
        
        % Decrease the number of neighbours for the next iteration
        K = round((1-0.15)*K);
    end

end