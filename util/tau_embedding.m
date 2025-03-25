function [tau,tau_indx] = tau_embedding(signal,t,max_tau,dimension)
% This function computes the time delay for the Takens embedding based on
% space filling.

    if dimension < 2 || dimension > 3
        disp('Selected dimension not implemented');
    end
    
    % Iteration parameters
    n_iter = 750;
    n_bins = round(0.75*length(signal));
    
    % Delay vector and space fill vector
    tau_vec = 0:max_tau/n_iter:max_tau;
    tau_vec = tau_vec(2:end);
    fill_factor_vec = zeros(1,length(tau_vec));
    
    % Space fill factor parameters
    epsilon = 0.1;
    half_box = max(abs(signal)) + epsilon;
    if dimension == 2
        delta_bin = 2*half_box/sqrt(n_bins);
    else
        delta_bin = 2*half_box/(n_bins^(1/3));
    end

    % Bin limits
    bin_edges = -half_box:delta_bin:half_box;

    for k=1:length(tau_vec)
        % Find index and perform time-delay
        indx = find(t >= (tau_vec(k)+t(1)), 1);
        if dimension == 2
            tmp1 = signal(1:end-indx);
            tmp2 = signal(indx+1:end);
            % Count number of bins of trajectory
            [N,~,~] = histcounts2(tmp1,tmp2,bin_edges,bin_edges);
            bin_count = nnz(N);
            
            % Count number of bins of convex hull of trajectory
            conv_hull_count = nnz(bwconvhull(N > 0));
        else
            dE = 3;
            emb{1} = signal((dE-1)*indx+1:end);
            for i = 2:dE
                emb{i} = signal((dE-i)*indx+1:end-(i-1)*indx);
            end
            tmp1 = emb{1};
            tmp2 = emb{2};
            tmp3 = emb{3};
            
            % Count number of bins of trajectory
            N = histcnd(tmp1,tmp2,tmp3,bin_edges,bin_edges,bin_edges);
            bin_count = nnz(N);

            for i=1:size(N,1)
                N(i,:,:) = bwconvhull(squeeze(N(i,:,:) > 0));
            end

            for i=1:size(N,2)
                N(:,i,:) = bwconvhull(squeeze(N(:,i,:) > 0));
            end

            for i=1:size(N,3)
                N(:,:,i) = bwconvhull(squeeze(N(:,:,i) > 0));
            end
            
            % Count number of bins of convex hull of trajectory
            conv_hull_count = nnz(N);

        end

        % Compute fill factor
        fill_factor_vec(k) = bin_count/conv_hull_count;
    end

    [~,indx_tau_min] = min(fill_factor_vec);
    tau = tau_vec(indx_tau_min);

    tau_indx = sum(t < (tau+t(1)));

end


