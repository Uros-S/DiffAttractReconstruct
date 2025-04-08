function figID = space_filling_Lorenz(y_true, y_measured, y_D, y_G, coord_upper, coord_lower,figID)
% This function computes the space-filling properties of the Lorenz system.
        
    % Bin width parameters    
    delta_max = 0.35;
    delta_min = 0.01;
    delta_step = 0.025;
    delta = delta_min:delta_step:delta_max;

    epsilon = 0.1;

    bincount1 = zeros(1,length(delta));
    bincount2 = zeros(1,length(delta));
    bincount3 = zeros(1,length(delta));
    if ~isnan(y_G)
        bincount4 = zeros(1,length(delta));
    end

    % Data coordinates
    X1 = y_true(:,1);
    Y1 = y_true(:,2);
    X2 = y_measured(:,1);
    Y2 = y_measured(:,2);
    X3 = y_D(:,1);
    Y3 = y_D(:,2);
    if ~isnan(y_G)
        X4 = y_G(:,1);
        Y4 = y_G(:,2);
    end
    
    % Box containing all the data points
    if ~isnan(y_G)
        X_lim = max([max(abs(X1)),max(abs(X2)),max(abs(X3)),max(abs(X4))])+epsilon;
        Y_lim = max([max(abs(Y1)),max(abs(Y2)),max(abs(Y3)),max(abs(Y4))])+epsilon;
    else
        X_lim = max([max(abs(X1)),max(abs(X2)),max(abs(X3))])+epsilon;
        Y_lim = max([max(abs(Y1)),max(abs(Y2)),max(abs(Y3))])+epsilon;
    end

    for i=1:length(delta)

        % Bin limits
        Xedges = -X_lim:delta(i):X_lim;
        Yedges = -Y_lim:delta(i):Y_lim;

        [N,~,~] = histcounts2(X1,Y1,Xedges,Yedges);
        bincount1(i) = nnz(N);

        [N,~,~] = histcounts2(X2,Y2,Xedges,Yedges);
        bincount2(i) = nnz(N);

        [N,~,~] = histcounts2(X3,Y3,Xedges,Yedges);
        bincount3(i) = nnz(N);
        
        if ~isnan(y_G)
            [N,~,~] = histcounts2(X4,Y4,Xedges,Yedges);
            bincount4(i) = nnz(N);
        end
    end

    figID = figID + 1;
    figure(figID);
    plot(delta,bincount1,'LineWidth',3);
    hold on;
    plot(delta,bincount2,'Color','#77AC30','LineWidth',3);
    hold on;
    plot(delta,bincount3,'Color','#D95319','LineWidth',3);
    if ~isnan(y_G)
        hold on;
        plot(delta,bincount4,'LineWidth',3);
    end
    hold on;
    yline(size(y_true,1),'--','LineWidth',3);
    hold off;
    xlabel('Bin width'); ylabel('Number of bins');
    title('Global space-filling');
    if ~isnan(y_G)
        legend('Noise-free data','Measured data','Differentiator','Grassberger','Number of data points');
    else
        legend('Noise-free data','Measured data','Differentiator','Number of data points');
    end
    xlim([0,max(delta)]);
    ylim([0,1.1*size(y_true,1)]);
    ax = gca;
    ax.FontSize = 35;

    % Space-filling metrics for holes
    if ~isnan(coord_upper)
        bincount_up1 = zeros(1,length(delta));
        bincount_up2 = zeros(1,length(delta));
        bincount_up3 = zeros(1,length(delta));
        if ~isnan(y_G)
            bincount_up4 = zeros(1,length(delta));
        end

        X1_upper = min(coord_upper(1,:));
        X2_upper = max(coord_upper(1,:));
        Y1_upper = min(coord_upper(2,:));
        Y2_upper = max(coord_upper(2,:));

        for i=1:length(delta)
            % Bin limits
            Xedges = X1_upper:delta(i):X2_upper;
            Yedges = Y1_upper:delta(i):Y2_upper;

            [N,~,~] = histcounts2(X1,Y1,Xedges,Yedges);
            bincount_up1(i) = nnz(N);

            [N,~,~] = histcounts2(X2,Y2,Xedges,Yedges);
            bincount_up2(i) = nnz(N);

            [N,~,~] = histcounts2(X3,Y3,Xedges,Yedges);
            bincount_up3(i) = nnz(N);
            
            if ~isnan(y_G)
                [N,~,~] = histcounts2(X4,Y4,Xedges,Yedges);
                bincount_up4(i) = nnz(N);
            end
        end

        figID = figID + 1;
        figure(figID);
        plot(delta,bincount_up1,'LineWidth',3);
        hold on;
        plot(delta,bincount_up2,'Color','#77AC30','LineWidth',3);
        hold on;
        plot(delta,bincount_up3,'Color','#D95319','LineWidth',3);
        if ~isnan(y_G)
            hold on;
            plot(delta,bincount_up4,'LineWidth',3);
        end
        xlabel('Bin width'); ylabel('Number of bins');
        title('Upper hole space-filling');
        if ~isnan(y_G)
            legend('Noise-free data','Measured data','Differentiator','Grassberger');
        else
            legend('Noise-free data','Measured data','Differentiator');
        end
        xlim([0,max(delta)]);
        ax = gca;
        ax.FontSize = 35;
    end

    if ~isnan(coord_lower)
        bincount_low1 = zeros(1,length(delta));
        bincount_low2 = zeros(1,length(delta));
        bincount_low3 = zeros(1,length(delta));
        if ~isnan(y_G)
            bincount_low4 = zeros(1,length(delta));
        end

        X1_lower = min(coord_lower(1,:));
        X2_lower = max(coord_lower(1,:));
        Y1_lower = min(coord_lower(2,:));
        Y2_lower = max(coord_lower(2,:));

        for i=1:length(delta)
            % Bin limits
            Xedges = X1_lower:delta(i):X2_lower;
            Yedges = Y1_lower:delta(i):Y2_lower;

            [N,~,~] = histcounts2(X1,Y1,Xedges,Yedges);
            bincount_low1(i) = nnz(N);

            [N,~,~] = histcounts2(X2,Y2,Xedges,Yedges);
            bincount_low2(i) = nnz(N);

            [N,~,~] = histcounts2(X3,Y3,Xedges,Yedges);
            bincount_low3(i) = nnz(N);
            
            if ~isnan(y_G)
                [N,~,~] = histcounts2(X4,Y4,Xedges,Yedges);
                bincount_low4(i) = nnz(N);
            end
        end

        figID = figID + 1;
        figure(figID);
        plot(delta,bincount_low1,'LineWidth',3);
        hold on;
        plot(delta,bincount_low2,'Color','#77AC30','LineWidth',3);
        hold on;
        plot(delta,bincount_low3,'Color','#D95319','LineWidth',3);
        if ~isnan(y_G)
            hold on;
            plot(delta,bincount_low4,'LineWidth',3);
        end
        xlabel('Bin width'); ylabel('Number of bins');
        title('Lower hole space-filling');
        if ~isnan(y_G)
            legend('Noise-free data','Measured data','Differentiator','Grassberger');
        else
            legend('Noise-free data','Measured data','Differentiator');
        end
        xlim([0,max(delta)]);
        ax = gca;
        ax.FontSize = 35;
    end

end