function figID = error_comparison_Lorenz(y_delayed,y_delayed_D,y_delayed_D_f,y_delayed_G,Nbins,figID)
% This function computes absolute and relative error for Lorenz system
% data.

    % Error computation differentiator
    abs_error_D = zeros(1,size(y_delayed_D,1));
    rel_error_D = zeros(1,size(y_delayed_D,1));
    for i=1:size(y_delayed_D,1)
        abs_error_D(i) = norm([y_delayed(i,1)-y_delayed_D(i,1);y_delayed(i,2)-y_delayed_D(i,2)]);
        rel_error_D(i) = abs_error_D(i)/norm([y_delayed(i,1);y_delayed(i,2)]);
    end

    % Error computation differentiator with SG filter
    abs_error_D_f = zeros(1,size(y_delayed_D_f,1));
    rel_error_D_f = zeros(1,size(y_delayed_D_f,1));
    for i=1:size(y_delayed_D_f,1)
        abs_error_D_f(i) = norm([y_delayed(i,1)-y_delayed_D_f(i,1);y_delayed(i,2)-y_delayed_D_f(i,2)]);
        rel_error_D_f(i) = abs_error_D_f(i)/norm([y_delayed(i,1);y_delayed(i,2)]);
    end
    
    if ~isnan(y_delayed_G)
        % Error computation Grassberger denoising
        abs_error_G = zeros(1,size(y_delayed_G,1));
        rel_error_G = zeros(1,size(y_delayed_G,1));
    
        for i=1:size(y_delayed_G,1)
            abs_error_G(i) = norm([y_delayed(i,1)-y_delayed_G(i,1);y_delayed(i,2)-y_delayed_G(i,2)]);
            rel_error_G(i) = abs_error_G(i)/norm([y_delayed(i,1);y_delayed(i,2)]);
        end
    end

    %% Plot histogram absolute error
    figID = figID+1;
    figure(figID);
    if ~isnan(y_delayed_G)
        subplot(1,2,1);
    end
    h1 = histogram(abs_error_D);
    h1.Normalization = 'probability';
    h1.NumBins = Nbins;
    h1.FaceColor = 'blue';
    hold on;
    h2 = histogram(abs_error_D_f);
    h2.Normalization = 'probability';
    h2.BinEdges = h1.BinEdges;
    h2.FaceColor = 'green';
    if ~isnan(y_delayed_G)
        dim = [0.3,0.35,0.3,0.3];
    else
        dim = [0.6,0.48,0.3,0.3];
    end
    str = {strcat(['Average error: ',num2str(round(mean(abs_error_D),2,'significant'))]), ...
           strcat(['Median error: ',num2str(round(median(abs_error_D),2,'significant'))]), ...
           strcat(['Worst-case error: ',num2str(round(max(abs_error_D),2,'significant'))])};
    annotation('textbox',dim,'String',str,'FitBoxToText','on', ...
               'BackgroundColor','blue','FaceAlpha',0.1,'FontSize',20);
    
    if ~isnan(y_delayed_G)
        dim = [0.3,0.21,0.3,0.3];
    else
        dim = [0.6,0.35,0.3,0.3];
    end
    str = {strcat(['Average error: ',num2str(round(mean(abs_error_D_f),2,'significant'))]), ...
           strcat(['Median error: ',num2str(round(median(abs_error_D_f),2,'significant'))]), ...
           strcat(['Worst-case error: ',num2str(round(max(abs_error_D_f),2,'significant'))])};
    annotation('textbox',dim,'String',str,'FitBoxToText','on', ...
               'BackgroundColor','green','FaceAlpha',0.1,'FontSize',20);

    title(['Absolute error',newline,'Differentiator']);
    legend('Differentiator','Differentiator with SG filter');
    ax = gca;
    ax.FontSize = 25;
    set(gca,'YScale','log');
    pbaspect([1.2,1,1])
    
    if ~isnan(y_delayed_G)
        subplot(1,2,2);
        h3 = histogram(abs_error_G);
        h3.Normalization = 'probability';   
        h3.NumBins = Nbins;
        h3.FaceColor = '#D95319';
    
        dim = [0.7,0.4,0.3,0.3];
        str = {strcat(['Average error: ',num2str(round(mean(abs_error_G),2,'significant'))]), ...
        strcat(['Median error: ',num2str(round(median(abs_error_G),2,'significant'))]), ...
        strcat(['Worst-case error: ',num2str(round(max(abs_error_G),2,'significant'))])};
        annotation('textbox',dim,'String',str,'FitBoxToText','on', ...
               'BackgroundColor','#D95319','FaceAlpha',0.1,'FontSize',20);
    
        title(['Absolute error',newline,'Schreiber-Grassberger']);
        ax = gca;
        ax.FontSize = 25;
        set(gca,'YScale','log');
        pbaspect([1.2,1,1])
    end

    %% Plot histogram relative error
    figID = figID+1;
    figure(figID);
    if ~isnan(y_delayed_G)
        subplot(1,2,1);
    end
    h1 = histogram(rel_error_D);
    h1.Normalization = 'probability';
    h1.NumBins = Nbins;
    h1.FaceColor = 'blue';
    hold on;
    h2 = histogram(rel_error_D_f);
    h2.Normalization = 'probability';
    h2.BinEdges = h1.BinEdges;
    h2.FaceColor = 'green';

    if ~isnan(y_delayed_G)
        dim = [0.3,0.35,0.3,0.3];
    else
        dim = [0.6,0.48,0.3,0.3];
    end
    str = {strcat(['Average error: ',num2str(round(mean(rel_error_D),2,'significant'))]), ...
           strcat(['Median error: ',num2str(round(median(rel_error_D),2,'significant'))]), ...
           strcat(['Worst-case error: ',num2str(round(max(rel_error_D),2,'significant'))])};
    annotation('textbox',dim,'String',str,'FitBoxToText','on', ...
               'BackgroundColor','blue','FaceAlpha',0.1,'FontSize',20);

    if ~isnan(y_delayed_G)
        dim = [0.3,0.21,0.3,0.3];
    else
        dim = [0.6,0.35,0.3,0.3];
    end
    str = {strcat(['Average error: ',num2str(round(mean(rel_error_D_f),2,'significant'))]), ...
           strcat(['Median error: ',num2str(round(median(rel_error_D_f),2,'significant'))]), ...
           strcat(['Worst-case error: ',num2str(round(max(rel_error_D_f),2,'significant'))])};
    annotation('textbox',dim,'String',str,'FitBoxToText','on', ...
               'BackgroundColor','green','FaceAlpha',0.1,'FontSize',20);

    title(['Relative error',newline,'Differentiator']);
    legend('Differentiator','Differentiator with SG filter');
    ax = gca;
    ax.FontSize = 25;
    set(gca,'YScale','log');
    pbaspect([1.2,1,1])
    
    if ~isnan(y_delayed_G)
        subplot(1,2,2);
        h3 = histogram(rel_error_G);
        h3.Normalization = 'probability';   
        h3.NumBins = Nbins;
        h3.FaceColor = '#D95319';
    
        dim = [0.7,0.4,0.3,0.3];
        str = {strcat(['Average error: ',num2str(round(mean(rel_error_G),2,'significant'))]), ...
        strcat(['Median error: ',num2str(round(median(rel_error_G),2,'significant'))]), ...
        strcat(['Worst-case error: ',num2str(round(max(rel_error_G),2,'significant'))])};
        annotation('textbox',dim,'String',str,'FitBoxToText','on', ...
               'BackgroundColor','#D95319','FaceAlpha',0.1,'FontSize',20);
    
        title(['Relative error',newline,'Schreiber-Grassberger']);
        ax = gca;
        ax.FontSize = 25;
        set(gca,'YScale','log');
        pbaspect([1.2,1,1])
    end

end