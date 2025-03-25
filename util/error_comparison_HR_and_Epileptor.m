function figID = error_comparison_HR_and_Epileptor(x,x_D,x_fD,Nbins,figID)
% This function computes absolute and relative error for the HR system
% data.

    % Error computation differentiator
    abs_error_D = zeros(1,size(x_D,1));
    rel_error_D = zeros(1,size(x_D,1));
    abs_error_fD = zeros(1,size(x_fD,1));
    rel_error_fD = zeros(1,size(x_fD,1));
    for i=1:size(x_D,1)
        abs_error_D(i) = norm([x(i,1)-x_D(i,1);x(i,2)-x_D(i,2);x(i,3)-x_D(i,3)]);
        rel_error_D(i) = abs_error_D(i)/norm(x(i,:)');
        abs_error_fD(i) = norm([x(i,1)-x_fD(i,1);x(i,2)-x_fD(i,2);x(i,3)-x_fD(i,3)]);
        rel_error_fD(i) = abs_error_fD(i)/norm(x(i,:)');
    end

    %% Plot histogram absolute error
    figID = figID+1;
    figure(figID);
    h1 = histogram(abs_error_D);
    h1.Normalization = 'probability';
    h1.NumBins = Nbins;
    h1.FaceColor = 'blue';

    hold on;

    h2 = histogram(abs_error_fD);
    h2.Normalization = 'probability';
    h2.BinEdges = h1.BinEdges;
    h2.FaceColor = 'green';

    dim = [0.5,0.45,0.3,0.3];
    str = {strcat(['Average error: ',num2str(round(mean(abs_error_D),2,'significant'))]), ...
    strcat(['Median error: ',num2str(round(median(abs_error_D),2,'significant'))]), ...
    strcat(['Worst-case error: ',num2str(round(max(abs_error_D),2,'significant'))])};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','blue','FaceAlpha',0.1,'FontSize',30);

    dim = [0.5,0.44,0.3,0.1];
    str = {strcat(['Average error: ',num2str(round(mean(abs_error_fD),2,'significant'))]), ...
    strcat(['Median error: ',num2str(round(median(abs_error_fD),2,'significant'))]), ...
    strcat(['Worst-case error: ',num2str(round(max(abs_error_fD),2,'significant'))])};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','green','FaceAlpha',0.1,'FontSize',30);

    title('Absolute error');
    legend('Differentiator','Differentiator with SG filter');
    ax = gca;
    ax.FontSize = 35;
    set(gca,'YScale','log');
    pbaspect([1.2,1,1])

    %% Plot histogram relative error
    figID = figID+1;
    figure(figID);
    h1 = histogram(rel_error_D);
    h1.Normalization = 'probability';
    h1.NumBins = Nbins;
    h1.FaceColor = 'blue';

    hold on;

    h2 = histogram(rel_error_fD);
    h2.Normalization = 'probability';
    h2.BinEdges = h1.BinEdges;
    h2.FaceColor = 'green';

    dim = [0.5,0.45,0.3,0.3];
    str = {strcat(['Average error: ',num2str(round(mean(rel_error_D),2,'significant'))]), ...
    strcat(['Median error: ',num2str(round(median(rel_error_D),2,'significant'))]), ...
    strcat(['Worst-case error: ',num2str(round(max(rel_error_D),2,'significant'))])};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','blue','FaceAlpha',0.1,'FontSize',30);
    
    dim = [0.5,0.44,0.3,0.1];
    str = {strcat(['Average error: ',num2str(round(mean(rel_error_fD),2,'significant'))]), ...
    strcat(['Median error: ',num2str(round(median(rel_error_fD),2,'significant'))]), ...
    strcat(['Worst-case error: ',num2str(round(max(rel_error_fD),2,'significant'))])};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','green','FaceAlpha',0.1,'FontSize',30);

    title('Relative error');
    legend('Differentiator','Differentiator with SG filter');
    ax = gca;
    ax.FontSize = 35;
    set(gca,'YScale','log');
    pbaspect([1.2,1,1])

end