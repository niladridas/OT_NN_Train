% Plot prior samples and mean
function pointcloudplot(X,xtruept)
    % Input data:
    % X: no_samples x state_size
    pointSize = 20;
    scatter(X(:,1).*cos(X(:,3)),X(:,1).*sin(X(:,3)),pointSize,'filled'); hold on;grid on;
%     b = gca;legend(b,'off');
    meanxy = mean(X);
    plot(meanxy(1,1)*cos(meanxy(1,3)),meanxy(1,1)*sin(meanxy(1,3)),'.c', 'MarkerSize',30);hold on;
    if isempty(xtruept)~=1
        plot(xtruept(1,1)*cos(xtruept(1,3)),xtruept(1,1)*sin(xtruept(1,3)),'.m', 'MarkerSize',30);hold on;
    end
    box on;
end