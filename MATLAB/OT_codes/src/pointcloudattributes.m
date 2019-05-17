function pointcloudattributes(xtext,ytext,legendtext)
    xlabel(xtext);ylabel(ytext);
    grid on;
    set(gca,'fontsize',10,'fontweight','bold');
    if isempty(legendtext)~=1
        legend(legendtext,...
            'Orientation','horizontal','Position',[0.43 0.02 0.1778 0.04],'FontSize',14);
    end
    hold off
end
% 'Location', 'southwestoutside',