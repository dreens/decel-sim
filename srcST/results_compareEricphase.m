%% Plots results
function results_compareEricphase(rs)
    l = length(rs);
    figure('position',[50,50,1000,1300])
    subplot(l/2,2,1)
    colors = get(gca,'ColorOrder');
    colors = [colors ; colors ; colors ; colors];
    ms = 1;
    finals = zeros(1,l);
    
    %rs = rs([5 6 3 4 1 2]);

    for i=1:length(rs)
        r = rs(i);
        subplot(l/2,2,i)
        m = trimmean(r.pos(:,3),95);
        keep = abs(r.pos(:,3) - m) < 5.5e-3 & abs(r.vel(:,3)-380)<45;
        finals(i) = sum(keep);
        plot((r.pos(keep,3)-m)*1e3,r.vel(keep,3),'b.',...
            'Color',colors(2-mod(i,2),:),...
            'MarkerSize',ms)
        ylim([340 420])
        xlim([-5.5 5.5])
        if ~mod(i,2)
            text(1,345,['Gain = ' num2str(finals(i)/finals(i-1))])
        end
        text(-5.4,412,sprintf('%d Stages\nS=%d',r.stages,max(r.modeseq)))
        if i==1
            title('Comparing S=1,3 Bunching, Gaussian Initial Phase Space')
        end
        if i==3
            ylabel('V_z (m/s)')
        end
        if i>4
            xlabel('z (mm)')
        end
    end
    
    
    
%     xlabel('Z Position Initial (mm)','FontSize',12)
%     ylabel('Z Velocity Initial (m/s)','FontSize',12)
%     title('Phase Space Z','FontSize',14)
%     set(gca,'FontSize',12)
%     grid on
end