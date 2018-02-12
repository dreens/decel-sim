%% Plots results
function results_initphase(rs)
    figure('position',[50,50,1300,600])
    subplot(1,2,1)
    colors = get(gca,'ColorOrder');
    colors = [colors ; colors ; colors ; colors];
    
    hold on
    for i=1:length(rs)
        r = rs(i);
        if i==1
            plot(r.posiall(:,3)*1e3,r.veliall(:,3),'b.',...
                'Color',colors(1,:),'MarkerSize',2);
        end
        keep = r.vel(:,3) < 60;
        
        plot(r.posi(keep,3)*1e3,r.veli(keep,3),'b.','Color',...
            colors(i+1,:),'MarkerSize',10);
    end
    
    xlabel('Z Position Initial (mm)','FontSize',12)
    ylabel('Z Velocity Initial (m/s)','FontSize',12)
    title('Phase Space Z','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    
    subplot(1,2,2)
    hold on
    for i=1:length(rs)
        r = rs(i);
        if i==1
            plot(r.posiall(:,1)*1e3,r.veliall(:,1),'b.',...
                'Color',colors(1,:),'MarkerSize',2);
        end
        keep = r.vel(:,3) < 60;
        
        plot(r.posi(keep,1)*1e3,r.veli(keep,1),'b.','Color',...
            colors(i+1,:),'MarkerSize',10);
    end
    
    xlabel('X Position Initial (mm)','FontSize',12)
    ylabel('X Velocity Initial (m/s)','FontSize',12)
    title('Phase Space X','FontSize',14)
    set(gca,'FontSize',12)
    grid on

end