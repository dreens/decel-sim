%% Plots results
function resultsdeceltrap(rs)
    figure('position',[100,100,1100,800])
    colors = get(gca,'ColorOrder');
    colors = [colors ; colors ; colors ; colors];
    for i=1:length(rs)
        r = rs(i);
        j = i;
        subplot(2,3,1); hold on
        plot(r.molnum,'b-','Color',colors(j,:));
        subplot(2,3,2); hold on
        plot(r.vels,'b-','Color',colors(j,:));
        subplot(2,3,3); hold on
        plot(diff([0 r.times])*1e6,'b-','Color',colors(j,:));
        subplot(2,3,4); hold on
        plot(r.pos(:,3)*1e3,r.vel(:,3),'b.','Color',colors(j,:),'MarkerSize',1);
        subplot(2,3,5); hold on
        plot(r.pos(:,1)*1e3,r.vel(:,1),'b.','Color',colors(j,:),'MarkerSize',1);
        subplot(2,3,6); hold on
        plot(r.pos(:,2)*1e3,r.vel(:,2),'b.',...
            'DisplayName',sprintf('vf=%dm/s, load=%d\mus',r.finalvz,r.loadtime),...
            'Color',colors(j,:),'MarkerSize',1);
    end
    subplot(2,3,1)
    xlabel('Stage Number','FontSize',12)
    ylabel('Number','FontSize',12)
    title('Molecules Remaining','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    
    subplot(2,3,2)
    xlabel('Stage Number','FontSize',12)
    ylabel('Velocity (m/s)','FontSize',12)
    titlestring = sprintf('Phase Angle=%2.1f, Final V=%2.1f\n%s',...
        rs(1).phase,rs(1).finalvz,'Synchronous Velocity');
    title(titlestring,'FontSize',14)
    set(gca,'FontSize',12)
    grid on

    subplot(2,3,3)
    xlabel('Stage Number','FontSize',12)
    ylabel('Time (\mus)','FontSize',12)
    title('Stage Time','FontSize',14)
    set(gca,'FontSize',12)
    grid on

    subplot(2,3,4)
    xlabel('Z Position (mm)','FontSize',12)
    ylabel('Z Velocity (m/s)','FontSize',12)
    title('Phase Space Z','FontSize',14)
    set(gca,'FontSize',12)
    grid on

    subplot(2,3,5)
    xlabel('X Position (mm)','FontSize',12)
    ylabel('X Velocity (m/s)','FontSize',12)
    title('Phase Space X','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    
    subplot(2,3,6)
    xlabel('Y Position (mm)','FontSize',12)
    ylabel('Y Velocity (m/s)','FontSize',12)
    title('Phase Space Y','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    legend('show')
    
%     %tim's figure
%     figure;
%     hold on
%     for i=1:length(rs)
%         r = rs(i);
%         semilogy(r.molnum/r.molnum(1),'b-','Color',colors(i,:)); hold on;
%         ylim([0.01,1])
%         xlim([0,r.stages])
%         xlabel('Stage Number','FontSize',12)
%         ylabel('Number','FontSize',12)
%         title('Molecules Remaining','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%     end
    
end