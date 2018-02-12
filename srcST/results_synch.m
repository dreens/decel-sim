%% Plots results
function results_synch(rs)
    figure('position',[50,50,1000,600])
    colors = get(gca,'ColorOrder');
    colors = [colors ; colors ; colors ; colors];
    for j=1:length(rs)
        r = rs(j);
        subplot(2,3,1); hold on
        plot(r.times*1e6,r.molnum,'b-','Color',colors(j,:));
        subplot(2,3,2); hold on
        plot(r.times*1e6,r.vzsynch,'b-','Color',colors(j,:));
        subplot(2,3,3); hold on
        plot(r.times*1e6,r.zsynch,'b-','Color',colors(j,:));
        subplot(2,3,4); hold on
        plot(r.pos(:,3)*1e3,r.vel(:,3),'b.','Color',colors(j,:),'MarkerSize',5);
        subplot(2,3,5); hold on
        plot(r.times*1e6,r.ksynch(:,1),'b-','Color',colors(2*j-1,:));
        plot(r.times*1e6,r.ksynch(:,2),'b-','Color',colors(2*j,:));
        subplot(2,3,6); hold on
        plot(r.zsynch,r.pesynch,'b-','Color',colors(j,:));
    end
    subplot(2,3,1)
    xlabel('Time (\mus)','FontSize',12)
    ylabel('Number','FontSize',12)
    title('Molecules Remaining','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    hold on
    ys = ylim;
    plot([r.decelofftime r.decelofftime]*1e6,ys,'r-')
    plot([r.loadingofftime r.loadingofftime]*1e6,ys,'r-')
    
    subplot(2,3,2)
    xlabel('Time (\mus)','FontSize',12)
    ylabel('Velocity (m/s)','FontSize',12)
    titlestring = sprintf('%s=%2.1f, %s=%2.1f, %s=%s, %s=%s\n%s',...
        'Phase Angle',rs(1).phase,...
        'Final V',rs(1).finalvz,...
        'Decel',r.decel,...
        'Trap',r.trapname,...
        'Synchronous Velocity');
    title(titlestring,'FontSize',14)
    set(gca,'FontSize',12)
    grid on
    hold on
    ys = ylim;
    plot([r.decelofftime r.decelofftime]*1e6,ys,'r-')
    plot([r.loadingofftime r.loadingofftime]*1e6,ys,'r-')

    subplot(2,3,3)
    xlabel('Time (\mus)','FontSize',12)
    ylabel('Position','FontSize',12)
    title('Synchronous Position','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    hold on
    ys = ylim;
    plot([r.decelofftime r.decelofftime]*1e6,ys,'r-')
    plot([r.loadingofftime r.loadingofftime]*1e6,ys,'r-')

    subplot(2,3,4)
    xlabel('Z Position (mm)','FontSize',12)
    ylabel('Z Velocity (m/s)','FontSize',12)
    title('Phase Space Z','FontSize',14)
    set(gca,'FontSize',12)
    grid on

    subplot(2,3,5)
    xlabel('Time (\mus)','FontSize',12)
    ylabel('Focusing (GHz/mm^2)','FontSize',12)
    title('Focusing','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    hold on
    ys = ylim;
    plot([r.decelofftime r.decelofftime]*1e6,ys,'r-')
    plot([r.loadingofftime r.loadingofftime]*1e6,ys,'r-')
    
    subplot(2,3,6)
    xlabel('Distance (m)','FontSize',12)
    ylabel('Potential Energy (GHz)','FontSize',12)
    title('Potential Energy','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    hold on

end