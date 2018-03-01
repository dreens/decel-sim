%% Plots results
function resultsdecel(rs)
    figure('Position',[0,0,2000,2000])
    colors = get(gca,'ColorOrder');
    colors = [colors ; (1-colors) ; colors ; colors];
    colors = jet(length(rs));
    for i=1:length(rs)
        r = rs(i);
        subplot(2,3,1); hold on
        n = sum(abs(rs(i).phases)==70);
        plot(r.molnum/r.molnum(1),'b-','DisplayName',['st=' num2str(n)],'Color',colors(i,:));
        subplot(2,3,2); hold on
        plot(0:length(r.modes),[r.initvz r.vels],'b-','Color',colors(i,:));
        subplot(2,3,3); hold on
        plot(diff([0 r.times])*1e6,'b-','Color',colors(i,:));
        subplot(2,3,4); hold on
        plot(r.pos(:,3)*1e3,r.vel(:,3),'b.','Color',colors(i,:),'MarkerSize',3);
        subplot(2,3,5); hold on
        plot(r.pos(:,1)*1e3,r.vel(:,1),'b.','Color',colors(i,:),'MarkerSize',3);
        subplot(2,3,6); hold on
        %plot(r.pos(:,2)*1e3,r.vel(:,2),'b.','Color',colors(i,:),'MarkerSize',3);
        %errorbar(n,r.numleft,sqrt(r.numleft),'bo','Color',colors(i,:),'MarkerSize',10);
    end
    subplot(2,3,1)
    xlabel('Stage Number','FontSize',12)
    ylabel('Population (arb)','FontSize',12)
    title('Molecules Remaining','FontSize',14)
    set(gca,'FontSize',12)
    set(gca,'YScale','log')
    grid on
    legend('show')
    
    subplot(2,3,2)
    xlabel('Stage Number','FontSize',12)
    ylabel('Velocity (m/s)','FontSize',12)
    titlestring = sprintf('Collision Hunt: Phi=%2.1f for 0-100 stages\n%s',...
        rs(1).phases(1),'Synchronous Velocity');
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

    
    
    
    subplot (2,3,5)
    diamL = 2.5e-3;
    distL = 8e-3 + r.pos(1,3);
    maxes = zeros(size(rs));
    for j=1:length(rs)
        r = rs(j);
        tof = zeros(1,1000);
        times = 1e-7:1e-7:1e-4;
        for i=1:1000
            t = times(i);
            zsq = (r.pos(:,3)-distL + t*r.vel(:,3)).^2;
            xsq = (r.pos(:,1)+t*r.vel(:,1)).^2;
            ya = abs(r.pos(:,2)+t*r.vel(:,2));
            tof(i) = sum( zsq + xsq < diamL^2/4 & ya < 2e-3);
        end
        plot(times*1e6,tof,'Color',colors(j,:))
        maxes(j) = max(tof);
    end
    xlabel('Time \mus')
    ylabel('Population (arb)')
    title('Time of Flight')
    grid on
    
    
    
    
    subplot(2,3,6)
    hold on
    errorbar(10:10:100,maxes,sqrt(maxes),'bx','Color',colors(end-1,:),'MarkerSize',10)
    errorbar(10:10:100,[rs.numleft],sqrt([rs.numleft]),'bo','Color',colors(2,:),'MarkerSize',10);
    xlabel('Initial Deceleration Stages','FontSize',12)
    ylabel('Remaining Population','FontSize',12)
    title('   Phase Space Acceptance','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    %legend('show')
end