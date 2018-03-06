%% Plots results
function resultsdecel(rs)
    figure('Position',[0,0,2000,2000])
    colors = get(gca,'ColorOrder');
    colors = [colors ; (1-colors) ; colors ; colors];
    %colors = jet(length(rs));
    for i=1:length(rs)
        r = rs(i);
        subplot(2,3,1); hold on
        plot(r.molnum/r.molnum(1),'b-','Color',colors(i,:));
        
        subplot(2,3,2); hold on
        plot(0:max(r.stages),[r.initvz r.vels],'b-','Color',colors(i,:));
        
        subplot(2,3,3); hold on
        plot(diff([0 r.times])*1e6,'b-','Color',colors(i,:));
        
        subplot(2,3,4); hold on
        plot(r.pos(:,3)*1e3,r.vel(:,3),'b.','Color',colors(i,:),'MarkerSize',3);
        
        subplot(2,3,5); hold on
        diamL = 2.5e-3; 
        distL = 8e-3 + r.pos(1,3);
        tof = zeros(1,1000);
        times = 1e-7:1e-7:1e-4;
        for j=1:1000
            t = times(j);
            zsq = (r.pos(:,3)-distL + t*r.vel(:,3)).^2;
            xsq = (r.pos(:,1)+t*r.vel(:,1)).^2;
            ya = abs(r.pos(:,2)+t*r.vel(:,2));
            tof(j) = sum( zsq + xsq < diamL^2/4 & ya < 2e-3);
        end
        plot(times*1e6,tof,'Color',colors(i,:))
        
        subplot(2,3,6); hold on
        plot(r.pos(:,2)*1e3,r.vel(:,2),'b.','Color',colors(i,:),'MarkerSize',3);
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
    titlestring = sprintf('Delay-Mode v Normal Switching\nPhi=%2d, v=%3.1f',...
        mode(abs(r.endphases)),r.vels(end));
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
    xlabel('Time after turn-off (\mus)','FontSize',12)
    ylabel('Population in Laser','FontSize',12)
    title('Time of Flight','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    
    
    
    subplot(2,3,6)
    xlabel('Y Position (mm)','FontSize',12)
    ylabel('Y Velocity (m/s)','FontSize',12)
    title('Phase Space Y','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    %legend('show')
end