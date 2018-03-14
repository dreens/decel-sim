%% Plots results
function resultsdecel(rs)
    figure('Position',[0,0,2000,2000])
    colors = get(gca,'ColorOrder');
    colors = [colors ; (1-colors)];
    lines = repmat({'--','-'},1,10);
    markers = repmat({'o','x'},1,10);
    sname = repmat({'DSwitch','Normal'},1,6);
    for i=1:length(rs)
        r = rs(i);
        i2 = round(i/2+.25);
        c = colors(i2,:);
        l = lines{i};
        m = markers{i};
        subplot(2,3,1); hold on
        h = plot(r.molnum/r.molnum(1),'Color',c,'LineStyle',l);
        h.DisplayName = sprintf('v=%2.1f, %s',r.vels(end),sname{i});
        plot(r.numstage,r.molnum(end)/r.molnum(1),'Color',c,'Marker',m)
        
        subplot(2,3,2); hold on
        plot(0:max(r.stages),[r.initvz r.vels],'b-','Color',c,'LineStyle',l);
        
        subplot(2,3,3); hold on
        plot(r.pos(:,2)*1e3,r.vel(:,2),'b.','Color',c,'MarkerSize',3);
        
        subplot(2,3,4); hold on
        xx = 3*mod(i-1,2)-1.5;
        yy = 25*(i2-2);
        plot((r.pos(:,3)-r.pos(1,3))*1e3+xx,r.vel(:,3)+yy-50,'b.','Color',c,'MarkerSize',3);
        
        subplot(2,3,5); hold on
        diamL = 2.5e-3; 
        distL = 8e-3 + r.pos(1,3);
        tof = zeros(1,4000);
        times = (1e-7)*(1:length(tof));
        for j=1:length(tof)
            t = times(j);
            zsq = (r.pos(:,3)-distL + t*r.vel(:,3)).^2;
            xsq = (r.pos(:,1)+t*r.vel(:,1)).^2;
            ya = abs(r.pos(:,2)+t*r.vel(:,2));
            tof(j) = sum( zsq + xsq < diamL^2/4 & ya < 2e-3);
        end
        rs(i).tofpeak = max(tof);
        rs(i).tofarea = trapz(times,tof);
        plot(times*1e6,tof,'Color',c,'LineStyle',l)
    end
    
    subplot(2,3,6); hold on
    nls = [rs.numleft];
    tps = [rs.tofpeak];
    tpa = [rs.tofarea];
    xax = [rs(1:2:end).numstage];
    plot(xax,nls(2:2:end)./nls(1:2:end),'DisplayName','Total Ratio','Marker','x','LineStyle','-');
    plot(xax,tps(2:2:end)./tps(1:2:end),'DisplayName','Peak Ratio','Marker','o','LineStyle','-');
    plot(xax,tpa(2:2:end)./tpa(1:2:end),'DisplayName','Area Ratio','Marker','sq','LineStyle','-');
    legend('show')
    legend('Location','NorthWest')
    
    subplot(2,3,1)
    xlabel('Stage Number','FontSize',12)
    ylabel('Population (arb)','FontSize',12)
    title('Molecules Remaining','FontSize',14)
    set(gca,'FontSize',12)
    set(gca,'YScale','log')
    grid on
    
    subplot(2,3,2)
    xlabel('Stage Number','FontSize',12)
    ylabel('Velocity (m/s)','FontSize',12)
    titlestring = sprintf('%s\n%s, Phi=%2.1f',...
        'Delay Switching v Normal v Length',...
        'Slowing to 50 m/s', rs(1).phase);
    title(titlestring,'FontSize',14)
    set(gca,'FontSize',12)
    grid on

    subplot(2,3,3)
    xlabel('Y Position (mm)','FontSize',12)
    ylabel('Y Velocity (m/s)','FontSize',12)
    title('Phase Space Y','FontSize',14)
    set(gca,'FontSize',12)
    grid on

    subplot(2,3,4)
    xlabel('Z Position (mm)','FontSize',12)
    ylabel('Z Velocity (m/s)','FontSize',12)
    ylim([-40 40])
    xlim([-4 4])
    title('Phase Space Z','FontSize',14)
    set(gca,'FontSize',12)
    grid on

    subplot(2,3,5)
    xlabel('Time after turn-off (\mus)','FontSize',12)
    ylabel('Population in Laser','FontSize',12)
    title('Time of Flight','FontSize',14)
    set(gca,'FontSize',12)
    %set(gca,'XScale','log')
    xlim([50,350])
    grid on
        
    subplot(2,3,6)
    xlabel('Stage Number','FontSize',12)
    ylabel('Molecules','FontSize',12)
    title('Phase Space Totals','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    %legend('show')
end