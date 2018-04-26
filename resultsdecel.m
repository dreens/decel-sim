%% Plots results
function resultsdecel(rs)
    figure('Position',[0,0,1000,1000])
    colors = get(gca,'ColorOrder');
    colors = [colors([3 2 1 4 5 6 7],:) ; (1-colors)];
    lines = repmat({'-','-'},1,10);
    sname = repmat({'DSwitch','Normal'},1,6);
    for i=1:length(rs)
        r = rs(i);
        i2 = i;%round(i/2+.25);
        c = colors(i2,:);
        l = lines{i};
        subplot(2,3,1); hold on
        h = plot(r.molnum/r.molnum(1),'Color',c,'LineStyle',l);
        h.DisplayName = sprintf('v=%2.1f, %s',r.vels(end),sname{i});
        
        subplot(2,3,2); hold on
        plot(0:max(r.stages),[r.initvz r.vels],'b-','Color',c,'LineStyle',l);
        
        subplot(2,3,3); hold on
        plot(r.pos(:,1)*1e3,r.vel(:,1),'b.','Color',c,'MarkerSize',3);
        
        subplot(2,3,4); hold on
        xx = 0*10*mod(i-1,3);
        plot(r.pos(:,3)*1e3+xx,r.vel(:,3),'b.','Color',c,'MarkerSize',3);
        
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
    nls(4) = nls(4)/2.5;
    tps = [rs.tofpeak];
    tpa = [rs.tofarea];
    plot(nls(4:-1:1),'DisplayName','Total','Marker','x','LineStyle','-');
    plot(tps(4:-1:1),'DisplayName','Peak','Marker','o','LineStyle','-');
    plot(tpa(4:-1:1)*1e4,'DisplayName','Area','Marker','sq','LineStyle','-');
    legend('show')
    legend('Location','South')
    
    subplot(2,3,1)
    xlabel('Stage Number','FontSize',12)
    ylabel('Population (arb)','FontSize',12)
    title('Molecules Remaining','FontSize',14)
    set(gca,'FontSize',12)
    set(gca,'YScale','log')
    grid on
    legend('VSF','SF','VSF2','VSF2.5')
    
    subplot(2,3,2)
    xlabel('Stage Number','FontSize',12)
    ylabel('Velocity (m/s)','FontSize',12)
    titlestring = sprintf('Alternate Charging Strategies\nSlowing to 50 m/s');
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
    %ylim([0 100])
    title('Phase Space Z','FontSize',14)
    set(gca,'FontSize',12)
    grid on

    subplot(2,3,5)
    xlabel('Time after turn-off (\mus)','FontSize',12)
    ylabel('Population in Laser','FontSize',12)
    title('Time of Flight','FontSize',14)
    set(gca,'FontSize',12)
    %set(gca,'XScale','log')
    xlim([0,50])
    grid on
        
    subplot(2,3,6)
    xlabel('Charging Type','FontSize',12)
    ylabel('Molecules','FontSize',12)
    title('Phase Space Totals','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    %legend('show')
end