%% Plots results
function resultsdecel(rs)
    figure('Position',[0,0,2000,2000])
    colors = get(gca,'ColorOrder');
    colors = [colors ; (1-colors)];
    lines = repmat({'-','-'},1,10);
    markers = repmat({'o','x','sq'},1,10);
    sname = repmat({'Triple','DSwitch','Normal'},1,6);
    for i=1:length(rs)
        r = rs(i);
        i2 = i;%round(i/2+.25);
        c = colors(i2,:);
        l = lines{i};
        m = markers{i};
        subplot(2,3,1); hold on
        if r.triplemode
            r.numstage = r.numstage/3;
            r.molnum = r.molnum(1:3:end);
        end
        h = plot(r.molnum/r.molnum(1),'Color',c,'LineStyle',l);
        h.DisplayName = sprintf('v=%2.1f, %s',r.vels(end),sname{i});
        g = plot(r.numstage,r.molnum(end)/r.molnum(1),'Color',c,'Marker',m);
        g.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        subplot(2,3,2); hold on
        if r.triplemode
            r.stages = 1:(length(r.stages)/3);
            r.vels = r.vels(1:3:end);
        end
        plot(0:max(r.stages),[r.initvz r.vels],'b-','Color',c,'LineStyle',l);
        
        subplot(2,3,3); hold on
        pos = r.pos; vel = r.vel;
        pos(:,3) = (pos(:,3)-pos(1,3))/2;
        vel(:,3) = vel(:,3) - vel(1,3);
        pos = pos/5e-4;
        vel = vel/5;
        r.PSN = pos(:,1).^2+vel(:,1).^2 < 1;
        r.PSN = r.PSN & (pos(:,2).^2+vel(:,2).^2 < 1);
        r.PSN = r.PSN & (pos(:,3).^2+vel(:,3).^2 < 1);
        r.PSN = sum(r.PSN);
        vol = pi^3*0.5^2*1*5^3;
        r.PSD = r.PSN / vol; % in [mm*m/s]^-3
        plot(r.pos(:,2)*1e3,r.vel(:,2),'b.','Color',c,'MarkerSize',3,...
            'DisplayName',sprintf('PSD=%0.1e',r.PSD));
        
        subplot(2,3,4); hold on
        xx = 3*mod(i-1,2)-1.5;
        yy = 30*(i2-2);
        plot((r.pos(:,3)-r.pos(1,3))*1e3+xx*0,r.vel(:,3)+yy-255,'b.','Color',c,'MarkerSize',3);
        
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
    xax = 1:length(rs);%[rs(1:2:end).numstage];
    plot(xax,nls(1:1:end),'DisplayName','Total','Marker','x','LineStyle','-');
    plot(xax,tps(1:1:end),'DisplayName','Peak','Marker','o','LineStyle','-');
    plot(xax,tpa(1:1:end)*1e5,'DisplayName','0.1*Area [N*us]','Marker','sq','LineStyle','-');
    legend('show')
    legend('Location','South')
    
    subplot(2,3,1)
    xlabel('Stage Number','FontSize',12)
    ylabel('Population (arb)','FontSize',12)
    title('Molecules Remaining','FontSize',14)
    set(gca,'FontSize',12)
    set(gca,'YScale','log')
    legend('Advanced','Delay','Normal')
    grid on
    
    subplot(2,3,2)
    xlabel('Stage Number','FontSize',12)
    ylabel('Velocity (m/s)','FontSize',12)
    titlestring = sprintf('%s\n%s, Phi=%2.1f',...
        'Advanced v Delayed Switch Modes',...
        'Slowing to 255 m/s in 111 stages', rs(1).phase);
    title(titlestring,'FontSize',14)
    set(gca,'FontSize',12)
    grid on

    subplot(2,3,3)
    xlabel('Y Position (mm)','FontSize',12)
    ylabel('Y Velocity (m/s)','FontSize',12)
    title('Phase Space Y','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    legend('show')

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
    xlim([10,50])
    grid on
        
    subplot(2,3,6)
    xlabel('Stage Number','FontSize',12)
    ylabel('Molecules','FontSize',12)
    title('Phase Space Totals','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    %legend('show')
end