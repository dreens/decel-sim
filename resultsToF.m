%% Plots results
function resultsToF(rs)
    figure
    hold on
    r = rs(1);
    colors = get(gca,'ColorOrder');
    colors = [colors ; colors ; colors ; colors];
    diamL = 2.5e-3;
    distL = 8e-3 + r.pos(1,3);
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
    end
    xlabel('Time \mus')
    ylabel('Population (arb)')
    title('Time of Flight of 50 m/s Molecules')
    legend('11 mm','22 mm','33 mm')
    grid on
end