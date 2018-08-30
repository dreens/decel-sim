%% Optimize Trapping
% Plot first optimization in tricycle trap with new, sleek load and trap sim
function resultsOptLoad(rs)
    
    p = zeros(size(rs));
    for i=1:length(rs)
        p(i) = str2double([rs(i).contname(end-1) rs(i).contname(end)]);
    end
    p = sort(unique(p));
    p = 56+p/100;
    %p = sort(unique([rs.phase]));
    l = zeros(size(rs));
    for i=1:length(rs)
        l(i) = rs(i).endphases(end-1);
    end
    l = sort(unique(l));
    l = l(l<40);
    l = l(l~=5e-3);
    [pp ll] = meshgrid(p,l);
    rr = reshape([rs.numleft],length(p),length(l))';
    tt = zeros(size(rs));
    for i=1:length(rs)
        rv = rs(i).vel;
        rv = rv(~isnan(rv));
        tt(i) = sum(rv(:).^2)*rs(i).mOH/(3*rs(i).k*rs(i).numleft);
    end   
    tt = reshape(tt,length(p),length(l))';
    figure;
    %surf(pp,ll,rr)
    bar3(l,rr)
    set(gca,'XTickLabel',p)
    %set(gca,'XTickLabel',p)
    title('Optimize Speed and Loading On-Time, Cryocycle VSF')
    xlabel('Slowing Phase Angle (deg)')
    ylabel('Loading Speed (m/s)')
    zlabel('Final Population')

    figure;
    %surf(pp,ll,rr)
    bar3(l,1e3*tt)
    set(gca,'XTickLabel',p)
    %set(gca,'XTickLabel',p)
    title('Optimize Loading for Temperature, Cryocycle VSF')
    xlabel('Slowing Phase Angle (deg)')
    ylabel('Loading Speed (m/s)')
    zlabel('Temperature (mK)')




end