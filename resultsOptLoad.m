%% Optimize Trapping
% Plot first optimization in tricycle trap with new, sleek load and trap sim
function resultsOptLoad(rs)
    p = sort(unique([rs.phase]));
    l = zeros(size(rs));
    for i=1:length(rs)
        l(i) = rs(i).endphases(end-1);
    end
    l = sort(unique(l));
    l = l(l>250);
    [pp ll] = meshgrid(p,l);
    rr = reshape([rs.numleft],length(p),length(l))';
    
    figure;
    surf(pp,ll,rr)
    
    title('Optimize Speed and Loading On-Time, Tricycle')
    xlabel('Slowing Phase Angle')
    ylabel('Loading ''Phase Angle''')
    zlabel('Final Population')




end