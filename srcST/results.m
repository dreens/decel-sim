%% Plots results
function results(rs)
    figure;hold on
    for i=1:length(rs)
        r = rs(i);
        num(i) = r.molnum(end);
        ake(i) = .5*r.mOH*sum(sum(r.vel.^2))/num(i);
        x(i) = r.loadtime * 1e6;
    end
    plot(x,num,'b*')
    xlabel('Stopping Time (\mus)','FontSize',12)
    ylabel('Population','FontSize',12)
    title('Loaded Number v Stopping Time','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    figure
    plot(x,ake,'r*')
    xlabel('Stopping Time (\mus)','FontSize',12)
    ylabel('Average Kinetic Energy','FontSize',12)
    title('Kinetic Energy v Stopping Time','FontSize',14)
    set(gca,'FontSize',12)
    grid on
end