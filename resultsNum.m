%% Plots results
function resultsNum(rs)
    figure
    hold on
    nums = [rs.numleft];
    temps = [rs.tempxy];
    decels = {rs.decel};
    pmpm = strcmp(decels,'pmpm_2mm_no-sym');
    numpm = nums(pmpm);
    numpp = nums(~pmpm);
    ratio = numpp./numpm;
    errorpercent = sqrt(1./numpm+1./numpp);
    err = ratio.*errorpercent;
    figure;
    errorbar(temps(pmpm),ratio,err)

    xlabel('Temperature (^\circ K)')
    ylabel('Ratio (arb)')
    set(gca,'XScale','log')
    xlim([0.08 15])
    title('Guiding +-+- v ++--, Symmetry Corrected')
    legend('6.5 kV, 10^5 Molecules')
    grid on
end