function results_stage_gap(r)
%% Plot the Final Molecule Number v Stage Gap

lr = length(r);
molnums = zeros(1,lr);
stagegaps = zeros(1,lr);

for i = 1:lr
    molnums(i) = r(i).molnum(end) - sum(abs(r(i).pos(:,3) - r(i).pos(1,3) ) > 3e-3);
    stagegaps(i) = 2+0.1*str2num(r(i).decel(end));

end
figure
plot(stagegaps,molnums)
xlabel('Interstage Gap','FontSize',12)
ylabel('Final Number','FontSize',12)
title('Choosing Stage Gap','FontSize',14)


end