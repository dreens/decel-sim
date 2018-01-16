function resultsvarylength(rs)
%% Trapped Number versus Decelerator Length

numstages = zeros(length(rs));
finalnum = zeros(length(rs));
phases = zeros(length(rs));

for i=1:length(rs)
    numstages(i) = rs(i).stages;
    finalnum(i) = rs(i).molnum(end);
    phases(i) = rs(i).phase;
end


figure
plot(numstages,finalnum)

figure
plot(numstages,phases)