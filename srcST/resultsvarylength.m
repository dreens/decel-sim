function resultsvarylength(rs)
%% Trapped Number versus Decelerator Length

rs = rs(1:1:end);

numstages = zeros(size(rs));
finalnum = zeros(size(rs));
phases = zeros(size(rs));
finalvz = zeros(size(rs));

for i=1:length(rs)
    numstages(i) = rs(i).stages;
    finalnum(i) = rs(i).molnum(end);
    phases(i) = rs(i).phase;
    finalvz(i) = rs(i).finalvz;
end

nums = unique(numstages);
fvzs = unique(finalvz);

numberlookupsingle = @(n,f) finalnum(n==numstages & f==finalvz);
numberlookup = @(n,f) arrayfun(numberlookupsingle,n,f);

[numg,fvg] = meshgrid(nums,fvzs);

figure(4)
subplot(2,1,1)
hold on
plot(nums',numberlookup(numg,fvg)')
%xlabel('Decel Length (Stages)')
ylabel('Molecules')
title('Molecules v Decelerator Length')


subplot(2,1,2)
hold on
plot(nums,phases)
xlabel('Decel Length (Stages)')
ylabel('Phase Angle')
title('Phase Angle v Decelerator Length')
hold on
%plot(nums,mean(numberlookup(numg,fvg)))
%plot(nums,min(numberlookup(numg,fvg)))


%figure
%plot(numstages,phases)