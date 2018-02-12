function molsatmaxinlaser = resultsdetectslowingpeaks(rs)
%% Plot molnum v molnumlaser
f1 = figure;

r = rs(1);

% Create axes
a1 = axes('Parent',f1,'FontSize',13);

timebegindetect = r.decelofftime;
locbegindetect = find(abs(r.times-timebegindetect)<1e-8);

xlim(a1,[timebegindetect+0e-4, timebegindetect+7e-4]);
%ylim(a1,[0 1.1*r.molnum(locbegindetect)/3+1]);

% line([timebegindetect timebegindetect],[0 1.1*r.molnum(locbegindetect)],...
%     'Color',[1 0 0],...
%     'LineWidth',1,...
%     'LineStyle','--');

box(a1,'on');
grid(a1,'on');
hold(a1,'on');

for i=1:length(rs)
    r = rs(i);
    %plot(r.times,r.molnum,'LineWidth',1);
    plot(r.times,r.molnumlaser,'LineWidth',1);
    detectpeak(i) = max(r.molnumlaser);
    s12(i) = r.s12;
    s34(i) = r.s34;
end

legend('Begin Detection','Total Number','Detected Number')

% Create xlabel
xlabel({'Time (ms)'},'FontSize',14);

% Create ylabel
ylabel({'Molecule Number'},'FontWeight','bold','FontSize',14);

% Create title
title({'Collection Peak'},'FontWeight','bold',...
    'FontSize',14);

% Create textbox
annotation(f1,'textbox',...
    [0.211469534050179 0.636820512112657 0.197729991782591 0.0512820519899068],...
    'Color',[1 0 0],...
    'String',{'Trap E-field Off'},...
    'LineStyle','none',...
    'FontSize',14);




%% Here I make a contour of the peaks

for i=1:length(rs)
    r = rs(i);
    detectpeak(i) = max(r.molnumlaser);
    s12(i) = r.s12;
    s34(i) = r.s34;
end


s12s = unique(s12);
s34s = unique(s34);

numberlookupsingle = @(d,l) detectpeak(d==s12 & l==s34);
numberlookup = @(d,l) arrayfun(numberlookupsingle,d,l);

minr = round(min(detectpeak),-2);
maxr = round(max(detectpeak),-2);

lines = minr:100:maxr;

[s12ss, s34ss] = meshgrid(s12s,s34s);

hh = figure('position',[100,100,800,800]);
[C, h] = contourf(s12ss,s34ss,numberlookup(s12ss,s34ss),lines);
xlabel('Scaling Rods 1&2','FontSize',12)
ylabel('Scaling Rods 3&4','FontSize',12)
%zlabel(,'FontSize',12)
title('Slowing Peaks, Tune Pin Gap, Timing Fixed','FontSize',14)
a = colorbar;
ylabel(a,'Peak in Slowing','FontSize',12)
set(a,'FontSize',12);
set(gca,'FontSize',12);

set(a,'Limits',[minr maxr]);
set(gca,'CLim',[minr maxr]);
pause(1)
set(hh,'Position',[100,100,800,800])
pause(1)


