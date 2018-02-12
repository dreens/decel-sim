function molsatmaxinlaser = resultsdetect(rs)
%% Plot molnum v molnumlaser
f1 = figure;

r = rs(1);
molnum = zeros(size(r.molnum));
molnumlaser = zeros(size(r.molnumlaser));
molsatmaxinlaser = [];
for i=1:length(rs)
    molnum = molnum + rs(i).molnum;
    molnumlaser = molnumlaser + rs(i).molnumlaser;
    molsatmaxinlaser = [molsatmaxinlaser; rs(i).molsatmaxinlaser];
end

% Create axes
a1 = axes('Parent',f1,'FontSize',13);

timebegindetect = r.loadingofftime + r.traptime;
locbegindetect = find(abs(r.times-(r.loadingofftime+r.traptime))<1e-8);

xlim(a1,[timebegindetect-1e-4, timebegindetect+3e-4]);
ylim(a1,[0 1.1*molnum(locbegindetect)]);

line([timebegindetect timebegindetect],[0 1.1*molnum(locbegindetect)],...
    'Color',[1 0 0],...
    'LineWidth',1,...
    'LineStyle','--');

box(a1,'on');
grid(a1,'on');
hold(a1,'on');

plot(r.times,molnum,'LineWidth',1);
plot(r.times,molnumlaser,'LineWidth',1);

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

