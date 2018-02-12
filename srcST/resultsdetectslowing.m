function molsatmaxinlaser = resultsdetectslowing(rs)
%% Plot molnum v molnumlaser
f1 = figure;

r = rs(1);

% Create axes
a1 = axes('Parent',f1,'FontSize',13);

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
    timebegindetect = r.decelofftime;
    loc = find(abs(r.times-timebegindetect)<1e-8);

    %plot(r.times,r.molnum,'LineWidth',1);
    times = (r.times(loc:end)-timebegindetect) * 1e6;
    plot(times,r.molnumlaser(loc:end)*1e-3,...
        'LineWidth',2,...
        'DisplayName',sprintf('skip %d, phi=%2.1f',r.skipnum,r.phase));
end

legend('show') 

% Create xlabel
xlabel({'Time (\mus)'},'FontSize',14);

% Create ylabel
ylabel({'Molecule Number (Thousands)'},'FontWeight','bold','FontSize',14);

% Create title
title({'Slowing Argon, 150 stages, v_f=51m/s'},'FontWeight','bold',...
    'FontSize',15);

%{
% Create textbox
annotation(f1,'textbox',...
    [0.211469534050179 0.636820512112657 0.197729991782591 0.0512820519899068],...
    'Color',[1 0 0],...
    'String',{'Trap E-field Off'},...
    'LineStyle','none',...
    'FontSize',14);
%}



