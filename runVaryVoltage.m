%% Low Velocity Breakdown
%
% In this script, I vary the voltage of the decelrator.
%

%% First for S=1, F Mode, to save times
ivz = 825;
dcs1 = struct('a','longdecel','b','longdecel');
dcsf = struct('a','longdecel','b','singlerod');
ct = ivz; rots = ivz; trans = ivz; ep = ivz;
n = 333;
ct = repmat('ba',1,n);
ct(1) = 'a';
rot = [0 0 90 90 180 180 270 270];
rot = repmat(rot,1,ceil(n/4));
rots = rot(1:2*n);

tran = [1 1 0 0];
tran = repmat(tran,1,ceil(n/2));
trans = tran(1:2*n);

p = 50;
ep = @(p) repmat([180-p 180+p 360-p p],1,n); % S=1

% get the masses
voltages = 9:.25:14;
m = num2cell((2.82328e-26)*12.5./voltages);
voltage_flat = [ones(1,5)*9.85 10:.2:14];
m1 = num2cell((2.82328e-26)*12.5./voltage_flat);

% get the times
%s1t = simdecel('initvz',ivz,'decels',dcs1,'chargetype',ct,'rot',rots,'turnon',1e-16,...
%    'trans',trans,'endphases',ep,'calctype',repmat('p',1,1000),'num',1,'mOH',m1);
%sft = simdecel('initvz',ivz,'decels',dcsf,'chargetype',ct,'rot',rots,'turnon',1e-16,...
%    'trans',trans,'endphases',ep,'calctype',repmat('p',1,1000),'num',1,'mOH',m1);
%tts1 = {};
%ttsf = {};
%for i=1:length(voltages)
%    tts1{i} = @(x) s1t(i).times;
%    ttsf{i} = @(x) sft(i).times;
%end

s1t3 = simdecel('initvz',ivz,'decels',dcs1,'chargetype',ct,'rot',rots,'turnon',300e-9,...
    'trans',trans,'endphases',tts1,'calctype',repmat('t',1,1000),'num',1000000,'mOH',m,...
    'tempz',5,'tempxy',5,'spreadz',10e-3,'spreadxy',3e-3,'dist','sphere','smallt',2e-7);
sft3 = simdecel('initvz',ivz,'decels',dcsf,'chargetype',ct,'rot',rots,'turnon',300e-9,...
    'trans',trans,'endphases',ttsf,'calctype',repmat('t',1,1000),'num',1000000,'mOH',m,...
    'tempz',5,'tempxy',5,'spreadz',10e-3,'spreadxy',3e-3,'dist','sphere','smallt',2e-7);

%% plot it up
s1t3 = resultsTOFprocess(s1t3);

sft3 = resultsTOFprocess(sft3);

s = 11.2/10.5;

figure
plot(voltages*s,[s1t3.tofpeak],'b-')
hold on
plot(voltages*s,[sft3.tofpeak],'r-')
legend('S=1','F')
title('Vary Voltage, Full Initial Distribution.')
xlim([9.9 13.1])
xlabel('Voltage (kV)')
ylabel('Population (ToF Peak)')

%% Next for VSF:

%{
%% Plotting
%
% Even for pathetically low speeds, I don't find much loss. I need to
% include some flight time out of the decelerator, which is universally
% required for use of the beam either for collisions or trapping.
%
% I'll start with say 5 mm flight.

l = length(rss1);
s1 = [rss1.numleft];
s1v = s1;
dd = 10e-3;
for i=1:l
    s = rss1(i);
    s1v(i) = s.vels(end);
    tt = dd/s1v(i);
    pf = s.pos + s.vel*tt;
    pf(:,3) = (pf(:,3) - pf(1,3))/2;
    pf = pf/1e-3;
    pf = sum(pf.^2,2);
    s1(i) = sum(pf<1.5);
end

figure; hold on
plot(s1v,s1,'DisplayName','S=1','LineWidth',2)
xlabel('Final Speed','FontSize',13)
ylabel('Number Remaining','FontSize',13)
title('Breakdown of Effective Trap','FontSize',14)
set(gca,'FontSize',13)
h = legend('show');
set(h,'FontSize',13)

l = length(rssf);
sf = [rssf.numleft];
sfv = sf;
dd = 10e-3;
for i=1:l
    s = rssf(i);
    sfv(i) = s.vels(end);
    tt = dd/sfv(i);
    pf = s.pos + s.vel*tt;
    pf(:,3) = (pf(:,3) - pf(1,3))/2;
    pf = pf/1e-3;
    pf = sum(pf.^2,2);
    sf(i) = sum(pf<1.5);
end

plot(sfv,sf,'DisplayName','SF','LineWidth',2)

%}

