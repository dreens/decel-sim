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
rot = [0 0 90 90 180 180 270 270];
rot = repmat(rot,1,ceil(n/4));
rots = rot(1:2*n);

tran = [1 1 0 0];
tran = repmat(tran,1,ceil(n/2));
trans = tran(1:2*n);

p = 57.415;
ep = repmat([180-p 180+p 360-p p],1,n); % S=1

%s1test = simdecel('initvz',ivz,'decels',dcs1,'chargetype',ct,'rot',rots,'turnon',1e-12,...
%    'trans',trans,'endphases',ep,'calctype',repmat('p',1,1000),'num',1);

tt = s1test.times;
m = num2cell((2.82328e-26)*12.5./(10:.25:13));
s1 = simdecel('initvz',ivz,'decels',dcs1,'chargetype',ct,'rot',rots,'everyother',true,...
    'trans',trans,'endphases',tt,'calctype',repmat('t',1,1000),'num',3000000,'mOH',m);
sf = simdecel('initvz',ivz,'decels',dcsf,'chargetype',ct,'rot',rots,'everyother',false,...
    'trans',trans,'endphases',tt,'calctype',repmat('t',1,1000),'num',3000000,'mOH',m);

%% plot it up
s1 = resultsdecel(s1);

sf = resultsdecel(sf);

figure
voltages = 10:.25:13;
plot(voltages,[s1.tofpeak],'b-')
hold on
plot(voltages,[sf.tofpeak],'r-')
legend('S=1','F')
title('Corrected phase variation issue')
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

