%% Low Velocity Breakdown
%
% In this script, I vary the voltage of the decelrator.
%

%% First for S=1, F Mode, to save times
ivz = 825;
dcs1 = struct('a','longdecel','b','longdecel');
dcf = struct('a','longdecel','b','singlerod');
ct = ivz; rots = ivz; trans = ivz; ep = ivz;
n = 333;
ct = repmat('ba',1,n);
rot = [0 0 90 90 180 180 270 270];
rot = repmat(rot,1,ceil(n/4));
rots = rot(1:2*n);

tran = [1 1 0 0];
tran = repmat(tran,1,ceil(n/2));
trans = tran(1:2*n);

velocities = 50:20:820;
phases = (825^2 - velocities.^2)/(825^2-50^2)*57.64;
tt = {};

for p=phases
    ep = repmat([180-p 180+p 360-p p],1,n); % S=1
    s = simdecel('initvz',ivz,'decels',dcs1,'chargetype',ct,'rot',rots,'turnon',1e-12,...
        'trans',trans,'endphases',ep,'calctype',repmat('p',1,1000),'num',1);
    tt{end+1} = s.times;
end

s1 = simdecel('initvz',ivz,'decels',dcs1,'chargetype',ct,'rot',rots,'everyother',true,...
    'trans',trans,'endphases',tt,'calctype',repmat('t',1,1000),'num',1000000);
f = simdecel('initvz',ivz,'decels',dcf,'chargetype',ct,'rot',rots,'everyother',false,...
    'trans',trans,'endphases',tt,'calctype',repmat('t',1,1000),'num',1000000);

%% plot it up
s1 = resultsTOFprocess(s1);

f = resultsTOFprocess(f);

for i=1:length(f)
    f(i).finalvel = f(i).vels(end);
end

figure
plot([f.finalvel],[f.tofpeak]./[s1.tofpeak],'b-')
hold on
legend('Enhancement Ratio')
title('Vary Final Speed Sim')
xlabel('Final Speed (m/s)')
ylabel('Population (ToF Peak)')

