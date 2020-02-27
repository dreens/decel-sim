%% Low Velocity Breakdown
%
% In this script, I vary the voltage of the decelrator.
%

%% First for S=1, F Mode, to save times
ivz = 825.1;
dcs1 = struct('a','longdecel','b','longdecel');
dcf = struct('a','longdecel','b','singlerod');
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

ep = @(p) repmat([180-p 180+p 360-p p],1,n); % S=1

vf = num2cell(800:-25:25);

n = 5000000;

f = simdecel('initvz',ivz,'decels',dcf,'chargetype',ct,'rot',rots,'finalvz',vf,...
    'trans',trans,'endphases',ep,'calctype',repmat('p',1,1000),'num',n,...
    'tempz',2,'tempxy',1,'spreadz',10e-3,'spreadxy',2e-3,'dist','flat','smallt',1e-7);

s1 = simdecel('initvz',ivz,'decels',dcs1,'chargetype',ct,'rot',rots,'finalvz',vf,...
    'trans',trans,'endphases',ep,'calctype',repmat('p',1,1000),'num',n,...
    'tempz',2,'tempxy',1,'spreadz',10e-3,'spreadxy',2e-3,'dist','flat','smallt',1e-7);
%% plot it up
s1 = resultsTOFprocess(s1);

f = resultsTOFprocess(f);

for i=1:length(f)
    f(i).finalvel = f(i).vels(end);
end
%%
figure
plot([f.finalvel],[s1.tofpeak],'b-')
hold on
plot([f.finalvel],[f.tofpeak],'r-')
legend('S=1','F')
title('Vary Final Speed Sim')
xlabel('Final Speed (m/s)')
ylabel('Population (ToF Peak)')

