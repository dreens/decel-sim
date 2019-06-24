%% On-axis Energy Diagramming Tool
% plotseqhao.m
% Hao Wu, 4/20/2019

%% Setup the figure
ff = figure;
[ha, pos] = tight_subplot(3,2,[.05 .02],[.15 .08],[.06 .02]);
set(ha([1 2 5 6]),'XTickLabel',{'  0','5','10','15','20   '})
set(ha([3 5]),'YTickLabelMode','auto')
axes(ha(1));plotTweak();
p = get(gca,'Position'); p(2) = .75; set(gca,'Position',p);
axes(ha(2));plotTweak()
p = get(gca,'Position'); p(2) = .75; set(gca,'Position',p);
doF = true;

%% S=1
axes(ha(3))
hold on;

c=sequence1('longdecel',45,225);
d=sequence1('longdecel',0,360);
period=3;
L=5;
for i=1:period

plot(c(1,:)+(i-1)*2*L,c(2,:),'k','LineWidth',2)
plot(d(1,:)+(i-1)*2*L,d(2,:),'k')


plot(c(1,:)+(i-1)*2*L+L,c(2,:),'k','LineWidth',2)
plot(d(1,:)+(i-1)*2*L+L,d(2,:),'k--')


vline=linspace(c(2,1),c(2,end),100);


plot(c(1,end)*ones(100,1)+(i-1)*L,vline,'k:','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L+L,vline,'k:','LineWidth',2)
end
plotTweak()
%% S=3
axes(ha(5))
hold on;

c=sequence1('longdecel',45,360);
d=sequence1('longdecel',0,225);
e=sequence1('longdecel',0,360);
period=2;
L=5;
for i=1:period
col=['k' 'k'];
plot(c(1,:)+(i-1)*3*L,c(2,:),col(2-mod(i,2)),'LineWidth',2)


plot(d(1,:)+(i-1)*3*L+c(1,end)-d(1,1),d(2,:),col(2-mod(i,2)),'LineWidth',2)

plot(e(1,:)+(i-1)*2*L,e(2,:),'k')
plot(e(1,:)+(i)*2*L,e(2,:),'k')

plot(e(1,:)+(i-1)*2*L+L,e(2,:),'k--')
plot(e(1,:)+(i)*2*L+L,e(2,:),'k--')

% plot(c(1,:)+(i-1)*3*L+L,c(2,:),'r','LineWidth',2)
% plot(d(1,:)+(i-1)*3*L,d(2,:),'r','LineWidth',2)
% plot(d(1,:)+(i-1)*2*L+L,d(2,:),'r-.')


vline=linspace(c(2,1),d(2,end),100);

%
plot(d(1,end)*ones(100,1)+(i-1)*3*L+2*L,vline,'k:','LineWidth',2)

end
plotTweak()
%% SF
if ~doF
axes(ha(4))
hold on;

c=sequence1('longdecel',45+90,135+90);
d=sequence1('longdecel',0,360);
e=sequence1('ppgg',45,135);
f=sequence1('ppgg',0,360);
period=4;
L=5;
for i=1:period

plot(c(1,:)+(i-1)*2*L,c(2,:),'k','LineWidth',2)
plot(e(1,:)+(i-1)*2*L,e(2,:),'r','LineWidth',2)
plot(c(1,:)+(i-1)*2*L+L,c(2,:),'k','LineWidth',2)
plot(e(1,:)+(i-1)*2*L+L,e(2,:),'r','LineWidth',2)
plot(d(1,:)+(i-1)*2*L,d(2,:),'k')
plot(f(1,:)+(i-1)*2*L,f(2,:),'r')
plot(d(1,:)+(i-1)*2*L+L,d(2,:),'k--')

vline1=linspace(e(2,end),c(2,1),100);
vline2=linspace(e(2,1),c(2,end),100);

plot(c(1,1)*ones(100,1)+(i-1)*L,vline1,'k:','LineWidth',2)
plot(c(1,1)*ones(100,1)+(i-1)*L+L,vline1,'k:','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L,vline2,'k:','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L+L,vline2,'k:','LineWidth',2)

end
plotTweak()
else
%% F
axes(ha(4))
hold on;

c=sequence1('longdecel',45+90,135+90);
d=sequence1('longdecel',0,360);
e=sequence1('singlerod',225,315);
f=sequence1('singlerod',0,360);
period=4;
L=5;
for i=1:period

plot(c(1,:)+(i-1)*2*L,c(2,:),'k','LineWidth',2)

plot(e(1,:)+(i-1)*2*L,e(2,:),'r','LineWidth',2)

plot(c(1,:)+(i-1)*2*L+L,c(2,:),'k','LineWidth',2)
plot(e(1,:)+(i-1)*2*L+L,e(2,:),'r','LineWidth',2)
plot(d(1,:)+(i-1)*2*L,d(2,:),'k-')
plot(f(1,:)+(i-1)*2*L,f(2,:),'r--')

plot(d(1,:)+(i-1)*2*L+L,d(2,:),'k--')
plot(f(1,:)+(i-1)*2*L+L,f(2,:),'r-')

vline1=linspace(e(2,end),c(2,1),100);
vline2=linspace(e(2,1),c(2,end),100);

plot(c(1,1)*ones(100,1)+(i-1)*L,vline1,'k:','LineWidth',2)
plot(c(1,1)*ones(100,1)+(i)*L,vline1,'k:','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L,vline2,'k:','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i)*L,vline2,'k:','LineWidth',2)
end
plotTweak()
end
%% VSF
axes(ha(6))
hold on;
phi1=-45;
phi2=45;
c=sequence1('longdecel',phi1+180,phi2+180);
d=sequence1('longdecel',0,360);
e=sequence1('ppmm_2mm',phi2,180+phi1);
f=sequence1('ppmm_2mm',0,360);
period=4;
L=5;
for i=1:period

plot(c(1,:)+(i-1)*2*L,c(2,:),'k','LineWidth',2)
plot(e(1,:)+(i-1)*2*L,e(2,:),'r','LineWidth',2)
plot(c(1,:)+(i-1)*2*L+L,c(2,:),'k','LineWidth',2)
plot(e(1,:)+(i-1)*2*L+L,e(2,:),'r','LineWidth',2)
plot(d(1,:)+(i-1)*2*L,d(2,:),'k')
plot(f(1,:)+(i-1)*2*L,f(2,:),'r')


plot(d(1,:)+(i-1)*2*L+L,d(2,:),'k--')

vline1=linspace(e(2,end),c(2,1),100);
vline2=linspace(e(2,1),c(2,end),100);

plot(c(1,1)*ones(100,1)+(i-1)*L,vline1,'k:','LineWidth',2)
plot(c(1,1)*ones(100,1)+(i-1)*L+L,vline1,'k:','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L,vline2,'k:','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L+L,vline2,'k:','LineWidth',2)

end



plotTweak()


%% Convert to PNG and output
ff.Position = [100 50 900 350];
ff.PaperPosition = [0 0 8.8889 3.5];
print(gcf,'2x2-Diagram-Compare','-dpng','-r300')
system('open 2x2-Diagram-Compare.png')



%% Helper Functions
% returns a 2 column array, first column gives z axis coordinates in mm,
% second column gives the potential energy in milliKelvin.
function out=sequence1(ss,phi0,phif)
% 1st colume: x; 2nd colume: 3rd: z; 4th: domain; 5th: E-field (V/m)

load(['Fields/' ss '.mat'],'aenergy');
L=5;
phi=-90:1:90;
z=phi/180*L;
energy=aenergy(-90:10:90);
phi_all=[0:1:360];

if ~any(isnan(energy))
    energy_all=[aenergy(0:1:90) aenergy(89:-1:-89)     aenergy(-90:1:0) ];     
else
    load(['Fields/' ss '.mat'],'vf');
    pa = [90:450];
    zall = pa/180*L*1e-3;
    energy_all = vf(zeros(1,length(zall)),zeros(1,length(zall)),zall);
end

energy_all=flip(energy_all);
amp_z0=interp1(phi_all,energy_all,phi0:1:phif);

out=[(phi0:1:phif)*L/180;amp_z0/3.33e-30/1e5];
end


%% Plot Helper
function plotTweak()
    xlim([5 25])
    ylim([0 120])
    set(gca,'FontSize',13)
    set(gca,'LineWidth',2)
end