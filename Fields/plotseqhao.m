%%s=1
figure(200);
subplot(4,1,1)
hold on;

c=sequence1('longdecel',45,225);
d=sequence1('longdecel',0,360);
period=3;
L=5;
for i=1:period

plot(c(1,:)+(i-1)*2*L,c(2,:),'k','LineWidth',2)
plot(d(1,:)+(i-1)*2*L,d(2,:),'k-.')


plot(c(1,:)+(i-1)*2*L+L,c(2,:),'r','LineWidth',2)
plot(d(1,:)+(i-1)*2*L+L,d(2,:),'r-.')


vline=linspace(c(2,1),c(2,end),100);


plot(c(1,end)*ones(100,1)+(i-1)*L,vline,'b','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L+L,vline,'b','LineWidth',2)

end

%% s=3
figure(200);
subplot(4,1,2)
hold on;

c=sequence1('longdecel',45,360);
d=sequence1('longdecel',0,225);
e=sequence1('longdecel',0,360);
period=2;
L=5;
for i=1:period
col=['k' 'r'];
plot(c(1,:)+(i-1)*3*L,c(2,:),col(2-mod(i,2)),'LineWidth',2)


plot(d(1,:)+(i-1)*3*L+c(1,end)-d(1,1),d(2,:),col(2-mod(i,2)),'LineWidth',2)

plot(e(1,:)+(i-1)*2*L,e(2,:),'k-.')
plot(e(1,:)+(i)*2*L,e(2,:),'k-.')

plot(e(1,:)+(i-1)*2*L+L,e(2,:),'r-.')
plot(e(1,:)+(i)*2*L+L,e(2,:),'r-.')

% plot(c(1,:)+(i-1)*3*L+L,c(2,:),'r','LineWidth',2)
% plot(d(1,:)+(i-1)*3*L,d(2,:),'r','LineWidth',2)
% plot(d(1,:)+(i-1)*2*L+L,d(2,:),'r-.')


vline=linspace(c(2,1),d(2,end),100);

%
plot(d(1,end)*ones(100,1)+(i-1)*3*L+2*L,vline,'b','LineWidth',2)

end
%% VSF

figure(200);
subplot(4,1,4)
hold on;

c=sequence1('longdecel',45+90,135+90);
d=sequence1('longdecel',0,360);
e=sequence1('ppgg',45,135);
f=sequence1('ppgg',0,360);
period=3;
L=5;
for i=1:period

plot(c(1,:)+(i-1)*2*L,c(2,:),'k','LineWidth',2)
plot(e(1,:)+(i-1)*2*L,e(2,:),'c','LineWidth',2)
plot(c(1,:)+(i-1)*2*L+L,c(2,:),'r','LineWidth',2)
plot(e(1,:)+(i-1)*2*L+L,e(2,:),'c','LineWidth',2)
plot(d(1,:)+(i-1)*2*L,d(2,:),'k-.')
plot(f(1,:)+(i-1)*2*L,f(2,:),'c-.')


plot(d(1,:)+(i-1)*2*L+L,d(2,:),'r-.')

vline1=linspace(e(2,end),c(2,1),100);
vline2=linspace(e(2,1),c(2,end),100);

plot(c(1,1)*ones(100,1)+(i-1)*L,vline1,'b','LineWidth',2)
plot(c(1,1)*ones(100,1)+(i-1)*L+L,vline1,'b','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L,vline2,'b','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L+L,vline2,'b','LineWidth',2)

end

%% SF

figure(200);
subplot(4,1,3)
hold on;

c=sequence1('longdecel',45+90,135+90);
d=sequence1('longdecel',0,360);
e=sequence2('singlerod2',45,135);
f=sequence2('singlerod2',0,180);
period=3;
L=5;
for i=1:period

plot(c(1,:)+(i-1)*2*L,c(2,:),'k','LineWidth',2)

plot(e(1,:)+(i-1)*2*L,e(2,:),'c','LineWidth',2)

plot(c(1,:)+(i-1)*2*L+L,c(2,:),'r','LineWidth',2)
plot(e(1,:)+(i-1)*2*L+L,e(2,:),'m','LineWidth',2)
plot(d(1,:)+(i-1)*2*L,d(2,:),'k-.')
plot(f(1,:)+(i-1)*2*L,f(2,:),'c-.')

plot(f(1,:)+(i-1)*2*L+L,f(2,:),'m-.')
plot(d(1,:)+(i-1)*2*L+L,d(2,:),'r-.')

vline1=linspace(e(2,end),c(2,1),100);
vline2=linspace(e(2,1),c(2,end),100);

plot(c(1,1)*ones(100,1)+(i-1)*L,vline1,'b','LineWidth',2)
plot(c(1,1)*ones(100,1)+(i)*L,vline1,'b','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L,vline2,'b','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i)*L,vline2,'b','LineWidth',2)
end
%% XSF
figure(201);
subplot(4,1,3)
hold on;
phi1=-25;
phi2=45;
c=sequence1('longdecel',phi1+180,phi2+180);
d=sequence1('longdecel',0,360);
e=sequence2('pmpm_2mm',phi2,180+phi1);
f=sequence2('pmpm_2mm',0,360);
period=3;
L=5;
for i=1:period

plot(c(1,:)+(i-1)*2*L,c(2,:),'k','LineWidth',2)
plot(e(1,:)+(i-1)*2*L,-e(2,:),'c','LineWidth',2)
plot(c(1,:)+(i-1)*2*L+L,c(2,:),'r','LineWidth',2)
plot(e(1,:)+(i-1)*2*L+L,-e(2,:),'c','LineWidth',2)
plot(d(1,:)+(i-1)*2*L,d(2,:),'k-.')
plot(f(1,:)+(i-1)*2*L,-f(2,:),'c-.')


plot(d(1,:)+(i-1)*2*L+L,d(2,:),'r-.')

vline1=linspace(-e(2,end),c(2,1),100);
vline2=linspace(-e(2,1),c(2,end),100);

plot(c(1,1)*ones(100,1)+(i-1)*L,vline1,'b','LineWidth',2)
plot(c(1,1)*ones(100,1)+(i-1)*L+L,vline1,'b','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L,vline2,'b','LineWidth',2)
plot(c(1,end)*ones(100,1)+(i-1)*L+L,vline2,'b','LineWidth',2)

end




%% Helper Functions
function out=sequence1(ss,phi0,phif)
% 1st colume: x; 2nd colume: 3rd: z; 4th: domain; 5th: E-field (V/m)

load(['Fields/' ss '.mat'],'aenergy');
L=5;
phi=-90:1:90;
z=phi/180*L;
energy=aenergy(-90:1:90);

phi_all=[0:1:360];
energy_all=[aenergy(0:1:90) aenergy(89:-1:-89)     aenergy(-90:1:0) ];     

energy_all=flip(energy_all);
amp_z0=interp1(phi_all,energy_all,phi0:1:phif);

out=[(phi0:1:phif)*L/180;amp_z0];
end
                  
function out=sequence2(ss,phi0,phif)
% 1st colume: x; 2nd colume: 3rd: z; 4th: domain; 5th: E-field (V/m)

load(['Fields/' ss '.mat'],'aenergy');
L=5;

for phi=-90:1:0
z=phi/180*L;

energy(phi+91)=aenergy(phi);
end
phi_all=[0:1:360];

energy_all=[flip(energy) energy(2:end-1) flip(energy) energy(2:end)];
amp_z0=interp1(phi_all,energy_all,phi0:1:phif);

out=[(phi0:1:phif)*L/180;amp_z0];
                  
end