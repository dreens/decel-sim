%% Guide Escape
% Here I wish to study the escape dynamics from the 2D trap generated by 
% Decelerator Guiding potentials.

% Get fine guiding grid
guide = efftrap3D('ppgg','ppgg',-180,0,180);
guide = mean(guide,3);
x = -.975:.025:.975;
x = x*1e-3;
sp = 1e-2;
xq = -1:sp:1;
xq = xq*1e-3;
[x, y] = meshgrid(x,x);
[xq, yq] = meshgrid(xq,xq);
guideq = interp2(x,y,guide,xq,yq,'spline');

% Get derivatives for acceleration
dVdxm = conv2(guideq,[-.5 0 .5],'same')/(sp*1e-3);
dvdx  = griddedInterpolant(xq',yq',dVdxm,'linear','none');

% Run the odesolver
m = 2.82328e-26; %binding energy accounted for.
times = 0:1e-6:5e-3;
escO = odeset('Events',@(t,y) deal([isnan(y(1))-.5,1,0]));
[t, p] = ode45(@(t,p) [p(3);p(4);dvdx(p(2),p(1))/m;dvdx(p(1),p(2))/m],...
    [0 5e-3],[0 .0009 -19 0],escO);
length(t)
figure; plot(p(:,1),p(:,2))


%% Run odesolver for collection of related initial positions, velocities
% added some loops to automatically get ratio of ratios and check a few
% scalings of the force fields.
allr = zeros(11,1);
rs = 0:.1111:1;
for j=1:length(rs)
temps = 0.3;%[0.3 0.6];
%for tt = 1:2
N = 2000;
ri =  sqrt(rand(N,1))*0.5e-3;
thi = rand(N,1)*2*pi;
k = 1.3806e-23;
T = temps(1);
y0 = [ri.*cos(thi), ri.*sin(thi) , randn(N,2)*sqrt(k*T/m)];
%y0 = [rand(N,2)*2e-3-1e-3 , randn(N,2)*sqrt(k*T/m)];
times = zeros(N,1);
ke = zeros(N,1);
pe = zeros(N,1);
pp = zeros(N,4);
r = rs(j);
fprintf('\n     ');
for i=1:N
    fprintf('\b\b\b\b\b\b\b%7d',i)
    [t, p] = ode45(@(t,p) [p(3);p(4);r*dvdx(p(2),p(1))/m;r*dvdx(p(1),p(2))/m],...
        [0 2.05e-3],y0(i,:),escO);
    times(i) = t(end);
    pp(i,:) = p(end,:);
end
%figure;hist(times,50)
%s40(tt) = sum(times>.00025);
allr(j) = sum(times==max(times));
%end
%allr = 100*(s40(2)*s333(1)/(s40(1)*s333(2))-1)
fprintf('\n     ')
end
allr
%% Results Readme
% So far I have checked the escape dynamics for molecules loaded into the
% transverse trap created by DC guiding with +/-6.25 kV applied in a ++--
% manner to the rods, with different initial condition assumptions:
% 0.5 K, filled 2x2mm square: 1243 @ 40 stages, 1206 @ 333 stages.
% 0.3 K, filled 2x2mm square: 1317 @ 40 stages, 1273 @ 333 stages.
% The ratio of ratios gives a 0.4% effect, (experiment shows 10%).
%
% 0.5 K, filled 1mm diameter circle: 1244 @ 40 stages, 1229 @ 333 stages.
% 0.3 K, filled 1mm diameter circle: 2017 @ 40 stages, 1982 @ 333 stages.
% ratio of ratios: 0.5% effect.
%
% 1.0 K, filled 2x2mm square: 964 @ 40 stages, 940 @ 333 stages.
% 1.3 K, filled 2x2mm square: 839 @ 40 stages, 816 @ 333 stages.
% ratio of ratios: 5% effect.
%
% 1.0 K, offset 1mm circle: 770 @ 40 stages, 751 @ 333 stages.
% 1.3 K, offset 1mm circle: 633 @ 40 stages, 614 @ 333 stages.
% ratio of ratios: 0.5% effect.

%% Plot what's going on
figure(100)
hist(times,[2.5e-5:5e-5:2.025e-3])
hold on
times10 = repmat(times,1,10);
times10=times10(times10>1.5e-4 & times10<2.5e-4);
%hist(times10,[1.75e-4:5e-5:2.75e-4])
times100 = repmat(times,1,1000);
times100=times100(times100>2.5e-4 & times100<2e-3);
hist(times100,[2.75e-4:5e-5:1.975e-3])

%% Figure out how temperature is evolving.
KE = sum(pp(:,[3 4]).^2,2)*m*.5;
PE = interp2(xq,yq,guideq,pp(:,1),pp(:,2));
NN = ~isnan(pp(:,4));
T = mean(KE(NN))/k
V = mean(PE(NN))/k
E = T+V

% Ti=.05K, 1mm circle, 28, 40, 68  mK (T,V,E)
% Ti=0.3K, 1mm circle, 37, 56, 93  mK (T,V,E)
% Ti=0.5K, 1mm circle, 39, 58, 97  mK (T,V,E)
% Ti=1.0K, 1mm circle, 39, 60, 98  mK (T,V,E)
% Ti=.05K, 2mm square, 55, 59, 114 mK (T,V,E)
% Ti=0.3K, 2mm square, 59, 71, 130 mK (T,V,E)
% Ti=0.5K, 2mm square, 61, 71, 132 mK (T,V,E)
% Ti=1.0K, 2mm square, 58, 74, 132 mK (T,V,E)
% Ti=.05K, 1mm circle offset, 41, 50, 90  mK (T,V,E)
% Ti=0.3K, 1mm circle offset, 50, 62, 112 mK (T,V,E)
% Ti=0.5K, 1mm circle offset, 51, 64, 115 mK (T,V,E)
% Ti=1.0K, 1mm circle offset, 53, 64, 117 mK (T,V,E)

%% Put the above results in a figure
% First of all, for thermal equilibrium, we should have kT of thermal, 2kT
% of potential. Not surprisingly, we don't find this, since there's no
% mechanism for thermal equilibrium, we're just exchanging potential and
% kinetic by orbiting in the trap.
tp05 = [68 114 90];
tp3 = [93 130 112];
tp5 = [97 132 115];
tp10 =[98 132 117];
all = [tp05;tp3;tp5;tp10];
all = all(:,[1 3 2]);
x = [1 2 3];
figure; hold on;
colors = get(gca,'ColorOrder');
bar(all')
ll = legend('50 mK','300 mK','500 mK','1000 mK');
set(ll,'Location','NorthWest');
ylabel('Total Energy (mK)')
xlabel('Initialization Strategy')
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',{'','1mm circle','','offset circle','','2mm square',''})
set(gca,'FontSize',12);
title('Temperature after Guiding, Different Initializations')
grid on
