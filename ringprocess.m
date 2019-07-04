%% Ring Decelerator Field Processor
%
allpots = cell(1,5);
for ii=5:-1:0

pot = importdata(['Ring Fields/ringdecel_efftrap' num2str(ii) '.dat'],' ',9);
all = pot.data;

r = all(:,1)*1e-3; %convert to meters.
z = all(:,2)*1e-3;
e = all(:,3)*(5/8);   %units of V/m % added 5/8 since trav wave devices only go to rest at 10 kV peak-peak thus far.

% find unique x,y,z coordinates
rs = sort(uniquetol(r,1e-6,'DataScale',1));
zs = sort(uniquetol(z,1e-6,'DataScale',1));

% get the datapoint spacing, used for taking derivatives later.
rsp = uniquetol(diff(rs));
zsp = uniquetol(diff(zs));

% create lookup functions that tell you 'n' for a given x value such
% that x is the nth x value.
r2i = @(rr) arrayfun(@(r) find(r==rs),rr);
z2i = @(zz) arrayfun(@(z) find(z==zs),zz);

% the size of the 3D data matrices to be filled
fullsize = [length(rs), length(zs)];

% for each x,y,z value in the COMSOL loaded x,y,z columns, check which
% linear index this corresponds to in a 3D matrix of size fullsize.
locs = sub2ind(fullsize,r2i(r),z2i(z));

ee = zeros(fullsize);
rr = zeros(fullsize);
zz = zeros(fullsize);
ee(locs) = e;
rr(locs) = r;
zz(locs) = z;

% create an anonymous function that can give OH doubly stretched state
% potential energy as a function of bfield, efield, and ebangle (t).
last = @(x) x(end);
energysingle = @(e) last(sort(eig(OH_Ham_Simple_SI(0,e,0))));
energy = @(e)  arrayfun(energysingle,e);

% get the potential energy
vv = energy(ee);

% convert for the 1 1 state of D2O.
aa = 200; b = 0.45; de = 3.3e-30; h = 6.626e-34; cc = 3e10;
fff = @(j) h*cc*b*( sqrt(aa^2 + (1e-5*j/de).^2) - aa ) / ( sqrt(aa^2 + 150^2) - aa );
vv = fff(vv);

allpots{ii+1} = vv;
end
%% Look things over
if false
for ii=1:6
    figure;
    surf(rr,zz,allpots{ii});
    title(['Traveling Wave, \phi=' num2str(ii-1) '\pi/24'])
    
end
end
%% Take the average
meanpot = zeros(size(allpots{1}));

for ii=1:6
    meanpot = meanpot + allpots{ii};
end
meanpot = meanpot/6;
    
if false
figure;
surf(rr,zz,meanpot);
title(['Traveling Wave Averaged'])
end
    
    
%% Consider different tilts
accs = 0:.25:3;
deps = [];
pot3D = zeros(2*size(meanpot,1)-1, size(meanpot,2), 3);
mOH = 2.82328e-26; % Accounts for Oxygen binding energy
mid = size(meanpot,1);
raccs = 0:.25:3;
tpots = cell(size(raccs));

for a=accs
    thisTiltPot = meanpot + zz * a * 1e3 * mOH;
    thisTiltPot = thisTiltPot - thisTiltPot(1,31);
    pot3D(mid:end,:,2) = thisTiltPot;
    pot3D(mid:-1:1,:,2) = thisTiltPot;
    pot3D(:,:,[1 3]) = max(thisTiltPot(:));
    deps(end+1) = effTrapMinDepth(pot3D);
end
for i=1:length(raccs)
    thisTiltPot = meanpot + zz * raccs(i) * 1e3 * mOH;
    thisTiltPot = thisTiltPot - thisTiltPot(1,31);
    tpots{i} = thisTiltPot;
end

    
%% Plot the ring decel trap depth v decel
%figure;
depsmK = deps / 1.38e-23 * 1e3;
plot(accs,depsmK,'DisplayName','TW','LineWidth',2)
%xlabel('Deceleration (km/s/s)','FontSize',13)
%ylabel('Worst Case Trap Depth (mK)','FontSize',13)
%title('Trap Depth v Deceleration, 10 kV (pp) T-Wave','FontSize',14)
%set(gca,'FontSize',13)




