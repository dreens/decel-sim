%% Ring Decelerator Field Processor
%
cd ..
allpots = cell(1,5);
for ii=0:4

pot = importdata(['Ring Fields/ringdecel_efftrap' num2str(ii) '.dat'],' ',9);
all = pot.data;

r = all(:,1)*1e-3; %convert to meters.
z = all(:,2)*1e-3;
e = all(:,3);   %units of V/m

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

allpots{ii+1} = vv;

end
%% Look things over
for ii=1:5
    figure;
    surf(rr,zz,allpots{ii});
    title(['Traveling Wave, \phi=' num2str(ii-1) '\pi/20'])
    
    
end


%% Take the average
