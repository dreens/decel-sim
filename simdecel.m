function rsf = simdecel()
%MATLAB Simulation of OH experiment
%   I had a classdef version of this before but it was brutally slow.
%   turns out the slowness wasn't due to it being a class, but I'm going to
%   stick with this struct implementation nonetheless. See simrun.m for the
%   legacy classdef version

    %cd('~/Documents/MATLAB')

    %% constants for general use
    r.k = 1.381e-23;
    r.mOH = 2.82328e-26; % Accounts for Oxygen binding energy
    r.uOH = 9.27401e-24 * 1.4;
    r.h = 6.62607e-34;
    r.hb = r.h/(2*pi);
    
    %% initialization of a run
    % Like the fortran sim, these can be put in braces to indicate several
    % runs over different parameter options.
    
    % variables for the initial distribution
    r.dname = 'delayswitch_speed_panel';
    r.num = 4e5;
    r.tempxy = 200e-3; %{100e-3 200e-3 400e-3 800e-3 1.6 3 6 12};
    r.spreadxy = 2e-3;
    r.tempz = 200e-3;
    r.spreadz = 5e-3;
    r.initvz = 820;
    r.dist = 'gaussian';
        
    % decelerator configuration variables
    r.vdd = 2e-3;
    
    % Choose from electrodering, uniformmagnet, normal, magneticpin,
    % varygap2pX, where X is from 0 to 5, 
    % ppmm_2mm, pmpm_2mm, pmpm_2mm_no-sym
    r.decels.a =  'longdecel'; %{'pmpm_2mm_no-sym','ppmm_2mm'};
    r.decels.b = 'singlerod';
    
    r.reloadfields = false;
    
    % Make sure are true except for guiding fields
    r.fieldsymmetryXY = true;
    r.fieldsymmetryZ = true;
    
    % decelerator timing variables
    p = 55;
    n = 333;
    r.chargetype = repmat({repmat('ab',1,n), repmat('a',1,n)},1,6);
    r.stages = repmat({floor((1:(2*n-1))/2+1), 1:n},1,6);
    r.rot180 = repmat({mod(floor((1:(2*n-1))/4),2), zeros(1,n)},1,6);
    r.endphases =   repmat({[p repmat([-p, p],1,n-1)], p*ones(1,n)},1,6);
    r.finalvz = {800 800 500 500 200 200 100 100 50 50 37 37};

    % simulation timing variables
    r.smallt = 1e-7;
    
    % laser beam variables
    r.lasertype = 'disk';
    r.LBD = 2.5e-3;
    
    % random number seed
    r.seed = 21112;
    
    %Unpack r into a struct of runs
    rs = unpacker(r,'linear');
    
    %% Here we just loop through the struct of runs, and run each one.
    for i=1:length(rs)
        rng(rs(i).seed) %seed the random number generator
        fprintf('run:%3d/%d\n ',i,length(rs))
        rr = init(rs(i));
        rsf(i) = run(rr);
        rsf(i).f = 0; %clear the fields, massive data sink.
    end
    % Save the struct of runs, just in case the fitsim2data or results
    % functions have petty errors. Wouldn't want to lost everything.
    t = datestr(now,'mmmm-dd-yyyy_HH-MM-SS');
    if ~exist('autosaves','dir')
        mkdir('autosaves')
    end
    save(['autosaves/rundecelstructs_' t '_' r.dname '.mat'],'rsf')
    system(['cp simdecel.m ./autosaves/simdecel_' t '_' r.dname '.m']);
   
    %disp(rsf(1).vels(end))
    %resultsdecel(rsf)
end

function r = init(r)
    r = initdecel(r);
    r = initmols(r);
    r = initvars(r);
end

function r = initdecel(r)
    % Load the mat file, or generate it from a COMSOL .dat file if it
    % doesn't exist yet.
    labels = fields(r.decels);
    for i=1:length(labels)
        d = r.decels.(labels{i});
        if exist(['Decels/' d '.mat'],'file') && ~r.reloadfields
            r.f.(labels{i}) = load(['Decels/' d '.mat']);
        elseif exist(['Decels/' d '.dat'],'file')
            r = processfields(r,d);
            r.f.(labels{i}) = load(['Decels/' d '.mat']);
        else
            error(['File ''Decels/' d '.dat'' not found']);
        end
    end
    
    % Choose the phase angle as a function of vfinal, vinitial, and stage
    % number.
    if r.finalvz
        energy = .5*r.mOH*(r.initvz^2 - r.finalvz^2);
        % changed the bounds for acceleration
        c = labels{1};
        r.phase = fminbnd(@(phi) (r.f.(c).renergy(phi)*max(r.stages) + ...
            r.f.(c).aenergy(-phi) - energy)^2,-90,90); 
        fprintf('Phase Angle: %2.3f\n',r.phase);
        curphase = mode(abs(r.endphases));
        r.endphases(r.endphases==curphase) = r.phase;
        r.endphases(r.endphases==-curphase) = -r.phase;
        r.endphases(r.endphases==curphase-180) = r.phase-180;
    end

end
%% Initialize Molecules
function r = initmols(r)
    spreads = repmat([r.spreadxy r.spreadxy r.spreadz],r.num,1);
    temps = repmat([r.tempxy r.tempxy r.tempz],r.num,1);
    if strcmpi(r.dist,'gaussian') || strcmpi(r.dist,'normal')
        r.vel = randn(r.num,3)*sqrt(r.k/r.mOH);
        r.vel = r.vel.*sqrt(temps);
        r.pos = randn(r.num,3)/sqrt(8*log(2));
        r.pos = r.pos.*spreads;
    elseif strcmpi(r.dist,'homogeneous') || strcmpi(r.dist,'flat')
        r.vel = (rand(r.num,3)-.5)*sqrt(r.k/r.mOH);
        r.vel = r.vel.*sqrt(temps);
        r.pos = (rand(r.num,3)-.5).*spreads;
    else
        error(['Distribution type ''' lower(r.dist) ''' not recognized.'...
            ' Capitalization not important.'])
    end

    %First molecule always 'perfect' for timing
    r.vel(1,:) = [0 0 0];
    r.pos(1,:) = [0 0 0];
    
    %Shift molecules in z, vz.
    r.vel(:,3) = r.vel(:,3) + r.initvz;
    r.pos(:,3) = r.pos(:,3) - r.vdd;
end

function r = initvars(r)
    % Stage Number
    r.numstage = 1;
    r.charge = r.chargetype(1);
    r.rot = 0;

    % Store the molecule number each decel stage.
    r.molnum = zeros(1,max(r.stages));
        
    % Timing
    r.time = 0;
end


%% This function does a run.
% It just times itself and calls the step function.
function r = run(r)
    r = tofirststage(r);
    while r.numstage <= max(r.stages)
        r = stage(r);
        fprintf('step:%3d/%d,\t%d\n',r.numstage,max(r.stages),r.molnum(r.numstage))
        r.numstage = r.numstage + 1;
    end
end
    
%% Propagate Molecules to first stage
% Eventually this could be modified to include hexapole focusing or maybe
% magnetic quadrupole focusing.
function r = tofirststage(r)
    %time = ((-90+r.phases(1))*r.f.zstagel/360-r.pos(1,3))/r.vel(1,3);
    % I used to artificially change the decelerator start time according to
    % the phase angle to make sure that even the first stage was equal
    % length with all others. This made phase angle calculations for a
    % given final velocity target easier. For collision hunting this is not
    % good. Instead:
    time = -r.pos(1,3)/r.vel(1,3);
    r.pos = r.pos + r.vel*time;
end

%% The stage function.
% Propagates one decel stage. Removes lost molecules, checks number, etc. A
% stage now includes variable numbers of sub-stages.
function r = stage(r)
    indices = find(r.stages==r.numstage);
    for i=1:length(indices)
        ind = indices(i);
        c = r.chargetype(ind);
        r.charge = c;
        r.rot = r.rot180(ind);
        while r.f.(c).phase(r.pos(1,3),r.numstage) < r.endphases(ind)
            r = smallstep(r,r.smallt);
            if r.vel(1,3) <= 0
                error('synchronous molecule reflected')
            end
        end
        undershoot = 1;
        while abs(undershoot) > 1e-11 && r.vel(1,3) > 0
            % get the undershoot in terms of phase
            undershoot = r.endphases(ind) - ...
                r.f.(c).phase(r.pos(1,3),r.numstage);

            % translate to time ignoring acceleration
            undershoot = undershoot *r.f.(c).zstagel/360/r.vel(1,3);

            % step the molecules according to this time
            r = smallstep(r,undershoot);
        end
        
    end
    
    r.pos(abs(r.pos(:,3)-r.pos(1,3))>7e-3,:)=nan;
    r.pos(r.vel(:,3)<0,:)=nan;
    r.lost    = isnan(sum(r.pos,2));
    r.pos     = r.pos(~r.lost,:);
    r.vel     = r.vel(~r.lost,:);
    r.numleft = size(r.pos,1);
    
    r.vels(r.numstage) = r.vel(1,3);
    r.times(r.numstage) = r.time;
    r.molnum(r.numstage) = r.numleft;
    
end

%% Small simulation step.
% Updates velocity based on acceleration and position based on velocity,
% offset by half of a timestep.
function r = smallstep(r,t)
    r.pos = r.pos + r.vel*t/2;
    r.vel = r.vel + acc(r)*t;
    r.pos = r.pos + r.vel*t/2;
            
    %update time.
    r.time = r.time + t;
end

%gets acceleration
function a = acc(r)
    % first rotate into the right frame:
    rad = pi/2*(mod(r.numstage,2) + 2*r.rot);
    c = cos(rad); s = sin(rad);
    pos = [ r.pos(:,1)*c - r.pos(:,2)*s , ...
            r.pos(:,1)*s + r.pos(:,2)*c , r.pos(:,3)];
    
    %just look up the force from the tables of dvdr (v as in potential
    %energy capital V.)
    ax = r.f.(r.charge).dvdx(pos,r.numstage);
    ay = r.f.(r.charge).dvdy(pos,r.numstage);
    az = r.f.(r.charge).dvdz(pos,r.numstage);
    a = [ax*c + ay*s , -ax*s + ay*c , az]/r.mOH;
end
    
function rs = unpacker(r,type)
    switch(type)
        case 'linear'
            fs = fields(r);
            l = 1;
            for i=1:length(fs)
                if iscell(r.(fs{i}))
                    l = length(r.(fs{i}));
                    break;
                end
            end
            
            rs = repmat(r,1,l);
            
            for i=1:length(fs)
                for j=1:l
                    tf = r.(fs{i});
                    if iscell(tf)
                        rs(j).(fs{i}) = tf{j};
                    else
                        rs(j).(fs{i}) = tf;
                    end
                end
            end
        case 'product'
            fs = fields(r);
            for i=1:length(fs)
                if ~iscell(r.(fs{i}))
                    r.(fs{i}) = {r.(fs{i})};
                    if i==1
                        continue
                    end
                end
                these = r.(fs{i});
                len = length(these);
                height = length(r.(fs{i-1}));
                for j=1:i-1
                    r.(fs{j}) = repmat(r.(fs{j}),1,len);
                end
                r.(fs{i}) = reshape(repmat(these,height,1),1,len*height);
            end
            structargs = [fs struct2cell(r)];
            structargs = structargs';
            rs = struct(structargs{:});
    end
end

%% Load in from COMSOL
% Takes a COMSOL .dat file containing fields in a deceleration stage
%
% Expectations:
% * z is the decelerator axis, ranging from phase angle of $-90^\circ$ to $+90^\circ$.
% * x and y are symmetric.
% * First data column is a geometry mask, then E-field, then B and angle if relevant.
% * The Geometry mask is 1 where obstacles exist.
% * E-field in V/m, B-field in Tesla, angle in radians, xyz in mm.
% * It is assumed that odd numbered stages have field symmetry with respect
% to even numbered stages, but rotated.
% * Data points are given on a rectangular, uniform grid, although the
% spacing of the three dimensions need not be identical.
function r = processfields(r,decel)
    
    % Announcement
    fprintf('%s\n','Processing Decelerator Fields from COMSOL...');

    % COMSOL files usually have 9 header lines.
    data = importdata(['Decels/' decel '.dat'],' ',9);
    
    % data is a struct with data and text header. We ignore the header and
    % take out the data.
    all = data.data;
    
    %% Convert data to 3D format.
    % After COMSOL export, the data are in list form. Each data row gives
    % the mask, efield, bfield, and angle evaluated at one point given by
    % x,y,z in the first three columns. 
    %
    % To get the data into 3D matrices, we first infer the 3D span of the
    % datapoints by checking for unique values in x,y,z. Then we get the
    % linear index of the 3D matrix corresponding to each row in the COMSOL
    % list. Finally we can fill the 3D matrix just be indexing into the 3D
    % matrix with this linear index.
    %
    
    % read x,y,z,m,e,b,t from each data column. 
    x = all(:,1)*1e-3; %convert to meters.
    y = all(:,2)*1e-3;
    z = all(:,3)*1e-3;
    m = all(:,4);
    e = all(:,5);   %units of V/m
    
    % not all decels will have b-field information.
    if size(all,2)>5
        b = all(:,6);   %units of T
        t = all(:,7);   %radians
        t(isnan(t))=0;
        t(imag(t)~=0)=0;
    else
        b = zeros(size(e));
        t = b;
    end
    
    % find unique x,y,z coordinates
    xs = sort(uniquetol(x,1e-6,'DataScale',1));
    ys = sort(uniquetol(y,1e-6,'DataScale',1));
    zs = sort(uniquetol(z,1e-6,'DataScale',1));
    
    % get the datapoint spacing, used for taking derivatives later.
    xsp = uniquetol(diff(xs));
    ysp = uniquetol(diff(ys));
    zsp = uniquetol(diff(zs));
    
    % z is assumed to run from $-90^\circ$ to $+90^\circ$, but no
    % assumptions are made about its actual coordinate range in COMSOL.
    % Thus it is shifted so that $+90^\circ$ phase is at zero.
    z = z - zs(end);
    zs = zs - zs(end);
    zstagel = -2*zs(1);
    
    % create lookup functions that tell you 'n' for a given x value such
    % that x is the nth x value.
    x2i = @(xx) arrayfun(@(x) find(x==xs),xx);
    y2i = @(yy) arrayfun(@(y) find(y==ys),yy);
    z2i = @(zz) arrayfun(@(z) find(z==zs),zz);
    
    % check for datapoint uniformity, and x,y starting from zero
    if length(xsp)+length(ysp)+length(zsp) > 3
        error('Non-uniform Datapoint Spacing');
    elseif abs(xs(1)) > 1e-6 && r.fieldsymmetryXY
        error('X data doesn''t begin at zero')
    elseif abs(ys(1)) > 1e-6 && r.fieldsymmetryXY
        error('Y data doesn''t begin at zero')
    end
    
    % the size of the 3D data matrices to be filled
    fullsize = [length(xs) length(ys) length(zs)];
    
    % for each x,y,z value in the COMSOL loaded x,y,z columns, check which
    % linear index this corresponds to in a 3D matrix of size fullsize.
    locs = sub2ind(fullsize,x2i(x),y2i(y),z2i(z));
    
    % for each single-letter variable in varname, create a new variable
    % given by a double-letter, which is a 3D matrix of size fullsize, and
    % fill it by indexing into it with the locs linear index column.
    for varname={'b','e','t','m','x','y','z'}
        eval([varname{1} varname{1} '=zeros(fullsize);']);
        eval([varname{1} varname{1} '(locs) = ' varname{1} ';']);
    end
    
    %% Calculate Stark-Zeemen Potential Energy
    
    % create an anonymous function that can give OH doubly stretched state
    % potential energy as a function of bfield, efield, and ebangle (t).
    last = @(x) x(end);
    energysingle = @(b,e,t) last(sort(eig(OH_Ham_Simple_SI(b,e,t))));
    energy = @(b,e,t)  arrayfun(energysingle,b,e,t);
    
    % get the potential energy
    vv = energy(bb,ee,tt);
    
    % smooth it for better derivatives. I tuned the standard deviation, 3,
    % until surfaces (run "figure;surf(squeeze(dvdzm(1,:,:)))" for example)
    % don't show waviness or chopiness.
    % vv = smooth3(vv,'gaussian',7,3);
    
    %% Take derivatives to get force fields
    
    % get symmetric derivative kernels so we can differentiate the
    % potential matrix via convolution
    xd = reshape([-.5 0 .5],3,1,1);
    yd = reshape(xd,1,3,1);
    zd = reshape(xd,1,1,3);
    
    % perform the derivatives. Convolution style differentiation requires
    % scaling by the matrix point spacing.
    dvdxu = convn(vv,xd,'same')/xsp;
    dvdyu = convn(vv,yd,'same')/ysp;
    dvdzu = convn(vv,zd,'same')/zsp;
    
    % zero the forces outside of the geometry mask. This ensures that
    % molecules aren't accidentally reflected off of pathologically large
    % field spikes that can occur near conductor or magnet surfaces in
    % COMSOL. Instead, molecules that hit geometry will continue through
    % until they fall out of the force field and are removed by the
    % molecule stepper which looks for this.
    mmx = convn(mm,abs(xd),'same') > 0;
    mmy = convn(mm,abs(yd),'same') > 0;
    mmz = convn(mm,abs(zd),'same') > 0;
    dvdxu(mmx)=0;
    dvdyu(mmy)=0;
    dvdzu(mmz)=0;
    if r.fieldsymmetryXY
        dvdxu(1,:,:)=0; 
        dvdyu(:,1,:)=0; 
        dvdxu(end,:,:)=dvdxu(end-1,:,:);
        dvdyu(:,end,:)=dvdyu(:,end-1,:);
    else
        dvdxu([1 end],:,:)=dvdxu([2,end-1],:,:);
        dvdyu(:,[1 end],:)=dvdyu(:,[2,end-1],:);
    end
    dvdzu(:,:,[1 end])=dvdzu(:,:,[2, end-1]);    
    dvdxu(mmx)=2e-19;
    dvdyu(mmy)=2e-19;
    dvdzu(mmz)=2e-19;
    dvdxm = shiftableBF3D(dvdxu,2,2e-20,1e-4,2e-19);
    dvdym = shiftableBF3D(dvdyu,2,2e-20,1e-4,2e-19);
    dvdzm = shiftableBF3D(dvdzu,2,2e-20,1e-4,2e-19);
    dvdxm(mmx)=0;
    dvdym(mmy)=0;
    dvdzm(mmz)=0;
        
    % zero the x,y force along their respective lines of symmetry. This
    % should already be the case but convolution based derivatives can
    % behave strangely near borders due to zero padding assumptions.
    if r.fieldsymmetryXY
        dvdxm(1,:,:)=0;
        dvdym(:,1,:)=0;
    end
    if r.fieldsymmetryZ
        dvdzm(:,:,[1 end])=0;
    end
    
    % finally, kill the fields inside obstacles:
    dvdxm(mmx)=nan;
    dvdym(mmy)=nan;
    dvdzm(mmz)=nan;
    dvdxm(~~mm)=nan;
    dvdym(~~mm)=nan;
    dvdzm(~~mm)=nan;
    
    %% Re-create Potential on Axis
    % After the smoothing, the forcefields are a bit different than the
    % original potential in terms of their integrated potential energy.
    % This creates timing discrepancies if not addressed.
    dvdzax = squeeze(dvdzm(1,1,:));
    zax = squeeze(zz(1,1,:));
    vax = zax;
    vax(1) = 0;
    for i=2:length(vax)
        vax(i) = -trapz(zax(1:i),dvdzax(1:i));
    end
    
    %% Create interpolants
    % These convenient datatypes can be evaluated directly as functions and
    % carry with them all of the gridded data used to instantiate them.
    dvdxg  = griddedInterpolant(xx,yy,zz,dvdxm,'linear','none');
    dvdyg  = griddedInterpolant(xx,yy,zz,dvdym,'linear','none');
    dvdzg  = griddedInterpolant(xx,yy,zz,dvdzm,'linear','none');
    vfg = griddedInterpolant(xx,yy,zz,vv);
    vaxg = griddedInterpolant(zax,vax,'linear','none');
    % bf = griddedInterpolant(xx,yy,zz,bb);
    % ef = griddedInterpolant(xx,yy,zz,ee);
    % mf = griddedInterpolant(xx,yy,zz,mm);
    % tf = griddedInterpolant(xx,yy,zz,tt);
    
    

    
    
    %% Create helpful lookup functions
    % Anonymous functions created here will be saved along with the
    % relevant workspace when defined. This is an automatic matlab feature.
    % This allows the intricacies of stage parity and reflection symmetry
    % to be encapsulated here and not dealt with in the functions that use
    % these lookup tables.
    
    
    % This phase lookup function gives the phase angle as a function of the
    % z coordinate and the parity of the stage number n. It gives the phase
    % in the range -270 to 90, which is convenient since the decelerator
    % will be run with at a phase angle close to +90 degrees, and thus
    % during a single stage most well-decelerated molecules won't wrap
    % their phase angles as returned by this lookup.
    phase = @(z,n) mod(z/zstagel+n/2,1)*360 - 270;
    
    % wrap returns a z-coordinate within the force lookup table given a
    % general z-coordinate and the stage parity. It achieves this in two
    % steps. First z is converted to a coordinate between -zstagel and 0
    % which is also the range -270 to +90 in phase angle. Then, the mirror
    % symmetry between the forces in -270--90 and -90-+90 is exploited. The
    % side lookup indicates whether this symmetry is exploited so the
    % z-forces can be inverted, since molecules in the -270 to -90 range
    % are accelerated, not decelerated.
    wrapc = @(z,n) (phase(z,n)-90)/360 * zstagel;
    wrap = @(z,n) abs(wrapc(z,n)+zstagel/2)-zstagel/2;
    side = @(z,n) (wrapc(z,n) > -zstagel/2)*2 - 1;
    if ~r.fieldsymmetryZ
        wrap = wrapc;
        side = @(z,n) 1;
    end
    
    % This returns the energy removed per stage as a function of phase
    % angle. Its inverse enables quickly choosing the phase angle given a
    % final velocity. 
    renergy = @(phi) vaxg((phi-90)/360 * zstagel) - ...
        vaxg((-phi-90)/360 * zstagel);
    
    % This is the potential energy at a given phase angle, measured
    % relative to the potential energy at phi=-90 degrees. One could
    % subtract this from itself reversed to get renergy above.
    aenergy = @(phi) vaxg((phi-90)/360 * zstagel) - vaxg(-zstagel/2);

    % These functions reference the gridded interpolants, but with the
    % coordinates appropriately wrapped.
    if r.fieldsymmetryXY
        dvdxa = @(x,y,z,n) dvdxg(abs(x),abs(y),wrap(z,n)).*sign(x);
        dvdya = @(x,y,z,n) dvdyg(abs(x),abs(y),wrap(z,n)).*sign(y);
        dvdza = @(x,y,z,n) dvdzg(abs(x),abs(y),wrap(z,n)).*side(z,n);
        vfa = @(x,y,z,n) vfg(abs(x),abs(y),wrap(z,n));
    else
        dvdxa = @(x,y,z,n) dvdxg(x,y,wrap(z,n));
        dvdya = @(x,y,z,n) dvdyg(x,y,wrap(z,n));
        dvdza = @(x,y,z,n) dvdzg(x,y,wrap(z,n)).*side(z,n);
        vfa = @(x,y,z,n) vfg(x,y,wrap(z,n));
    end
    dvdx = @(xyz,n) dvdxa(xyz(:,1),xyz(:,2),xyz(:,3),n);
    dvdy = @(xyz,n) dvdya(xyz(:,1),xyz(:,2),xyz(:,3),n);
    dvdz = @(xyz,n) dvdza(xyz(:,1),xyz(:,2),xyz(:,3),n);
    vf = @(xyz,n) vfa(xyz(:,1),xyz(:,2),xyz(:,3),n);
    

    % Save in a file for loading and propagating during decel simulation.
    save(['Decels/' decel '.mat'],'dvdx','dvdy','dvdz',...
        'vf','phase','zstagel','renergy','aenergy');

    %% Produce Output Figure for Debugging
    % There are many potential errors that could be made in the COMSOL
    % output, data input processing. 
    cc = 4e-20;
    
    figure('position',[50,50,1100,1100])
    subplot(2,2,1)
    surf(cap(squeeze(dvdzm(:,1,:)),cc));
    title(['Decelerator dvdz, X-Z plane, ' decel '.dat']);

    % You might think the labels are backwards, but they're not. surf uses
    % the second index (the column of the matrix) as the x-axis and the
    % first index (the row) as the y-axis. It makes sense if you think of
    % matrices as oriented with x going left-right and y going up-down, but
    % its crazy when working in 3D.
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    subplot(2,2,2)
    surf(cap(squeeze(dvdzm(1,:,:)),cc));
    title(['Decelerator dvdz, Y-Z plane, ' decel '.dat']); 
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    subplot(2,2,3)
    surf(cap(squeeze(dvdxm(:,1,:)),cc));
    title(['Decelerator dvdx, X-Z plane, ' decel '.dat']);
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    subplot(2,2,4)
    surf(cap(squeeze(dvdym(1,:,:)),cc));
    title(['Decelerator dvdy, Y-Z plane, ' decel '.dat']);
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    %{
    figure('position',[100,200,400,400])
    plot(abs(zs),vf([zeros(length(zs),2) zs(:)],1))
    title(['Decel Potential along Z-axis, ' r.trapname '.dat'])
    xlabel('Z axis')
    ylabel('Potential Energy (J)')
    %}
end

function x = cap(x,varargin)
    if length(varargin)==1
        c = varargin{1};
        x(x>c) = c;
        x(x<-c) = -c;
    else
        [cl, ch] = varargin{:};
        x(x>ch) = ch;
        x(x<cl) = cl;
    end
end







