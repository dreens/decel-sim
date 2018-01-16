function rsf = simdecel()
%MATLAB Simulation of OH experiment
%   I had a classdef version of this before but it was brutally slow.
%   turns out the slowness wasn't due to it being a class, but I'm going to
%   stick with this struct implementation nonetheless. See simrun.m for the
%   legacy classdef version

    %cd('~/Documents/MATLAB')

    %% constants for general use
    r.k = 1.381e-23;
    r.mOH = 17*1.6726e-27;
    r.uOH = 9.27401e-24 * 1.4;
    r.h = 6.62607e-34;
    r.hb = r.h/(2*pi);
    
    %% initialization of a run
    % Like the fortran sim, these can be put in braces to indicate several
    % runs over different parameter options.
    
    % variables for the initial distribution
    r.dname = 'Detection_Details';
    r.num = 1e4;
    r.tempxy = 200e-3;
    r.spreadxy = 3e-3;
    r.tempz = 200e-3;
    r.spreadz = 4e-3;
    r.initvz = 810;
    r.dist = 'homogeneous';
    
    % decelerator configuration variables
    r.stages = 333;%{100,125,150,175,200,225,250,275,300};      
    r.vdd = 10e-3;
    
    % Choose from electrodering, uniformmagnet, normal, magneticpin,
    % varygap2pX, where X is from 0 to 5, 
    r.decel = 'normal';%{'varygap2p0','varygap2p1','varygap2p2','varygap2p3','varygap2p4','varygap2p5'};
    
    % decelerator timing variables
    r.deceltiming = 'finalvz';
    r.finalvz = 50;%{10 20 30 40 50 75 100 150};
    r.phase = 75;%{0 10 20 30 30.39};
    
    % simulation timing variables
    r.smallt = 1e-7;
    
    % laser beam variables
    r.lasertype = 'disk';
    r.LBD = 2.5e-3;
    
    % random number seed
    r.seed = 21112;
    
    %Unpack r into a struct of runs
    rs = unpacker(r,'product');
    
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
    resultsdecel(rsf)
end

function r = init(r)
    r = initdecel(r);
    r = initmols(r);
    r = initvars(r);
end

function r = initdecel(r)
    % Load the mat file, or generate it from a COMSOL .dat file if it
    % doesn't exist yet.
    if exist(['Decels/' r.decel '.mat'],'file');
        r.f = load(['Decels/' r.decel '.mat']);
    elseif exist(['Decels/' r.decel '.dat'],'file');
        r = processfields(r);
        r.f = load(['Decels/' r.decel '.mat']);
    else
        error(['File ''Decels/' r.decel '.dat'' not found']);
    end
    
    % Choose the phase angle as a function of vfinal, vinitial, and stage
    % number.
    if strcmp(r.deceltiming,'finalvz')
        energyper = .5*r.mOH*(r.initvz^2 - r.finalvz^2)/r.stages;
        r.phase = fminbnd(@(phi) (r.f.renergy(phi)-energyper)^2,0,90);
        fprintf('Phase Angle: %2.3f\n',r.phase);    
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
        error(['Distribution type ''' lower(r.dist) ''' not recognized. Capitalization not important.'])
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

    % Store the molecule number each decel stage.
    r.molnum = zeros(1,r.stages);
    
    % Store the full phase space each stage.
    %r.phasestage = zeros(r.num,6,r.stages);

    % other variables used by different runs
    %r.lost = false(r.num,1);
    
    % Timing
    r.time = 0;
end


%% This function does a run.
% It just times itself and calls the step function.
function r = run(r)
    r = tofirststage(r);
    while r.numstage <= r.stages
        r = stage(r);
        fprintf('step:%3d/%d,\t%d\n',r.numstage,r.stages,r.molnum(r.numstage))
        r.numstage = r.numstage + 1;
        %r.pos(:,1:2) = r.pos(:,[2 1]); r.pos(:,1) = -r.pos(:,1);
        %r.vel(:,1:2) = r.vel(:,[2 1]); r.vel(:,1) = -r.vel(:,1);
    end
end

%% Propagate Molecules to first stage
% Eventually this could be modified to include hexapole focusing or maybe
% magnetic quadrupole focusing.
function r = tofirststage(r)
    time = ((-90+r.phase)*r.f.zstagel/360-r.pos(1,3))/r.vel(1,3);
    r.pos = r.pos + r.vel*time;
end

%% The stage function.
% Propagates one decel stage. Removes lost molecules, checks number, etc.
function r = stage(r)
    % Step until the synchronous molecule is past the required phase angle.
    while r.f.phase(r.pos(1,3),r.numstage) < r.phase  && r.vel(1,3) > 0
        r = smallstep(r,r.smallt);
    end
    
    % Iterate forward/backward in time until synchronous molecule is right
    % on top of the correct phase angle.
    undershoot = 1;
    while abs(undershoot) > 1e-9 && r.vel(1,3) > 0
        % get the undershoot in terms of phase
        undershoot = r.phase - r.f.phase(r.pos(1,3),r.numstage);
        
        % translate to time ignoring acceleration
        undershoot = undershoot *r.f.zstagel/360/r.vel(1,3);
        
        % step the molecules according to this time
        r = smallstep(r,undershoot);
    end
    
    r.pos(abs(r.pos(:,3)-r.pos(1,3))>15e-3,:)=nan;
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
    
    %just look up the force from the tables of dvdr (v as in potential
    %energy capital V.)
    if mod(r.numstage,2)
        ax = r.f.dvdx(r.pos,r.numstage);
        ay = r.f.dvdy(r.pos,r.numstage);
        az = r.f.dvdz(r.pos,r.numstage);
        a = [ax ay az]/r.mOH;
    else
        ay = r.f.dvdx(r.pos(:,[2 1 3]),r.numstage);
        ax = r.f.dvdy(r.pos(:,[2 1 3]),r.numstage);
        az = r.f.dvdz(r.pos(:,[2 1 3]),r.numstage);
        a = [ax ay az]/r.mOH;
    end
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
function r = processfields(r)
    
    % COMSOL files usually have 9 header lines.
    data = importdata(['Decels/' r.decel '.dat'],' ',9);
    
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
    x2i = @(x) int16(1+x/xsp);
    y2i = @(y) int16(1+y/ysp);
    z2i = @(z) int16(length(zs)+z/zsp);
    
    % check for datapoint uniformity, and x,y starting from zero
    if length(xsp)+length(ysp)+length(zsp) > 3
        error('Non-uniform Datapoint Spacing');
    elseif abs(xs(1)) > 1e-6
        error('X data doesn''t begin at zero')
    elseif abs(ys(1)) > 1e-6
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
    vv = smooth3(vv,'gaussian',7,3);
    
    %% Take derivatives to get force fields
    
    % get symmetric derivative kernels so we can differentiate the
    % potential matrix via convolution
    xd = reshape([-.5 0 .5],3,1,1);
    yd = reshape(xd,1,3,1);
    zd = reshape(xd,1,1,3);
    
    % perform the derivatives. Convolution style differentiation requires
    % scaling by the matrix point spacing.
    dvdxm = convn(vv,xd,'same')/xsp;
    dvdym = convn(vv,yd,'same')/ysp;
    dvdzm = convn(vv,zd,'same')/zsp;
    
    % zero the forces outside of the geometry mask. This ensures that
    % molecules aren't accidentally reflected off of pathologically large
    % field spikes that can occur near conductor or magnet surfaces in
    % COMSOL. Instead, molecules that hit geometry will continue through
    % until they fall out of the force field and are removed by the
    % molecule stepper which looks for this.
    dvdxm(mm~=0)=0;
    dvdym(mm~=0)=0;
    dvdzm(mm~=0)=0;
    
    % zero the x,y force along their respective lines of symmetry. This
    % should already be the case but convolution based derivatives can
    % behave strangely near borders due to zero padding assumptions.
    dvdxm(1,:,:)=0; dvdxm(end,:,:)=0;
    dvdym(:,1,:)=0; dvdym(:,end,:)=0;
    dvdzm(:,:,1)=0; dvdzm(:,:,end)=0;
    
    %% Create interpolants
    % These convenient datatypes can be evaluated directly as functions and
    % carry with them all of the gridded data used to instantiate them.
    dvdxg  = griddedInterpolant(xx,yy,zz,dvdxm,'linear','none');
    dvdyg  = griddedInterpolant(xx,yy,zz,dvdym,'linear','none');
    dvdzg  = griddedInterpolant(xx,yy,zz,dvdzm,'linear','none');
    vf = griddedInterpolant(xx,yy,zz,vv);
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
    
    % This returns the energy removed per stage as a function of phase
    % angle. Its inverse enables quickly choosing the phase angle given a
    % final velocity. 
    renergy = @(phi) vf(0,0,(phi-90)/360 * zstagel) - ...
        vf(0,0,(-phi-90)/360 * zstagel);

    % These functions reference the gridded interpolants, but with the
    % coordinates appropriately wrapped.
    dvdxa = @(x,y,z,n) dvdxg(abs(x),abs(y),wrap(z,n)).*sign(x);
    dvdya = @(x,y,z,n) dvdyg(abs(x),abs(y),wrap(z,n)).*sign(y);
    dvdza = @(x,y,z,n) dvdzg(abs(x),abs(y),wrap(z,n)).*side(z,n);
    dvdx = @(xyz,n) dvdxa(xyz(:,1),xyz(:,2),xyz(:,3),n);
    dvdy = @(xyz,n) dvdya(xyz(:,1),xyz(:,2),xyz(:,3),n);
    dvdz = @(xyz,n) dvdza(xyz(:,1),xyz(:,2),xyz(:,3),n);
    
    % Save in a file for loading and propagating during decel simulation.
    save(['Decels/' r.decel '.mat'],'dvdx','dvdy','dvdz','phase','zstagel','renergy')
end




