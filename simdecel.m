function rsf = simdecel(varargin)
%MATLAB Simulation of OH experiment
    %% constants for general use
    ri.k = 1.381e-23;
    ri.mOH = 2.82328e-26; % Accounts for Oxygen binding energy
    ri.uOH = 9.27401e-24 * 1.4;
    ri.h = 6.62607e-34;
    ri.hb = ri.h/(2*pi);
    
    ri.verbose = true;
    
    %% initialization of a run
    % Like the fortran sim, these can be put in braces to indicate several
    % runs over different parameter options.
    
    % variables for the initial distribution
    ri.dname = 'Vary_Initial_Speed_Length';
    
    % initial number:
    ri.num = 1e4;
    
    % initial temperature, spatial distribution:
    ri.tempxy = 5; %{100e-3 200e-3 400e-3 800e-3 1.6 3 6 12};
    ri.spreadxy = 3e-3;
    ri.tempz = 5;
    ri.spreadz = 8e-3;
    ri.initvz = 900;
    ri.dist = 'flat'; % or gaussian, spherical, other options.
    ri.vdd = 1e-3; % valve decelerator distance
    
    % a few variables for messing with a turnon time.
    ri.turnon = 400e-9;
    ri.everyother = false;

    
    % set this to load fields from a previous run. Useful for separating
    % deceleration from trap loading for example.
    ri.continue = false;
    phaserange = 54:2:84; iii=1;
    phaserange = 0;
    
    % you can just name the file to load the previous run:
    ri.contname = cell(1,length(phaserange));
    
    % or specify stubs to be prepended to the filenames in contname:
    ri.contfillstub = 'VSF';% {'S1','SF','VSF'};
    for n=phaserange
        ri.contname{iii} = sprintf('W0%%s56p%d',n); iii = iii + 1;
    end
    ri.contname = 'none';
    
    %% Sequence Specification
    %
    % Load ri.decels with one field for each charge configuration. The
    % field needs to be a single letter. The value of the field is the file
    % name without '.dat', or occasionally an arraylist where the first
    % field is the filename and the second contains some processing
    % instructions.
    %
    % Choose from electrodering, uniformmagnet, normal, magneticpin,
    % varygap2pX, where X is from 0 to 5, 
    % ppmm_2mm, pmpm_2mm, pmpm_2mm_no-sym, singlerod, ppgg
    %ri.decels = struct('a','longdecel','b','longdecel');
    ri.decels = struct('a','longdecel','b',{'singlerod','noXY'});
    %ri.decels = struct('a','longdecel','b','ppgg');
    %ri.decels = struct('a','longdecel','b','ppmm_2mm');
    
    
    ri.reloadfields = false;
        
    % decelerator timing variables
    ri.phase = 50;
    ri.phi2off = 0;
    
    ri.initvz = 900;
    i = 1;
    for n = 300
        ri.chargetype{i} = repmat('ba',1,n);
        rt = [0 0 90 90 180 180 270 270];
        rt = repmat(rt,1,ceil(n/4));
        ri.rot{i} = rt(1:2*n);

        tran = [1 1 0 0];
        tran = repmat(tran,1,ceil(n/2));
        ri.trans{i} = tran(1:2*n);
    
    
        % the inf flag will get replaced with the ri.phase variable
        %ri.endphases{i} = repmat([125 235 305 55],1,n); % S=1
        ri.endphases{i} = [repmat([125 235 305 55],1,n)]; % SF
        %ri.endphases{i} = [repmat([145 234.3 325 54.3],1,n)]; % VSF
        %ri.endphases{i} = [repmat([150 229.35 330 49.35],1,n)]; % XSF
        ri.finalvz = 0; % leaving the 'inf' flag allows phase tuning for final vz.

        % each stage has its endpoint calculated by phase, velocity, or time.
        % specify this as 'p', 'v', or 't':
        ri.calctype{i} = [repmat('pp',1,n)];

        i = i + 1;
    end
    
    % simulation timing variables
    ri.smallt = 2e-7;
    ri.reflectEnd = true; % end if the synch molecule is reflected.
    
    % laser beam variables
    ri.lasertype = 'disk';
    ri.LBD = 2.5e-3;
    
    % random number seed
    ri.seed = 21112;
    
    % overwrite directly passed variables
    lv = length(varargin);
    assert(~mod(lv,2),'Arguments should be specified in pairs, but %d arguments passed.',lv);
    for i=1:lv/2
        fn = varargin{2*i-1};
        if ~isfield(ri,fn)
            warning('simdecelWarn:notfield','Field %s does not appear to be a valid simdecel specification variable',fn)
            warning('off','simdecelWarn:notfield')
        end
        ri.(varargin{2*i-1}) = varargin{2*i};
    end
    
    %Unpack r into a struct of runs
    rs = unpacker(ri,'linear');
    
    % Making r global doesn't save time, but its silly passing it back and
    % forth all the time.
    global r 
    
    %% Here we just loop through the struct of runs, and run each one.
    for i=1:length(rs)
        rng(rs(i).seed) %seed the random number generator
        if rs(i).verbose
            fprintf('run:%3d/%d\n ',i,length(rs))
        end
        r = rs(i);
        init();
        run();
        r.f = 0; %clear the fields, massive data sink.
        rsf(i) = r;
        if r.verbose
            fprintf('speed:%3.1f\ttime:%1.3f\n',r.vels(end),r.time*1e3)
        end
    end
    % Save the struct of runs, just in case the fitsim2data or results
    % functions have petty errors. Wouldn't want to lost everything.
    t = datestr(now,'mmmm-dd-yyyy_HH-MM-SS');
    if ~exist('autosaves','dir')
        mkdir('autosaves')
    end
    save(['autosaves/rundecelstructs_' t '_' r.dname '.mat'],'rsf')
    system(['cp simdecel.m ./autosaves/simdecel_' t '_' r.dname '.m']);
    
    % save in Partials folder for reloading to simulate load, trap.
    %save('Partials/endBetweenLast4.mat','rsf')
    
    %disp(rsf(1).vels(end))
    %rsf(end).time*1e3
    %resultsdecel(rsf)
end

function init()
    global r
    %% Initialize fields, 
    % Load the mat file, or generate it from a COMSOL .dat file if it
    % doesn't exist yet.

    % r.decels will be a struct with a field for each charge configuration
    labels = fields(r.decels);
    
    % for each charge configuration...
    for i=1:length(labels)
        % d is the file name + any other processing insturctions for that
        % configuration.
        d = {r.decels.(labels{i})};
        
        % deal with the no instruction case, get filename as dname.
        if iscell(d)
            dname = d{1};
        else
            dname = d;
            d = {d};
        end
        
        % send messages and pass off to processfields if not already done.
        if ~exist(['Fields/' dname '.mat'],'file') || r.reloadfields
            if exist(['Fields/' dname '.dat'],'file')
                processfields(d{:});
            else
                error(['File ''Fields/' dname '.dat'' not found']);
            end
        end
        
        % add a field to the r.f structure with the data for the
        % corresponding charge configuration.
        r.f.(labels{i}) = load(['Fields/' dname '.mat']);
    end
    
    % Insert phase wherever 'inf' flag is found:
    if r.verbose
        fprintf('%s%.2f\n','Phase: ',r.phase)
    end
    r.endphases(r.endphases==inf) = r.phase;
    r.endphases(r.endphases==-inf) = -r.phase;
    
    % Choose the phase angle as a function of vfinal, vinitial, and stage
    % number. Don't do this if r.finalvz is set to zero.
    if r.finalvz 
        energy = .5*r.mOH*(r.initvz^2 - r.finalvz^2);
        % changed the bounds for acceleration
        c = labels{1};
        d = labels{2};
        r.phase = fminbnd(@(phi) (r.f.(c).aenergy(phi)*(max(r.stages)-1) ...
            + r.f.(d).aenergy(-phi+r.phi2off)*(max(r.stages)-1) ...
            - r.f.(c).aenergy(-phi+r.phi2off)*(max(r.stages)-2) ...
            - r.f.(d).aenergy(-phi)*(max(r.stages)-1) ...
            - r.f.(c).aenergy(-90) - energy)^2,-90,90); 
        if r.verbose
            fprintf('Phase Angle: %2.3f\n',r.phase);
        end
        p1 = max(r.endphases);
        p2 = min(r.endphases);
        r.endphases(r.endphases==p1) = r.phase;
        r.endphases(r.endphases==p2) = -r.phase+r.phi2off;
    end
    

    %% Initialize Molecules
    % as long as we aren't loading molecules from a previous run:
    if ~r.continue
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
        elseif strcmpi(r.dist,'sphere') || strcmpi(r.dist,'spherical')
            ap = spheredist(r.num,6);
            ap = ap.*[spreads sqrt(temps*r.k/r.mOH)]/2;
            r.vel = ap(:,4:6);
            r.pos = ap(:,1:3);
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
    else
        name = sprintf(r.contname,r.contfillstub);
        rsl = load(['Partials/' name]);
        r.vel = rsl.r.vel;
        r.pos = rsl.r.pos;
    end

    %% Initialize other variables
    % Stage Number
    r.numstage = 1;
    r.numstages = length(r.chargetype);
    r.charge = r.chargetype(1);
    r.rot = r.rot * pi / 180;

    % Store variables each decel stage.
    r.molnum = zeros(1,r.numstages);
    %r.vels = zeros(1,r.numstages);
    %r.times = zeros(1,r.numstages);
        
    % Timing
    r.time = 0;
    
    % Field Translations
    r.istrans = false;
    
    % Track Geometry
    r.xx = zeros(30000,1);
    r.smallnum = 1;
end


%% This function does a run.
% It just times itself and calls the step function.
function run()
    global r
    
    % Go to the first stage:
    if ~r.continue
        time = -r.pos(1,3)/r.vel(1,3);
        r.pos = r.pos + r.vel*time;    
    end

    while r.numstage <= r.numstages
        if stage()
            fprintf('all gone.\n')
            break
        end
        if ~mod(r.numstage,10) || r.numstages - r.numstage < 10
            if r.verbose
                fprintf('step:%3d/%d,\tN:%d\tv:%.0f\n',r.numstage,r.numstages,r.molnum(r.numstage),r.vel(1,3))
            end
        end
        r.numstage = r.numstage + 1;
    end
end


%% Translator
function checktrans()
global r
    if xor(r.istrans,r.trans(r.numstage))
        c = r.chargetype(r.numstage);
        if r.istrans
            r.pos(:,3) = r.pos(:,3) - 5e-3;%r.f.(c).zstagel/2;
        else
            r.pos(:,3) = r.pos(:,3) + 5e-3;%r.f.(c).zstagel/2;
        end
        r.istrans = ~r.istrans;
    end
end

%% The stage function.
% Propagates one decel stage. Removes lost molecules, checks number, etc.
function gone = stage()
    global r

    c = r.chargetype(r.numstage);
    r.charge = c;

    % Handle 'translations' of the potentials by artifically futzing with
    % the z coordinates.
    %checktrans()
    
    % Cluge in "turn-on" time by reducing acceleration over a certain
    % timeframe.
    r.turnontime = 0;
    if r.everyother && mod(r.numstage,2)
        r.turnontime = 1; % essentially instant on
    end
        
    % Make sure the synch molecule is "in" the fields. If not assume
    % loading and redefine z-coordinates
    z = r.pos(1,3);
    if isnan(r.f.(c).dvdz(0,0,z))
        finalz = mod(z,10e-3)-10e-3;
        trans = z - finalz;
        r.pos(:,3) = r.pos(:,3) - trans;
        if isnan(r.f.(c).dvdz(0,0,z))
            % alignment of loading with decel depends on how the last pin
            % pair is used.
            r.pos(:,3) = r.pos(:,3) + 10e-3;
        end
    end

    
    warned = false;
    while ~done()
        smallstep(r.smallt);
        if r.vel(1,3) <= 0
            if r.reflectEnd
                error('synchronous molecule reflected')
            elseif ~warned
                fprintf('synchronous molecule reflected\n')
                warned = true;
            end
        elseif isnan(r.vel(1,3))
            if ~strcmp(r.calctype(r.numstage),'t') % don't end sim if a timed stage (e.g. trapping)
                r.pos(:,:) = nan;
                fprintf('synchronous molecule lost\n')
            end
        end
    end
    
    function b = done()
        switch r.calctype(r.numstage)
            case 'p'
                difference = mod(r.f.(c).phase(r.pos(1,3)),360) - ...
                    r.endphases(r.numstage);
                if difference < 0
                    b = false;
                elseif difference > 0
                    b = true;
                    if difference > 180
                        b = false;
                    end
                end
                b = b || isnan(r.vel(1,3));
            case 'v'
                b = r.vel(1,3) <= r.endphases(r.numstage);
                b = b || isnan(r.vel(1,3));
            case 't'
                b = r.time >= r.endphases(r.numstage);
        end
    end
    
    if strcmp(r.calctype(r.numstage),'p')
        undershoot = 1;
        while abs(undershoot) > 1e-11 && r.vel(1,3) > 0
            % get the undershoot in terms of phase
            undershoot = r.endphases(r.numstage) - ...
                mod(r.f.(c).phase(r.pos(1,3)),360);
            
            % if significant undershoot, throw error
            if undershoot > 5
                message = sprintf('%s%s%d%s%d.',...
                    'Too much overshoot, stage definition issue.',...
                    ' Synch molecule at phase angle ',...
                    mod(r.f.(c).phase(r.pos(1,3)),360),...
                    ', but target phase was ',r.endphases(r.numstage));
                error(message)
            end

            % translate to time ignoring acceleration
            undershoot = undershoot *r.f.(c).zstagel/360/r.vel(1,3);

            % step the molecules according to this time
            smallstep(undershoot);
        end
    end
    
    %r.pos(abs(r.pos(:,3)-r.pos(1,3))>10e-3,:)=nan;
    %r.pos(r.vel(:,3)<0,:)=nan;
    r.lost    = isnan(sum(r.pos,2));
    r.pos     = r.pos(~r.lost,:);
    r.vel     = r.vel(~r.lost,:);
    r.numleft = size(r.pos,1);
    
    r.times(r.numstage) = r.time;
    r.molnum(r.numstage) = r.numleft;
    gone = ~r.numleft;
    if ~gone
        r.vels(r.numstage) = r.vel(1,3);
    else
        r.vels(r.numstage) = nan;
    end
    
end

%% Small simulation step.
% Updates velocity based on acceleration and position based on velocity,
% offset by half of a timestep.
function smallstep(t)
    global r
    r.pos = r.pos + r.vel*t/2;
    r.vel = r.vel + acc()*t;
    r.pos = r.pos + r.vel*t/2;
            
    %update time.
    r.time = r.time + t;
    r.turnontime = r.turnontime + t;
end

%gets acceleration
function a = acc()
    global r

    % mess around to check geometry
    %r.pos(2,:) = r.pos(1,:);
    %r.pos(3,:) = r.pos(1,:);
    %r.pos(2,1) = 1.9e-3;
    %r.pos(3,2) = 1.9e-3;

    
    % first rotate, translate into the right frame:
    rad = r.rot(r.numstage);
    c = cos(rad); s = sin(rad);
    x = r.pos(:,1)*c - r.pos(:,2)*s;
    y = r.pos(:,1)*s + r.pos(:,2)*c;
    z = r.pos(:,3) + r.trans(r.numstage)*(5e-3); %translate
  
    
    %just look up the force from the tables of dvdr (v as in potential
    %energy capital V.)
    ax = r.f.(r.charge).dvdx(x,y,z);
    ay = r.f.(r.charge).dvdy(x,y,z);
    az = r.f.(r.charge).dvdz(x,y,z);
    a = [ax*c + ay*s , -ax*s + ay*c , az]/r.mOH;
    a = a * (1 - exp(-r.turnontime/r.turnon));
    
    %r.xx(r.smallnum) = any(isnan(a(2,:)));
    %r.yy(r.smallnum) = any(isnan(a(3,:)));
    %r.smallnum = r.smallnum + 1;   
    %a(2:3,:) = 0;
    
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
function processfields(varargin)
    
    % Get Instructions from varargin
    fileN = varargin{1};
    if nargin>1
        directions = varargin{2};
    else
        directions = 'none';
    end

    % Announcement
    fprintf('Processing %s Fields from COMSOL...\n',fileN);
    
    % COMSOL files usually have 9 header lines.
    data = importdata(['Fields/' fileN '.dat'],' ',9);
    
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
    
    % depracated mask convention error.
    if strcmp(fileN,'tricycleload') || strcmp(fileN,'tricycletrap')
        m = ~m;
    end
    
    % not all decels will have b-field information.
    if size(all,2)>5
        b = all(:,6);   %units of T
        t = all(:,7);   %radians
        t(isnan(t))=0;  
        b(isnan(b))=0;  % b, e, m shouldn't have nans unless there are some mesh errors.
        e(isnan(e))=0;
        m(isnan(m))=0;
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
    
    % shift z coordinates depending on type of fields
    switch(directions)
        case 'load'
            z = z - zs(1) - 5e-3; % start fields at -5e-3 like last stage.
            zs = sort(uniquetol(z,1e-6,'DataScale',1));
            zstagel = zs(end);
        case 'trap'
            % here the coordinates need to be specified. h gives the
            % distance from the center of the last pin pair in m to the
            % beginning of the trapping fields.
            if nargin > 2
                h = varargin{3};
            else
                h = 0; 
            end
            z = z - zs(1) + h;
            zs = sort(uniquetol(z,1e-6,'DataScale',1));
            zstagel = zs(end);
        otherwise % i.e. deceleration
            % z is assumed to run from $-90^\circ$ to $+90^\circ$, but no
            % assumptions are made about its actual coordinate range in COMSOL.
            % Thus it is shifted so that $+90^\circ$ phase is at zero.
            z = z - zs(end);
            zs = zs - zs(end);
            zstagel = -2*zs(1);            
    end

    
    % create lookup functions that tell you 'n' for a given x value such
    % that x is the nth x value.
    x2i = @(xx) arrayfun(@(x) find(x==xs),xx);
    y2i = @(yy) arrayfun(@(y) find(y==ys),yy);
    z2i = @(zz) arrayfun(@(z) find(abs(z-zs)<1e-6),zz);
    
    % check for datapoint uniformity
    assert(length(xsp)+length(ysp)+length(zsp) == 3,...
        'Non-uniform Datapoint Spacing');
    
    % the size of the 3D data matrices to be filled
    fullsize = [length(xs) length(ys) length(zs)];
    
    % for each x,y,z value in the COMSOL loaded x,y,z columns, check which
    % linear index this corresponds to in a 3D matrix of size fullsize.
    locs = sub2ind(fullsize,x2i(x),y2i(y),z2i(z));
    
    % for each single-letter variable in varname, create a new variable
    % given by a double-letter, which is a 3D matrix of size fullsize, and
    % fill it by indexing into it with the locs linear index column.
    bb=0; ee=0; tt=0; mm=0; xx=0; yy=0; zz=0;
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
    switch(directions)
        case 'noXY'
            fieldsymmetryXY = false;
        otherwise
            fieldsymmetryXY = true; % right now this is broken, set by hand
    end
    
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
    if fieldsymmetryXY
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
    if fieldsymmetryXY
        dvdxm(1,:,:)=0;
        dvdym(:,1,:)=0;
    end
    
    % zero the z force along reflection symmetry in the decelerator:
    if ~strcmp(directions,'load') && ~strcmp(directions,'trap')
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
    function p = getphase(z)
        ii = z/zstagel;
        p = (ii - fix(ii))*360 - 270;
    end
    phase = @getphase;
    
    % If loading or trapping, no wrapping.
    if strcmp(directions,'load') || strcmp(directions,'trap')
        phase = @(z) z/zstagel*360;
    end
    
    % wrap returns a z-coordinate within the force lookup table given a
    % general z-coordinate and the stage parity. It achieves this in two
    % steps. First z is converted to a coordinate between -zstagel and 0
    % which is also the range -270 to +90 in phase angle. Then, the mirror
    % symmetry between the forces in -270--90 and -90-+90 is exploited. The
    % side lookup indicates whether this symmetry is exploited so the
    % z-forces can be inverted, since molecules in the -270 to -90 range
    % are accelerated, not decelerated.
    %wrapc = @(z,n) (phase(z,n)-90)/360 * zstagel;
    %wrapc = @(z) (phase(z)-90)/360 * zstagel;
    %wrap = @(z,n) abs(wrapc(z,n)+zstagel/2)-zstagel/2;
    %side = @(z,n) (wrapc(z,n) > -zstagel/2)*2 - 1;
    
    % let's write a better wrap for shorter runtimes:
    function w = wrapf(z)%,n)
        ii = z/zstagel;%+n/2;
        jj = ii - fix(ii);
        w = (abs(jj-0.5)-0.5);
        w = w*zstagel;
    end
    %wrap = @(z,n) abs( mod(z+zstagel2*n,zstagel)-zstagel2)-zstagel2;
    wrap = @wrapf;
    function s = sidef(z)%,n)
        ii = z/zstagel;%+n/2;
        jj = ii - fix(ii);
        s = (jj > 0.5)*2 - 1;
    end
    %side = @(z,n) (mod(z+n*zstagel2,zstagel) > zstagel2)*2 - 1;
    side = @sidef;
    
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
    if fieldsymmetryXY
        dvdx = @(x,y,z) dvdxg(abs(x),abs(y),wrap(z)).*sign(x);
        dvdy = @(x,y,z) dvdyg(abs(x),abs(y),wrap(z)).*sign(y);
        dvdz = @(x,y,z) dvdzg(abs(x),abs(y),wrap(z)).*side(z);
        vf = @(x,y,z) vfg(abs(x),abs(y),wrap(z));
    else
        dvdx = @(x,y,z) dvdxg(x,y,wrap(z));
        dvdy = @(x,y,z) dvdyg(x,y,wrap(z));
        dvdz = @(x,y,z) dvdzg(x,y,wrap(z)).*side(z);
        vf = @(x,y,z) vfg(x,y,wrap(z));
    end
    if strcmp(directions,'load') || strcmp(directions,'trap')
        dvdx = @(x,y,z) dvdxg(abs(x),abs(y),z).*sign(x);
        dvdy = @(x,y,z) dvdyg(abs(x),abs(y),z).*sign(y);
        dvdz = @(x,y,z) dvdzg(abs(x),abs(y),z);
        vf = @(x,y,z) vfg(abs(x),abs(y),z);
    end    

    % Save in a file for loading and propagating during decel simulation.
    save(['Fields/' fileN '.mat'],'dvdx','dvdy','dvdz',...
        'vf','phase','zstagel','renergy','aenergy');

    %% Produce Output Figure for Debugging
    % There are many potential errors that could be made in the COMSOL
    % output, data input processing. 
    cc = 4e-20;
    
    titleLabel = 'Decelerator ';
    if strcmp(directions,'load')
        titleLabel = 'Loading ';
    elseif strcmp(directions,'trap')
        titleLabel = 'Traping ';
    end
    
    if fieldsymmetryXY
        p = 1;
    else
        p = ceil(size(dvdzm,1)/2);
    end
    
    figure('position',[50,50,1100,1100])
    subplot(2,2,1)
    surf(cap(squeeze(dvdzm(:,p,:)),cc));
    title([titleLabel 'dvdz, X-Z plane, ' fileN '.dat']);

    % You might think the labels are backwards, but they're not. surf uses
    % the second index (the column of the matrix) as the x-axis and the
    % first index (the row) as the y-axis. It makes sense if you think of
    % matrices as oriented with x going left-right and y going up-down, but
    % its crazy when working in 3D.
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    subplot(2,2,2)
    surf(cap(squeeze(dvdzm(p,:,:)),cc));
    title([titleLabel 'dvdz, Y-Z plane, ' fileN '.dat']); 
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    subplot(2,2,3)
    surf(cap(squeeze(dvdxm(:,p,:)),cc));
    title([titleLabel 'dvdx, X-Z plane, ' fileN '.dat']);
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    subplot(2,2,4)
    surf(cap(squeeze(dvdym(p,:,:)),cc));
    title([titleLabel 'dvdy, Y-Z plane, ' fileN '.dat']);
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    figure('position',[100,200,400,400])
    plot(abs(zs),vf(zeros(length(zs),1),zeros(length(zs),1),zs(:)))
    title([titleLabel 'Potential along Z-axis, ' fileN '.dat'])
    xlabel('Z axis')
    ylabel('Potential Energy (J)')
    
    % convenient subfunction for plotting capped surfaces
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
end



