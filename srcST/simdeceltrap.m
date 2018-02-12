function rsf = simdeceltrap()

    cd('~/Documents/MATLAB/slowANDtrap/src')

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
    r.dname = 'wide 311';
    r.num = 1e4;
    r.tempxy = 3000e-3;
    r.spreadxy = 3e-3;
    r.tempz =  4000e-3; 
    r.spreadz =  4e-3;
    r.initvz = 460;%{480 620 760};
    r.dist = 'homogeneous';
    r.rememberwho = true;
    r.rememberwhoall = false;
    
    % set this to force field reloading. 
    r.reloadfields = false;
    
    % Set these functions to plot or fit sim results.
    r.plotdata = true;
    r.plotfunc = 'results_synch';
    r.fitdata  = false;
    r.fitfunc = 'fitsim2data';
    
    % decelerator configuration variables
    r.stages = 143;%{600,600,300,300,143,143};
    r.vdd = 20e-3;
    r.decel = {'normalwide'};%('electrodering','uniformmagnet','normal','magneticpin','normalwide'};
    %r.loadname = 'tricycleload14p5';%{'tricycleload14p5','tricycleload14p5g'};
    r.loadname = 'magneticpinload';
    %r.trapname = 'tricycletrap_z0';
    r.trapname = 'magneticpintrap';
    
    % decelerator timing variables
    % Choose finalvz, magdecelspecial, phaseangle, phaseramp, 
    %     switchtimes, skipload, bunchatfirst, mixedmode;
    r.deceltiming = 'mixedmode';
    r.decelvoltage = 12.5;
    r.switchtimes = 0; %st;  % used for loading other switchtimes
    r.phase = 0; % used if timing by phase angle
    r.finalvz = 35;%repmat({35 40 45},1,5); % not used for magdecelspecial or phaseramp
    %r.finalvz = {480 426 365 289 40}; %for 142 stages 480m/s, 0-10-20-30-working phaseangle
    %r.finalvz = {380 20 0};
    %r.finalvz = {480 390 273 40};%for 227 stages 480m/s, 0-10-20-working phaseangle
    %r.finalvz = {620 554 478 385.5 40}; %for 227 stages 620m/s, 0-10-20-30-working phaseangle
    r.rampN = 1;%{1 2 3 4 6 10 15 25 50}; % How many stages to ramp the phase
    r.skipnum = 0;%{0 50 100 150};
    %r.modeseq = {ones(1,600), 3*ones(1,200)};
    %r.modeseq(3:4) = {ones(1,300),3*ones(1,100)};
    r.modeseq = {ones(1,143), [repmat([1 1 3],1,20) ones(1,43)]};
    
    % Decelerator Tweak Testing
    r.tweak = false;
    r.s12 = 1;%{.7 .8 .9 1 1.1 1.2 1.3};
    r.s34 = r.s12;
    
    % Set these to run only part of the slow & trap simulation
    r.slowonly = true;
    r.loadonly = false;  % useful for creating loaded molecule datasets.
    
    % loading configuration variables
    r.traploadstyle = 'loadvelz'; %{'timed','loadvelz','distance','decel'}
    r.loadtime = 420e-6;%num2cell((340:100:400)*1e-6);
    r.loaddist = 8e-3;
    r.loadvelz = 5;%{0 5 10 15 20 25 30 35 40 45};
    r.loadphase = 75;%{50 60 70};%{50 55 60 65 70 75 80 85 90 92 94 96 98 100};
    r.addflip = true;%{false true};
    r.voltagescaling = 1;
    
    % simulation timing variables
    r.tinyt = 1e-10;        % accuracy of decel switching
    r.smallt = 2e-7;        % propagator timescale
    r.decelt = 5e-8;        % data update timescale during deceleration
    r.loadt = 1e-6;         % ditto, during loading
    r.trapt = 1e-3;         % ditto, during trapping
    r.traptime = 200e-3;    % total time in the trap
    r.detecttime = 500e-6;  % total time detecting
    r.detectt = 1e-5;       % data update during detection
    r.decelofftime = 0;     % will be set when decel turns off
    r.loadingofftime = 1;   % will be set when trapping begins
    
    % laser beam variables
    r.lasertype = 'disk'; % may implement gaussian or elliptical
    r.LBD = 2.2e-3;
    r.detection = 'slowing'; %{'in situ', 'slowing', 'free flight'}
    
    % This is used for slowing detection
    r.laserpos = .79;

    
    % random number seed.
    r.doseed = false;
    r.seed = 100;
        
    % Set for cluster submissions
    r.parallel = false;
    r.subparallel = true;
    r.cores = 16;
    
    % Choose unpacking type (must precede overwriting, but unpacking itself
    % must follow overwriting)
    r.packtype = 'product';

    % Unpack r into a struct of runs
    rs = unpacker(r);


    % Loop through each run, save data call processing functions
    if ~r.parallel
        rsf = runall(rs);
    else
        jobs = 4;
        runs = length(rs);
        eachjob = runs/jobs;
        for i=1:jobs
            these = rs(1+(i-1)*eachjob:i*eachjob);
            mqsub('runall',{these},'Name',[r.dname '_' num2str(i)],...
                'Cores',16,'Queue','hexadec');
        end
        rsf = 1;
    end    
    
end

