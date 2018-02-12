function r = init(r)
    r = initdecel(r);
    r = initmols(r);
    r = initvars(r);
end

function r = initdecel(r)

    % Load the mat file, or generate it from a COMSOL .dat file if it
    % doesn't exist yet.
    if ~r.reloadfields && exist(['../Decels/' r.decel '.mat'],'file');
        r.f = load(['../Decels/' r.decel '.mat']);
    elseif exist(['../Decels/' r.decel '.dat'],'file');
        processdecelfields(r);
        r.f = load(['../Decels/' r.decel '.mat']);
    else
        error(['File ''slowANDtrap/Decels/' r.decel '.dat'' not found']);
    end
    
    %{
    % z position management variables. The main management principle is
    % that z=0 changes according to what is most convenient for the three
    % different force fields: deceleration, loading, and trapping, but the
    % simulation always maintains a reference to the original zero of the
    % simulation which is the beginning of the first decelerator stage.
    %
    % More specifically, the deceleration zero is the center of the first
    % pin pair, which is a phase angle of -90 degrees with respect to the
    % first stage, defined as the region between the first two pin pairs.
    % In principle the region before the first pin pair could also be used
    % for deceleration, but this would require an extra force lookup table
    % which would be way too much effort for something with minimal impact
    % on the experiment.
    %}
    r.decellength = 0;      % set after decel fields loaded
    r.loadlength = 0;       % set after all fields have been loaded
    r.loadzero = 0;         % offset between decel and loading zeros
    r.trapzero = 0;         % determined from zero definition in COMSOL
    r.offsetz = 0;

    % Make sure this variable exists regardless which kind of timing is
    % used for struct agreement purposes
    r.phasestep = 0;
    
    % If decel style loading, trick timing into thinking there's one fewer
    % stage:
    if strcmp(r.traploadstyle,'decel')
        r.stages = r.stages - 1;
    end
    
    switch r.deceltiming
        case 'finalvz'
            % Choose the phase angle as a function of vfinal, vinitial, and 
            % stage number.
            energyper = .5*r.mOH*(r.initvz^2 - r.finalvz^2)/r.stages;
            
            % Get the xtrascaling in case that is being used
            xs = r.decelvoltage / 12.5 ;

            % Don't let the phase go above 80 to avoid the stage function 
            % missing the synchronous molecule and letting it wrap to -270.
            r.phase = fminbnd(@(phi) (r.f.renergy(phi)*xs-energyper)^2,0,80);
        case 'magdecelspecial'
            % Here we account for the fact that the phase angle chosen in a
            % magnetic pin decelerator will influence the energy that can
            % be removed during loading, because of how close the trap is
            % to the decelerator.
            totenergy = @(phi) r.f.renergy(phi)*(r.stages+1) + ...
                r.f.aenergy(r.loadphase) - r.f.aenergy(phi);
            remenergy = .5*r.mOH*(r.initvz^2 - r.loadvelz^2);
            r.phase = fminbnd(@(phi) (totenergy(phi)-remenergy).^2,0,80);
        case 'phaseangle'
            % Don't do anything. The phase is already set.
        case 'phaseramp'
            if r.loadphase > 87
                lphasecap = 87;
            else
                lphasecap = r.loadphase;
            end
            phases = @(phi) ...
                [phi*ones(1,r.stages-r.rampN), ...
                linspace(phi,lphasecap,r.rampN+1)];
            rensev = @(phi) sum(arrayfun(@(p) r.f.renergy(p),phi));
            totenergy = @(phi) rensev(phases(phi)) - ...
                r.f.aenergy(-phi) + r.f.aenergy(-lphasecap);
            remenergy = .5*r.mOH*(r.initvz^2 - r.loadvelz^2);
            r.phase = fminbnd(@(phi) (totenergy(phi)-remenergy).^2,0,80);
            r.phasestep = (lphasecap - r.phase) / r.rampN ;
        case 'switchtimes'
            if ~isfield(r,'phase')
                error('Must specify the phase of the run whose switchtimes are being used');
            end
            
        case 'bunchatfirst'
            energyper = .5*r.mOH*(r.initvz^2 - r.finalvz^2)/(r.stages - r.skipnum);
            r.phase = fminbnd(@(phi) (r.f.renergy(phi)-energyper)^2,0,80);
        case 'mixedmode'
            if sum(r.modeseq)~=r.stages
                error('The Mode Sequence doesn''t match Decelerator stages')
            end
            energyper = .5*r.mOH*(r.initvz^2 - r.finalvz^2)/length(r.modeseq);
            r.phase = fminbnd(@(phi) (r.f.renergy(phi)-energyper)^2,0,80);
    end
    
    % Remember the total decelerator length
    r.decellength = r.f.zstagel * r.stages / 2;
    r.loadzero = r.decellength;
    
    % If the decelerator is doing the loading, correct stage number
    if strcmp(r.traploadstyle,'decel')
        r.stages = r.stages + 1;
    end
    
    % Remember stage length as well
    r.zstagel = r.f.zstagel;
end

function r = initmols(r)
%% Initialize Molecules
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
    
    if r.rememberwho
        r.veli = r.vel;
        r.posi = r.pos;
    end
    if r.rememberwhoall
        r.veliall = r.veli;
        r.posiall = r.posi;
    end
    
    %Shift molecules in z, vz.
    r.vel(:,3) = r.vel(:,3) + r.initvz;
    r.pos(:,3) = r.pos(:,3) - r.vdd;
end

function r = initvars(r)
    % Stage Number
    r.numstage = 1;

    % Store the molecule number each decel stage.
    r.molnum = [];
    
    % Store the number in the laser- applicable only during detection.
    r.molnumlaser = [];
    r.molsatmaxinlaser = [];
    
    % Store the synchronous molecule's information.
    r.zsynch = [];
    r.vzsynch = [];
    r.ksynch = zeros(0,2);
    r.pesynch = [];
    
    % This one is used for removing a set of synch molecule info.
    r.undo = false;
    
    % Extra scaling can be applied to simulate tuning load voltage.
    r.xtrascale = r.decelvoltage/12.5;
    
    % Information can be stored at different time spacings depending on
    % which part of the simulation, deceleration, loading, or trapping, is
    % active. For this reason we also note the time of each storage so that
    % data can be plot against time accurately.
    r.times = [];
    
    % This prevents erors in the first number storage
    r.numleft = r.num;
    
    
    % Create the switchtime matrix if it isn't set in simdeceltrap
    if ~isfield(r,'switchtimes')
        r.switchtimes = [];
    end
    
    % Store the full phase space each stage.
    %r.phasestage = zeros(r.num,6,r.stages);

    % other variables used by different runs
    %r.lost = false(r.num,1);
    
    % The data name shouldn't have spaces in it.
    
    
    % Timing
    r.time = 0;
    r.decelofftime = inf;
    r.loadingofftime = inf;
end
