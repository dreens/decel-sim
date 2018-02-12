function r = trapload(r)
%% Molecule Loading
% During loading, the molecules are stepped according to the force-fields
% loaded by processloadingfields.m and initialized by initloading.m.
% Different conditions may be set to control the endpoint of the loading
% phase.

    fprintf('Loading\n')
    
    % s is the step number, which counts loading time steps of length
    % r.loadt. This time-scale controls updating of the logging arrays used
    % for plotting.
    s = 1;
    isdone = false;
    while ~isdone
        [isdone, r] = done(r);
        r = smallstep(r,r.smallt);
        if mod(r.time,r.loadt)<r.smallt
            r = update(r);
            fprintf('Step %d,\t%dmols\n',s,r.numleft); s = s+1;
            if isnan(r.pos(1,3))
                break
            end
        end
    end
end

function [isdone, r] = done(r)
%% Style Selection
% Loading can be based on absolute time, on distance traveled by the
% synchronous molecule, or by final z velocity of that molecule.
    switch(r.traploadstyle)
        case 'decel'
            isdone = true; % in this case 'loading' already achieved.
            if r.loadphase > 87
                r.traploadstyle = 'timed';
                timeadd = (r.loadphase - 87)*r.zstagel/360/20;
                r.loadtime = r.time - r.decelofftime + timeadd;
                isdone = false;
            end            
        case 'timed'
            isdone = r.time >= r.decelofftime + r.loadtime;
        case 'distance'
            isdone = r.pos(1,3) >= r.loaddist;
        case 'loadvelz'
            isdone = abs(r.vel(1,3) - r.loadvelz) < .1;
            
            % This condition implements a kind of analytic continuation of
            % the notion of loading phase to include loading timings where
            % the synchronous molecule is allowed to go past the stopping
            % peak, i.e. to phase angles higher than 90 degrees.
            if isdone && r.loadphase > 87
                r.traploadstyle = 'timed';
                timeadd = (r.loadphase - 87)*r.zstagel/360/20;
                r.loadtime = r.time - r.decelofftime + timeadd;
                isdone = false;
            end
    end
end