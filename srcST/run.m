function r = run(r)
%% This function does a run.
% It times itself and calls the step function.
    r = tofirststage(r);
    
    if strcmp(r.deceltiming,'mixedmode')
        limit = @(num) num<=length(r.modeseq);
    else
        limit = @(num) num<=r.stages;
    end
        
    while limit(r.numstage)
        if ~mod(r.numstage,2)
            r.acc = @accn;
        else
            r.acc = @accr;
        end
        r = stage(r);
        r.numstage = r.numstage + 1;
    end
        
    %{
    % Best practice is to ensure that the loading fields, which include the
    % final stage of deceleration, are oriented in the same manner as the
    % deceleration fields when generating in COMSOL. If this is done, then 
    % no flip will be needed. To be more clear, COMSOL is used to generate
    % a "stage" defined as -90 degrees to +90 degrees of phase, or between
    % the centers of two pin pairs. The loading fields include this same
    % stretch, with a trap attached, and with different fields applied. The
    % two pin pairs in loading should be oriented the same way as the two
    % pairs in deceleration.
    %}
    if xor( r.addflip , strcmp(r.traploadstyle,'decel') )
        if ~mod(r.numstage,2)
            r.acc = @accn;
        else
            r.acc = @accr;
        end
    end
    
    % If we're looking at slowing traces, just run detection at this point.
    if r.slowonly
        r.decelofftime = r.time;
        r.laserpos = r.pos(1,3) + 9e-3;
        r = detect(r);
        return
    end
    
    r = initloading(r);
    r = trapload(r);
    r = inittrap(r);
    
    % If we're only doing loading so as to get molecule locations, let
    % inittrap run first so as to get sensible z position offsets.
    if r.loadonly
        return
    end
    
    r = trap(r);
    r = detect(r);
end
