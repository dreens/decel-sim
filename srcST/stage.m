function r = stage(r)
%% The stage function.
% Propagates one decel stage. Removes lost molecules, checks number, etc.

% Update phase in case running in phase ramp mode
    nr = r.numstage - r.stages - 1 + r.rampN + strcmp(r.traploadstyle,'decel');
    thisphase = r.phase + nr * r.phasestep * heaviside(nr) * ...
        strcmp(r.deceltiming,'phaseramp');
    
    if strcmp(r.deceltiming,'bunchatfirst') && r.numstage <= r.skipnum
        thisphase = 0;
    end
        
    if r.tweak
        if mod(r.numstage,2)
            r.xtrascale = r.s12;
        else
            r.xtrascale = r.s34;
        end
    end
    
    if strcmp(r.deceltiming,'mixedmode')
        xtratime = (r.modeseq(r.numstage)-1)*r.f.zstagel/2/r.vel(1,3);
        t=0;
        while t<xtratime
            r = smallstep(r,r.smallt);
            if mod(r.time,r.decelt) < r.smallt
                r = update(r);
            end
            t = t + r.smallt;
        end
    end
    
    if ~strcmp(r.deceltiming,'switchtimes')
    % Step until the synchronous molecule is past the required phase angle.
        while r.f.phase(r.pos(1,3),r.numstage) < thisphase  && r.vel(1,3) > 0
            r = smallstep(r,r.smallt);
            if mod(r.time,r.decelt) < r.smallt
                r = update(r);
            end
        end

        % Iterate forward/backward in time until synchronous molecule is right
        % on top of the correct phase angle, up to r.tinyt.
        undershoot = 1;
        while abs(undershoot) > r.tinyt && r.vel(1,3) > 0

            % get the undershoot in terms of phase
            undershoot = thisphase - r.f.phase(r.pos(1,3),r.numstage);

            % translate to time ignoring acceleration
            undershoot = undershoot *r.f.zstagel/360/r.vel(1,3);

            % step the molecules according to this time
            r = smallstep(r,undershoot);
        end
        r.switchtimes(r.numstage) = r.time;
    else
        while r.time < r.switchtimes(r.numstage)-r.smallt
            r = smallstep(r,r.smallt);
            if mod(r.time,r.decelt) < r.smallt
                r = update(r);
            end
        end
        r = smallstep(r,r.switchtimes(r.numstage)-r.time);
        r = update(r);
    end        
    % If the iteration undoes a datapoint, erase it:
    if r.time < r.times(end)
        r.undo = true;
        r = update(r);
    end
    
    r = update(r);
    if ~r.parallel
        fprintf('Stage %d,\t%dmols\n',r.numstage,r.numleft);
    end
end
