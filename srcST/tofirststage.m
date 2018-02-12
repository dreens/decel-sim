function r = tofirststage(r)
%% Propagate Molecules to first stage
% Eventually this could be modified to include hexapole focusing or maybe
% magnetic quadrupole focusing.

    % z=0 for deceleration is -90 degrees of the first stage. Since the
    % first stage exists between the first two pin pairs, -90 is the
    % geometry center of the very first pin pair.
    %
    % Thus one might think the distance to propagate the molecules should
    % be to z=0. However, to simplify calculations of the total energy
    % removed by the decelerator, it is useful to turn on the
    % decelerator exactly when the synchronous molecule is where it would
    % be if the nonexistent previous stage had just turned off. Otherwise
    % the first stage actually removes a tad bit more energy than all the
    % others.
    %
    % The z-coordinate when the decelerator is first turned on should thus
    % be the phase angle, -180 degrees because its the previous stage, +90
    % degrees because z=0 is actually at -90. Scaled by r.f.zstagel, the
    % length of 2 full stages, divided by 360 degrees.
    %
    % Note also that in init.m the molecules are started off at a z
    % coordinate of -r.vdd (valve-decel-distance).
    
    distance = (r.phase - 180 + 90)*r.f.zstagel/360-r.pos(1,3);
    time = distance/r.vel(1,3);
    r.pos = r.pos + r.vel*time;
end
