function r = smallstep(r,t)
%% Small simulation step.
% Updates velocity based on acceleration and position based on velocity,
% offset by half of a timestep.
    r.pos = r.pos + r.vel*t/2;
    r.vel = r.vel + r.xtrascale*r.acc(r)*t;
    r.pos = r.pos + r.vel*t/2;
            
    %update time.
    r.time = r.time + t;
end
