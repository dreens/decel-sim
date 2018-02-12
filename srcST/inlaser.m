function in = inlaser(r)
%% Checks if molecules are in the laser

if r.time > r.loadingofftime + r.traptime
    in = r.pos(:,2).^2 + (r.pos(:,3) - 5.461e-3/2).^2 < (r.LBD/2)^2;

    in = in | (r.pos(:,2).^2 + (r.pos(:,3) + 5.461e-3/2).^2 < (r.LBD/2)^2);
elseif r.slowonly && r.time > r.decelofftime
    in = r.pos(:,2).^2 + (r.pos(:,3) - r.laserpos).^2 < (r.LBD/2)^2;
else
    in = false;
end

end