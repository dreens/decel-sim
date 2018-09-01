% get final speeds from files labeled by phase angle
vels = zeros(1,16);
for f = 50:2:80
    name = sprintf('Partials/Bnorm56p%d',f);
    load(name);
    vels(f/2-24)=r.vel(1,3);
end
vels