%%
%
% effTrapMinDepth.m
%
% Find the minimum of an effective moving trap.
function minD = effTrapMinDepth(potential3D)

% first find trap center. Should be only local min along z axis.
c1 = (size(potential3D,1)+1)/2;
c2 = (size(potential3D,2)+1)/2;
z = potential3D(c1,c2,:);

basin = watershed(potential3D);
[m l] = min(potential3D(:));
[a b c] = ind2sub(size(potential3D),l);
l = basin(a,b,c);
end