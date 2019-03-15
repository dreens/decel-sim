%%
%
% effTrapMinDepth.m
%
% Find the minimum of an effective moving trap.
function varargout = effTrapMinDepth(potential3D)

% first find trap center. Should be only local min along z axis.
c1 = (size(potential3D,1)+1)/2;
c2 = (size(potential3D,2)+1)/2;
z = potential3D(c1,c2,:);
z = squeeze(z);
[~, l] = findpeaks(-z);

potD = potential3D;
mm = min(potD(:));
potD = potD - mm;
MM = max(potD(:));
potD = potD/MM;
potD = potD * 65535;

basin = watershed(potD);

trapLabel = basin(c1,c2,l);

potFindMin = potD;
potFindMin(basin~=trapLabel) = 65535;
[~, indTrueMin] = min(potFindMin(:));
trueMin = potential3D(indTrueMin);

trapBW = false(size(basin));
trapBW(basin==trapLabel) = true;
d3 = false(3,3,3);
d3(2,2,[1 3]) = true;
d3(2,[1 3],2) = true;
d3([1 3],2,2) = true;
strd3 = strel('arbitrary',d3);
trapBWE = imerode(trapBW,strel('cube',3));
boundary = trapBW & ~trapBWE;

potential3D(~boundary) = max(potential3D(:));
[minD,loc] = min(potential3D(:));
if nargout <= 1
    varargout = {minD - trueMin};
elseif nargout == 2    
    [a b c] = ind2sub(size(potential3D),loc);
    [~, lr] = max(abs([a b c] - [c1 c2 l]));
    if lr==3
        minLoc = 'longitudinal';
    else
        minLoc = 'transverse';
    end
    varargout = {minD - trueMin,minLoc};
else
    error('Too many output arguments');
end
end