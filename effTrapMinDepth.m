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
i = 1;
while length(l)>1
    [~, l] = findpeaks(-z,'MinPeakProminence',i);
    i = i + 1;
    if length(l)==0
        [~, l] = findpeaks(-z,'MinPeakProminence',i-2);
        l = l(1);
    end
end
%figure(123)
%hold on
%plot(z)
%drawnow

potD = potential3D;
mm = min(potD(:));
potD = potD - mm;
MM = max(potD(:));
potD = potD/MM;
potD = potD * 65535;

basin = watershed(potD);

trapLabel = basin(c1,c2,l);

potFindMin = potD;
if any(basin(:)~=trapLabel(1))
    potFindMin(basin~=trapLabel) = 65535;
end
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
edge = false(size(trapBW));
edge(:,:,[1 end]) = true;
edge(:,[1 end],:) = true;
edge([1 end],:,:) = true;
boundary = trapBW & (~trapBWE | edge);

restrictPot = potential3D;
restrictPot(~boundary) = max(restrictPot(:));
[minD,loc] = min(restrictPot(:));

% Let's get volume and PSV too.
trapPot = false(size(potential3D));
trapPot(trapBW) = true;
trapPot(potential3D > minD) = false;
vol = sum(trapPot(:));
sp = 25e-6;
volume = vol*(sp^3);

% And now phase space volume.
% get a matrix with max velocity at each point in the trap.
mOH = 17*1.67e-27;
trapPSV = sqrt(abs(potential3D-minD)*2/mOH);
% replace max velocity with 3-volume in velocity space of all allowed vel's
trapPSV = (4/3)*pi*trapPSV.^3;
trapPSV(~trapPot) = 0;
psVol = sum(trapPSV(:));
psVolume = psVol*(sp^3);

    if nargout <= 1
        varargout = {minD - trueMin};
    elseif nargout > 1    
        [a b c] = ind2sub(size(potential3D),loc);
        [~, lr] = max(abs([a b c] - [c1 c2 l]));
        if lr==3
            minLoc = 'longitudinal';
        else
            minLoc = 'transverse';
        end
        varargout = {minD - trueMin,minLoc};
        if nargout > 2
            varargout{end+1} = volume;
        end
        if nargout ==4
            varargout{end+1} = psVolume;
        end
        if nargout > 4
            error('Too many output arguments');
        end
    end
end