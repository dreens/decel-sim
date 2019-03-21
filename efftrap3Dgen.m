% Here I modify efftrap3D to accept a more generalized list of charge
% configurations and phase angles for switching between them. This could be
% used to get effective traps for DongDong's advanced switching schemes or
% for checking alternatives to VSF mode, etc.

function varargout = efftrap3Dgen(decels,phis,primes,varargin)

if ~isempty(varargin)
    acc = varargin{1};
    warning('4th Argument Found, Forcing a Deceleration rate.')
    assert(isreal(acc),'Acceleration must be a real number.')
else
    acc = nan;
end

assert(length(decels)+1==length(phis),'You need one phase for each charge configuration plus an extra.');
assert(length(primes)==length(decels),'Specify a prime for each charge config');

% Shift phi a tiny bit to avoid strange problem that occurs if phi/90 is a
% simple rational.
phis = phis + .01;

% Load the fields that will be combined in different phase angle ranges
% to give the effective moving trap.
ds(1:length(decels)) = struct();
for i=1:length(decels)
    if ~strcmp(decels{i},'None')
        ds(i).ff = load(['Decels/' decels{i} '.mat']);
    else
        ds(i).ff.vf = @(x,~) zeros(max(size(x)),1);
    end
end

% These will be used to convert the energy to temperature units.
mOH = 2.82328e-26;
kB = 1.381e-23;

% This is the range of coordinates that will be included. The z range is
% large because different chunks get integrated over to get the effective
% moving trap.
z=(-25:.025:10)*1e-3;

x = (-.975:.025:.975)*1e-3;
[xx,yy,zz] = ndgrid(x,x,z);

% lfgrids will contain large x-z grids giving the lab-fixed potential
% energy over a few stages under the various charge configurations
% specified.
lfgrids(1:length(decels)) = struct('vv',zeros(size(xx)));
for i=1:length(lfgrids)
    if primes(i)
        lfgrids(i).vv(:) = ds(i).ff.vf([xx(:) yy(:) zz(:)],2);
    else
        lfgrids(i).vv(:) = ds(i).ff.vf([yy(:) -xx(:) zz(:)],1);
    end
end

% Now we make a new z coordinate centered on the synchronous molecule. We
% will fill the variable vv with effective potential energy relative to
% this molecule.
zphi = (-3:.025:3)*1e-3;
[xp,yp,zp] = ndgrid(x,x,zphi);
vv = zeros(size(xp));

% For each z effective coordinate, we will integrate over lab space z for a
% 5 mm range with differing start and end points. The idea is that when a
% molecule is ahead of the synchronous one in z, it experiences the
% relevant charge configurations over a different range of lab space
% coordinates, different in start and end point but not length. All
% assuming the z velocity is large relative to the stage spacing.
for i=1:length(zphi)
    ps = -2.5e-3 + phis/90*2.5e-3;

    as = ps;
    for j=1:length(phis)
        [c, as(j)] = min((zphi(i)+ps(j)-z).^2);
        if c > (mode(diff(z))/2)^2
            error('Add more z-points. angles out of range.')
        end
        if j>1
            vv(:,:,i) = vv(:,:,i) + sum(lfgrids(j-1).vv(:,:,as(j-1):as(j)-1),3);
        end
    end
    vv(:,:,i) = vv(:,:,i)/(as(end)-as(1));
end

% Symmetrize for pggg. In principle I should add in the voltages
% corresponding to gmgg etc, but I know that the net result will just be
% symmetrizing over the axis, so why not do it directly.
vv = (vv + vv(end:-1:1,end:-1:1,:))/2;

% Symmetrize LR
vv = (vv + permute(vv,[2 1 3]))/2;

% Now velocity compensation. First get coordinates of what should be the
% trap center.
mZ = (length(zphi)+1)/2;
mX = (length(x)+1)/2;

% Find the slope of the potential energy at the "trap center". Now add a
% fictitious force corresponding to acceleration enough so this slope is
% zeroed.
slope = (vv(mX,mX,mZ+1)-vv(mX,mX,mZ-1))/(zp(mX,mX,mZ+1)-zp(mX,mX,mZ-1));
vv = vv - zp.*slope;

% We can also get the acceleration implied by that slope.
accel = slope/mOH;
    
if ~isnan(acc)
    vv = vv + zp.*slope;
    vv = vv - zp.*acc*mOH;
    accel = acc;
end

% Finally redefine zero energy.
vv = vv - vv(mX,mX,mZ);

if nargout == 1
    varargout = {vv};
elseif nargout == 2
    varargout = {vv, accel};
else
    error('Too many output arguments expected')
end

end