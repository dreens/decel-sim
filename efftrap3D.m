% Here I attempt to create the effective moving trap potential under
% different alternate charging sequences of interest.

function vv = efftrap3D(d,e,phiL,phiM,phiH)

% Shift phi a tiny bit to avoid strange problem that occurs if phi/90 is a
% simple rational.
phiL = phiL + .01;
phiM = phiM+ .01;
phiH = phiH+ .01;

% This allows this code to be used for a different case:
if strcmp(d,'None')
    d = e;
    zero = true;
else
    zero = false;
end

% Load the two fields that will be combined in different phase angle ranges
% to give the effective moving trap.
r = load(['Decels/' d '.mat']);
s = load(['Decels/' e '.mat']);

% These will be used to convert the energy to temperature units.
mOH = 2.82328e-26;
kB = 1.381e-23;

% This is the range of coordinates that will be included. The z range is
% large because different chunks get integrated over to get the effective
% moving trap.
z=(-15:.025:5)*1e-3;

if phiH-phiL>180
    z=(-20:.025:5)*1e-3;
end

x = (-.975:.025:.975)*1e-3;
[xx,yy,zz] = ndgrid(x,x,z);

% ff, gg, hh, ii will be large x-z grids giving the lab-fixed potential
% energy over a few stages under the two different charge configurations
% and under their two different possible translations. (e.g. rods 1/2
% charged verses rods 3/4 charged)
ff = zeros(size(xx));
gg = zeros(size(xx));
hh = zeros(size(xx));
ii = zeros(size(xx));

% Of course the separate possible translations are achieved by changing x
% for y as appropriate.
ff(:) = r.vf([yy(:) -xx(:) zz(:)],2);
gg(:) = r.vf([xx(:) yy(:) zz(:)],2);
hh(:) = s.vf([yy(:) -xx(:) zz(:)],2);
ii(:) = s.vf([xx(:) yy(:) zz(:)],2);

% Clear second fields
if zero
    ff = ff*0;
    gg = gg*0;
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
    pH = -2.5e-3 + phiH/90*2.5e-3;
    pM = -2.5e-3 + phiM/90*2.5e-3;
    pL = -2.5e-3 + phiL/90*2.5e-3;

    [~, a] = min((zphi(i)+pL-z).^2);
    [~, b] = min((zphi(i)+pM-z).^2);
    [~, c] = min((zphi(i)+pH-z).^2);
    vv(:,:,i) = (sum(ff(:,:,b:(c-1)),3) + sum(gg(:,:,b:(c-1)),3) + ...
        sum(hh(:,:,a:(b-1)),3) + sum(ii(:,:,a:(b-1)),3))/(c-a+1)/2; %beware sensitivity on 
    
end
    pH = -2.5e-3 + phiH/90*2.5e-3;
    pM = -2.5e-3 + phiM/90*2.5e-3;
    [~, bb] = min((z-pM).^2);
    [~, cc] = min((z-pH).^2);
slope2=(ff(40,40,cc)-ff(40,40,bb))/(5e-3+1e-8);
acel2=slope2/mOH;
% Symmetrize for pggg. In principle I should add in the voltages
% corresponding to gmgg etc, but I know that the net result will just be
% symmetrizing over the axis, so why not do it directly.
vv = (vv + vv(end:-1:1,end:-1:1,:))/2;

% Now velocity compensation. First get coordinates of what should be the
% trap center.
mZ = (length(zphi)+1)/2;
mX = (length(x)+1)/2;

% pHh = -2.5e-3 + phiH/90*2.5e-3;
%  [~, mZ] = min((zphi-pHh).^2);
% mZ



% Find the slope of the potential energy at the "trap center". Now add a
% fictitious force corresponding to acceleration enough so this slope is
% zeroed.
slope = (vv(mX,mX,mZ+1)-vv(mX,mX,mZ-1))/(zp(mX,mX,mZ+1)-zp(mX,mX,mZ-1));
vv = vv - zp.*slope2;
   % vv = vv - (zp-2.5e-3).*slope;
% We can also get the acceleration implied by that slope.
accel = slope/mOH;

% Finally redefine zero energy.
vv = vv - vv(mX,mX,mZ);

end