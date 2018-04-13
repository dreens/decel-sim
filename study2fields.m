% Here I attempt to create the effective moving trap potential under
% different alternate charging sequences of interest.

function h = study2fields(d,e,phiL,phiM,phiH)

% Shift phi a tiny bit to avoid strange problem that occurs if phi/90 is a
% simple rational.
phiL = phiL + .01;
phiM = phiM + .01;
phiH = phiH + .01;

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
z=(-15:.1:5)*1e-3;
x = (-.9:.1:.9)*1e-3;
[xx,zz] = meshgrid(x,z);

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
ff(:) = r.vf([zeros(length(xx(:)),1) xx(:) zz(:)],2);
gg(:) = r.vf([xx(:) zeros(length(xx(:)),1) zz(:)],2);
hh(:) = s.vf([zeros(length(xx(:)),1) xx(:) zz(:)],2);
ii(:) = s.vf([xx(:) zeros(length(xx(:)),1) zz(:)],2);

% Clear second fields
if zero
    ff = ff*0;
    gg = gg*0;
end

% Now we make a new z coordinate centered on the synchronous molecule. We
% will fill the variable vv with effective potential energy relative to
% this molecule.
zphi = (-2:.1:2)*1e-3;
[xp,zp] = meshgrid(x,zphi);
vv = zeros(length(zphi),length(x));

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
    vv(i,:) = (sum(ff(b:(c-1),:)) + sum(gg(b:(c-1),:)) + ...
        sum(hh(a:(b-1),:)) + sum(ii(a:(b-1),:)))/(c-a+1)/2;
    
end

% Symmetrize for pggg. In principle I should add in the voltages
% corresponding to gmgg etc, but I know that the net result will just be
% symmetrizing over the axis, so why not do it directly.
vv = (vv + vv(:,end:-1:1))/2;

% Now velocity compensation. First get coordinates of what should be the
% trap center.
mZ = (length(zphi)+1)/2;
mX = (length(x)+1)/2;

% Find the slope of the potential energy at the "trap center". Now add a
% fictitious force corresponding to acceleration enough so this slope is
% zeroed.
slope = (vv(mZ+1,mX)-vv(mZ-1,mX))/(zp(mZ+1,mX)-zp(mZ-1,mX));
vv = vv - zp.*slope;

% We can also get the acceleration implied by that slope.
accel = slope/mOH;

% Finally redefine zero energy.
vv = vv - vv(mZ,mX);
m = max(vv)/2;
h = figure('Position',[50 50 500 800]);
contourf(xp*1e3,zp*1e3,1e3*vv/kB,0:10:800);
caxis([0 500])
g = colorbar;
g.YLabel.String = 'Energy (mK)';
g.YLabel.Rotation = 270;
g.YLabel.Position = g.YLabel.Position + [1 0 0];
g.YLabel.FontSize = 12;
xlabel('Transverse Position (mm)','FontSize',12)
ylabel('Longitudinal Position (mm)','FontSize',12)
%title({'Moving Trap Depth',[' Phi=' num2str(phiH-.01)]},'FontSize',14)
title({'Moving Trap Depth',[' a = ' num2str(round(accel/1000)) ' km/s/s, p = ' num2str(phiM+180+phiL-.02) ' deg']},'FontSize',14)
%title({'Moving Trap Depth',[' a = ' num2str(round(accel/1000)) ' km/s/s, p = ' num2str(phiH-.01) ' deg']},'FontSize',14)
set(gca,'FontSize',12)

end