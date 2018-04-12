% Here I attempt to create the effective moving trap potential under
% different alternate charging sequences of interest.

function h = study2fields(d,e,phi)

r = load(['Decels/' d '.mat']);
s = load(['Decels/' e '.mat']);
mOH = 2.82328e-26;
kB = 1.381e-23;

z=(-15:.1:5)*1e-3;
x = (-.9:.1:.9)*1e-3;
[xx,zz] = meshgrid(x,z);
ff = zeros(size(xx));
gg = zeros(size(xx));
hh = zeros(size(xx));
ii = zeros(size(xx));
ff(:) = r.vf([zeros(length(xx(:)),1) xx(:) zz(:)],2);
gg(:) = r.vf([xx(:) zeros(length(xx(:)),1) zz(:)],2);
hh(:) = s.vf([zeros(length(xx(:)),1) xx(:) zz(:)],2);
ii(:) = s.vf([xx(:) zeros(length(xx(:)),1) zz(:)],2);

zphi = (-2:.1:2)*1e-3;
[xp,zp] = meshgrid(x,zphi);
vv = zeros(length(zphi),length(x));
for i=1:length(zphi)
    phiH = -2.5e-3 + phi/90*2.5e-3;
    phiM = -5e-3 + (1-phi/90)*2.5e-3;
    phiL = phiH - 5e-3;
    [~, a] = min((zphi(i)+phiL-z).^2);
    [~, b] = min((zphi(i)+phiM-z).^2);
    [~, c] = min((zphi(i)+phiH-z).^2);
    vv(i,:) = (sum(ff(b:(c-1),:)) + sum(gg(b:(c-1),:)) + ...
        sum(hh(a:(b-1),:)) + sum(ii(a:(b-1),:)))/(c-a+1)/2;
    
end

% Symmetrize for pggg
vv = (vv + vv(:,end:-1:1))/2;

% Now velocity compensation
mZ = (length(zphi)+1)/2;
mX = (length(x)+1)/2;
vv = vv - zp.*(vv(mZ+1,mX)-vv(mZ-1,mX))/(zp(mZ+1,mX)-zp(mZ-1,mX));
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
title({'Moving Trap Depth',[' Phi=' num2str(phi)]},'FontSize',14)
set(gca,'FontSize',12)

end