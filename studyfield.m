function h = studyfield(d,phi)

r = load(['Decels/' d '.mat']);
mOH = 2.82328e-26;
kB = 1.381e-23;

z=(-15:.1:5)*1e-3;
x = (-.9:.1:.9)*1e-3;
[xx,zz] = meshgrid(x,z);
ff = zeros(size(xx));
gg = zeros(size(xx));
ff(:) = r.vf([zeros(length(xx(:)),1) xx(:) zz(:)],2);
gg(:) = r.vf([xx(:) zeros(length(xx(:)),1) zz(:)],2);
zphi = (-2:.1:2)*1e-3;
[xp,zp] = meshgrid(x,zphi);
vv = zeros(length(zphi),length(x));
for i=1:length(zphi)
    phiL = -7.5e-3 + phi/90*2.5e-3;
    phiH = phiL + 5e-3;
    [~, a] = min((zphi(i)+phiL-z).^2);
    [~, b] = min((zphi(i)+phiH-z).^2);
    b = b-1;
    vv(i,:) = mean(ff(a:b,:)) + mean(gg(a:b,:));
end

% Now velocity compensation
mZ = (length(zphi)+1)/2;
mX = (length(x)+1)/2;
vv = vv - zp.*(vv(mZ+1,mX)-vv(mZ-1,mX))/(zp(mZ+1,mX)-zp(mZ-1,mX));
vv = vv - vv(mZ,mX);
m = max(vv)/2;
h = figure('Position',[50 50 500 800]);
contourf(xp*1e3,zp*1e3,1e3*vv/kB,0:10:200);
g = colorbar;
g.YLabel.String = 'Energy (K)';
g.YLabel.Rotation = 270;
xlabel('Transverse Position (mm)')
ylabel('Longitudinal Position (mm)')
title({'Moving Trap Depth',[' Phi=' num2str(phi)]})

end