function plotTrapSlice(vvv,cut,selectPlot)

% These will be used to convert the energy to temperature units.
mOH = 2.82328e-26;
kB = 1.381e-23;

accel = 0;

x = (-.95:.05:.95)*1e-3;
zphi = (-3:.05:3)*1e-3;
[xp,yp,zp] = ndgrid(x,x,zphi);

if strcmp(selectPlot,'yz')
    xp = squeeze(yp(cut,:,:));
    zp = squeeze(zp(cut,:,:));
    vv = squeeze(vvv(cut,:,:));
elseif strcmp(selectPlot,'xy')
    xp = squeeze(xp(:,:,cut));
    zp = squeeze(yp(:,:,cut));
    vv = squeeze(vvv(:,:,cut));
end

if length(selectPlot)==2
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
    %title({'Moving Trap Depth',[' a = ' num2str(round(accel/1000)) ' km/s/s, p = ' num2str(phiM+180+phiL-.02) ' deg']},'FontSize',14)
    %title({'Moving Trap Depth',[' a = ' num2str(round(accel/1000)) ' km/s/s, p = ' num2str(phiH-.01) ' deg']},'FontSize',14)
    set(gca,'FontSize',12)
else
    figure
    isosurface(yp,xp,zp,vvv*1e3/kB,200)
    
    
    
    
end