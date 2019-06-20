% Select Plot can select a 2D contour if in the form of a 2 digit string
% such as 'xy'. It can also be 'contourArray' or 'contours'. More rectnly,
% 1D options were implemented under '1D'.
function h = plotTrapSlice(vvv,cut,selectPlot)

% These will be used to convert the energy to temperature units.
mOH = 2.82328e-26;
kB = 1.381e-23;

accel = 0;

x = (-.975:.025:.975)*1e-3;
zphi = (-3:.025:3)*1e-3;
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

if length(selectPlot)==2 && ~strcmp(selectPlot,'1D')
    h = figure('Position',[50 50 500 800]);
    contourf(xp*1e3,zp*1e3,1e3*vv/kB,0:10:800);
    caxis([0 500]);
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
elseif strcmp(selectPlot,'1D')
    h = figure('Position',[50 50 800 800]);
    cx = (length(x)+1)/2;
    cz = (length(zphi)+1)/2;
    v = 1e3*vvv/kB;
    plot(x,squeeze(v(:,cx,cz)),'b','DisplayName','x-axis')
    hold on
    vv = x;
    for i=1:length(x)
        vv(i) = v(i,i,cz);
    end
    plot(x*sqrt(2),vv,'k--','DisplayName','w-axis')
    grid on
    plot(zphi,squeeze(v(cx,cx,:)),'r','DisplayName','z-axis')
    xlabel('Position on Axis through Center (mm)')
    ylabel('Energy (mK)')
    title('Effective Trap Potential along Axes')
    
elseif strcmp(selectPlot,'contours')
    f = figure('Position',[300 50 300 400]);
    a = axes('Parent',f);
    levels = [10,50,150,250];
    colors = jet(length(levels));
    for i=1:length(levels)
        makepatch(levels(i),colors(i,:),97);
    end
    view(a,[15 15]);
elseif strcmp(selectPlot,'contourArray')
    
    x = (-.975:.025:.975)*1e-3;
    zphi = (-3:.025:3)*1e-3;
    [zp,yp,xp] = ndgrid(zphi,x,x);
    

    
    f = figure('Position',[10 10 1116 570]);
    a = axes('Parent',f);
    vvvfam = vvv;
    levels = [10,55,100,200];
    %colors = jet(length(levels));
    %colors = flipud(colors);
    % these are the fixed brightness new matlab colors
    colors = [0 0.447 0.741 ; 0.85 0.325 0.098 ; 0.929 0.694 0.125 ; 0.494 0.184 0.556 ; 0.446 0.674 0.188];
    colors = flipud(colors);
    for types=1:length(vvv)
        vvv = vvvfam{types};
        vvv = permute(vvv,[3 2 1]);
        xp = xp + ones(size(yp))*2.5e-3;
        for views = 1:length(levels)
            zp = zp + ones(size(zp))*4.5e-3*views;
            if types==4
                cut = 303;
            else
                cut = 193;
            end
            makepatch(levels(views),colors(types,:),cut);
            zp = zp - ones(size(zp))*4.5e-3*views;
        end
    end
    xlim([-2e-3 2e-3]);
    view(a,[80.2 5.2]);
    light('Position',[-5 5 -1]);
    light('Position',[0 5 0]);
    light('Position',[-5 -2 -5]);
   % light('Position',[-5 -10 -1]);
end

function h = makepatch(t,c,cut)
    cut = min(cut,max(size(xp)));
    if t<=100
        cut = 191;
    end
    if ~strcmp(selectPlot,'contourArray')
        h = patch(isosurface(yp(:,:,1:cut),xp(:,:,1:cut),zp(:,:,1:cut),vvv(:,:,1:cut)*1e3/kB,t));
    else
        h = patch(isosurface(yp(1:cut,:,:),zp(1:cut,:,:),xp(1:cut,:,:),vvv(1:cut,:,:)*1e3/kB,t));
        isonormals(yp(1:cut,:,:),zp(1:cut,:,:),xp(1:cut,:,:),vvv(1:cut,:,:)*1e3/kB, h);
    end
    h.FaceColor = c;
    h.EdgeAlpha = 0;
    h.FaceAlpha = 1;
    h.AmbientStrength = 0.1;
    h.DiffuseStrength = .5;
    h.SpecularStrength = .25;
    h.BackFaceLighting = 'unlit';
end

end