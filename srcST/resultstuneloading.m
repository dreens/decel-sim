%% Plots results
function resultstuneloading(rs)
    loadz =  zeros(length(rs),1);
    loadvz = zeros(length(rs),1);
    rampN =  zeros(length(rs),1);
    number = zeros(length(rs),1);
    
    for i=1:length(rs)
        r = rs(i);
        loadz(i) =  r.loadphase;
        loadvz(i) = r.loadvelz;
        rampN(i) =  r.rampN;
        number(i) = r.molnum(end);
    end
    
    loadzs = unique(loadz);
    loadvzs = unique(loadvz);
    rampNs = unique(rampN);
    
    numberlookupsingle = @(d,l,t) number(d==loadz & l==loadvz & t==rampN);
    numberlookup = @(d,l,t) arrayfun(numberlookupsingle,d,l,t);
    
    minr = min(number); maxr = max(number);
    
    steporder = round(log10(maxr-minr))-1;
    
    minr = round(minr,-steporder);
    maxr = round(maxr,-steporder);
    
    lines = minr:10^steporder:maxr;
    
    [zs, vzs] = meshgrid(loadzs,loadvzs);
    
    mov = zeros(800,800,1,length(rampNs),'uint8');
    
    for i=1:length(rampNs)
        hh = figure('position',[100,100,800,800]);
        thiscontour = numberlookup(zs,vzs,ones(size(zs))*rampNs(i));
        fprintf('Max at Ramp=%d is %d\n',rampNs(i),max(max(thiscontour)))
        [C, h] = contourf(zs,vzs,thiscontour,lines);
        xlabel('Final Phase (deg)','FontSize',12)
        ylabel('Final v_f (m/s)','FontSize',12)
        %zlabel(,'FontSize',12)
        title(['Optimize Loading, 300 Stages, Argon, RampN=' num2str(rampNs(i))],'FontSize',14)
        a = colorbar;
        ylabel(a,'Pop Remaining after 5ms Trapping','FontSize',12)
        set(a,'FontSize',12);
        set(gca,'FontSize',12);
        
        set(a,'Limits',[minr maxr]);
        set(gca,'CLim',[minr maxr]);
        pause(1)
        set(hh,'Position',[100,100,800,800])
        pause(1)
        
        frame = getframe(hh);
        im = frame2im(frame);
        if i==1
            [A,map] = rgb2ind(im,256);
        else
            A = rgb2ind(im,map);
        end

        mov(:,:,1,i) = A;
    end
    
    file = 'ramptestargon160_broad.gif';
    path = '~/Documents/MATLAB/slowANDtrap/Figures/';
    
    if exist([path file])
        file = input('File name exists. Enter a new one:');
    end
    
    imwrite(mov, map, [path file] , 'DelayTime', 1, 'LoopCount', inf);
    
end