%% Plots results
function resultsgenerateloaddist(rs)
    loadt =  zeros(length(rs),1);
    loadvz = zeros(length(rs),1);
    number = zeros(length(rs),1);
    
    for i=1:length(rs)
        r = rs(i);
        loadt(i) =  r.loadtime;
        loadvz(i) = r.finalvz;
        number(i) = r.molnum(end);
    end
    
    loadts = unique(loadt);
    loadvzs = unique(loadvz);
    
    iff = @(varargin) varargin{varargin{3}+1};
    clean = @(x) iff(x,0,isempty(x));
    numberlookupsingle = @(t,v) clean(number(t==loadt & v==loadvz));
    numberlookup = @(t,v) arrayfun(numberlookupsingle,t,v);
    
    minr = min(number); maxr = max(number);
    
    steporder = round(log10(maxr-minr))-1;
    
    minr = round(minr,-steporder);
    maxr = round(maxr,-steporder);
    
    lines = minr:10^steporder:maxr;
    
    [ts, vzs] = meshgrid(loadts,loadvzs);
    
    
    hh = figure('position',[100,100,800,800]);
    thiscontour = numberlookup(ts,vzs);
    fprintf('Max is %d\n',max(max(thiscontour)))
    [C, h] = contourf(ts,vzs,thiscontour,lines);
    xlabel('Loading Time (s)','FontSize',12)
    ylabel('Final v_z (m/s)','FontSize',12)
    title('Loading Distribution Pool, 143 Stages','FontSize',14)
    a = colorbar;
    ylabel(a,'Pop Remaining after 5ms Trapping','FontSize',12)
    set(a,'FontSize',12);
    set(gca,'FontSize',12);

    set(a,'Limits',[minr maxr]);
    set(gca,'CLim',[minr maxr]);


end