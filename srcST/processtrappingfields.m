function processtrappingfields(r)
%% Load trapping fields from COMSOL
% Expectations:
% * z is still decel axis.
% * z=0 is the trap center.
% * Column order, data point gridding as above.
    
    fprintf('%s\n','Processing Trapping Fields from COMSOL...');
    
    % COMSOL files usually have 9 header lines.
    data = importdata(['../Traps/' r.trapname '.dat'],' ',9);
    
    % data is a struct with data and text header. We ignore the header and
    % take out the data.
    all = data.data;
    
    % end of the cell
    %% Convert data to 3D format.
    %{
    % After COMSOL export, the data are in list form. Each data row gives
    % the mask, efield, bfield, and angle evaluated at one point given by
    % x,y,z in the first three columns.
    %
    % To get the data into 3D matrices, we first infer the 3D span of the
    % datapoints by checking for unique values in x,y,z. Then we get the
    % linear index of the 3D matrix corresponding to each row in the COMSOL
    % list. Finally we can fill the 3D matrix just be indexing into the 3D
    % matrix with this linear index.
    %}
    
    % read x,y,z,m,e,b,t from each data column. 
    x = all(:,1)*1e-3; %convert to meters.
    y = all(:,2)*1e-3;
    z = all(:,3)*1e-3;
    m = all(:,4);
    e = all(:,5);   %units of V/m
    
    % processtrappingfields can be called for processing in a detection
    % situation. We ascertain whether this is the case by checking the
    % simulation timing.
    if r.time > r.loadingofftime
        if strcmp(r.detection,'free flight')
            e = e*0;
        end
    end
    
    % I made an error generating the magnetic pin data from COMSOL. The
    % trap center should have been at z=0, but instead it was one rod pair
    % away.
    if strcmp(r.trapname,'magneticpintrap')
        z = z - 5.461e-3;
    end
    
    % not all traps will have b-field information.
    if size(all,2)>5
        b = all(:,6);   %units of T
        t = all(:,7);   %radians
        t(isnan(t))=0;
        t(imag(t)~=0)=0;
    else
        b = zeros(size(e));
        t = b;
    end
    
    % find unique x,y,z coordinates
    xs = sort(uniquetol(x,1e-6,'DataScale',1));
    ys = sort(uniquetol(y,1e-6,'DataScale',1));
    zs = sort(uniquetol(z,1e-6,'DataScale',1));
    
    % get the datapoint spacing, used for taking derivatives later.
    xsp = uniquetol(diff(xs));
    ysp = uniquetol(diff(ys));
    zsp = uniquetol(diff(zs));
    
    % It is assumed that the z-coordinate of the trapping fields is the
    % most sensible z-coordinate for the simulation to use. In order to
    % translate between this most useful coordinate definition and the one
    % used previously during deceleration and loading, we make sure to save
    % the maximum zvalue, since both loading and trapping are guaranteed to
    % have the same endpoint.
    zmaxtrap = max(zs);
    
    % create lookup functions that tell you 'n' for a given x value such
    % that x is the nth x value.
    x2i = @(x) int16(1+x/xsp);
    y2i = @(y) int16(1+y/ysp);
    z2i = @(z) int16(1+(z-zs(1))/zsp);
    
    % check for datapoint uniformity, and x,y starting from zero
    if length(xsp)+length(ysp)+length(zsp) > 3
        error('Non-uniform Datapoint Spacing');
    elseif abs(xs(1)) > 1e-6
        error('X data doesn''t begin at zero')
    elseif abs(ys(1)) > 1e-6
        error('Y data doesn''t begin at zero')
    end
    
    % the size of the 3D data matrices to be filled
    fullsize = [length(xs) length(ys) length(zs)];
    
    % for each x,y,z value in the COMSOL loaded x,y,z columns, check which
    % linear index this corresponds to in a 3D matrix of size fullsize.
    locs = sub2ind(fullsize,x2i(x),y2i(y),z2i(z));
    
    %{
    % for each single-letter variable in varname, create a new variable
    % given by a double-letter, which is a 3D matrix of size fullsize, and
    % fill it by indexing into it with the locs linear index column.
    %}
    for varname={'b','e','t','m','x','y','z'}
        eval([varname{1} varname{1} '=zeros(fullsize);']);
        eval([varname{1} varname{1} '(locs) = ' varname{1} ';']);
    end
    
    %% Calculate Stark-Zeemen Potential Energy
    
    % create an anonymous function that can give OH doubly stretched state
    % potential energy as a function of bfield, efield, and ebangle (t).
    last = @(x) x(end);
    energysingle = @(b,e,t) last(sort(eig(OH_Ham_Simple_SI(b,e,t))));
    energy = @(b,e,t)  arrayfun(energysingle,b,e,t);
    
    
    % get the potential energy
    vv = energy(bb,ee,tt);
    
    % smooth it for better derivatives. I tuned the standard deviation, 3,
    % until surfaces (run "figure;surf(squeeze(dvdzm(1,:,:)))" for example)
    % don't show waviness or chopiness.
    vv = smooth3(vv,'gaussian',7,3);
    
    %% Take derivatives to get force fields
    
    % get symmetric derivative kernels so we can differentiate the
    % potential matrix via convolution
    xd = reshape([-.5 0 .5],3,1,1);
    yd = reshape(xd,1,3,1);
    zd = reshape(xd,1,1,3);
    
    % perform the derivatives. Convolution style differentiation requires
    % scaling by the matrix point spacing.
    dvdxm = convn(vv,xd,'same')/xsp;
    dvdym = convn(vv,yd,'same')/ysp;
    dvdzm = convn(vv,zd,'same')/zsp;
    
    % zero the forces outside of the geometry mask. This ensures that
    % molecules aren't accidentally reflected off of pathologically large
    % field spikes that can occur near conductor or magnet surfaces in
    % COMSOL. Instead, molecules that hit geometry will continue through
    % until they fall out of the force field and are removed by the
    % molecule stepper which looks for this.
    % 
    % In tricycletrap the mask is set incorrectly, so that it is everywhere
    % 1 instead of 0. The correct convention is that the mask is nonzero
    % wherever obstacles exist. The condition accounts for this issue in
    % that particular case.
    %
    if any(any(any(mm==0)))
        dvdxm(mm~=0)=0;
        dvdym(mm~=0)=0;
        dvdzm(mm~=0)=0;
    end
    
    % zero the x,y force along their respective lines of symmetry. This
    % should already be the case but convolution based derivatives can
    % behave strangely near borders due to zero padding assumptions.
    dvdxm(1,:,:)=0; dvdxm(end,:,:)=0;
    dvdym(:,1,:)=0; dvdym(:,end,:)=0;
    dvdzm(:,:,1)=0; dvdzm(:,:,end)=0;
    
    %% Create interpolants
    % These convenient datatypes can be evaluated directly as functions and
    % carry with them all of the gridded data used to instantiate them.
    dvdxg  = griddedInterpolant(xx,yy,zz,dvdxm,'linear','none');
    dvdyg  = griddedInterpolant(xx,yy,zz,dvdym,'linear','none');
    dvdzg  = griddedInterpolant(xx,yy,zz,dvdzm,'linear','none');
    vfg = griddedInterpolant(xx,yy,zz,vv,'linear','none');
    bfg = griddedInterpolant(xx,yy,zz,bb);
    efg = griddedInterpolant(xx,yy,zz,ee);
    mfg = griddedInterpolant(xx,yy,zz,mm);
    tfg = griddedInterpolant(xx,yy,zz,tt);
    
    %% Create actual lookup functions
    % These functions reference the gridded interpolants, but with the
    % coordinates appropriately wrapped.
    dvdxa = @(x,y,z) dvdxg(abs(x),abs(y),z).*sign(x);
    dvdya = @(x,y,z) dvdyg(abs(x),abs(y),z).*sign(y);
    dvdza = @(x,y,z) dvdzg(abs(x),abs(y),z);
    dvdx = @(xyz,n) dvdxa(xyz(:,1),xyz(:,2),xyz(:,3));
    dvdy = @(xyz,n) dvdya(xyz(:,1),xyz(:,2),xyz(:,3));
    dvdz = @(xyz,n) dvdza(xyz(:,1),xyz(:,2),xyz(:,3));
    vf = @(xyz,n) vfg(abs(xyz(:,1)),abs(xyz(:,2)),xyz(:,3));
    bf = @(xyz) bfg(abs(xyz(:,1)),abs(xyz(:,2)),xyz(:,3));
    ef = @(xyz) efg(abs(xyz(:,1)),abs(xyz(:,2)),xyz(:,3));
    mf = @(xyz) mfg(abs(xyz(:,1)),abs(xyz(:,2)),xyz(:,3));
    tf = @(xyz) tfg(abs(xyz(:,1)),abs(xyz(:,2)),xyz(:,3));
    
    
    %% Save in a file for loading and propagating during decel simulation.
    sname = ['../Traps/' r.trapname '.mat'];
    if r.time > r.loadingofftime
        if strcmp(r.detection,'free flight')
            sname = ['../Traps/' r.trapname '_free.mat'];
        end
    end
    save(sname,'dvdx','dvdy','dvdz','vf','bf','ef','mf','tf','zmaxtrap');
    
    %% Produce Output Figure for Debugging
    % There are many potential errors that could be made in the COMSOL
    % output, data input processing. 
    figure('position',[50,400,1100,400])
    subplot(1,2,1)
    surf(squeeze(dvdzm(:,1,:)));
    title(['Trapping dvdz, X-Z plane, ' r.trapname '.dat']);

    % You might think the labels are backwards, but they're not. surf uses
    % the second index (the column of the matrix) as the x-axis and the
    % first index (the row) as the y-axis. It makes sense if you think of
    % matrices as oriented with x going left-right and y going up-down, but
    % its crazy when working in 3D.
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);

    subplot(1,2,2)
    surf(squeeze(dvdzm(1,:,:)));
    title(['Trapping dvdz, Y-Z plane, ' r.trapname '.dat']);
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    
    figure('position',[100,200,400,400])
    plot(zs,vf([zeros(length(zs),2) zs(:)],1))
    title(['Trapping Potential along Z-axis, ' r.trapname '.dat'])
    xlabel('Z axis')
    ylabel('Potential Energy (J)')
end
