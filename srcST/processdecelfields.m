function processdecelfields(r)
%% Load in from COMSOL
% Takes a COMSOL .dat file containing fields in a deceleration stage
%
% Expectations:
% * z is the decelerator axis, ranging from phase angle of $-90^\circ$ to $+90^\circ$.
% * x and y are symmetric.
% * First data column is a geometry mask, then E-field, then B and angle if relevant.
% * The Geometry mask is 1 where obstacles exist.
% * E-field in V/m, B-field in Tesla, angle in radians, xyz in mm.
% * It is assumed that odd numbered stages have field symmetry with respect
% to even numbered stages, but rotated.
% * Data points are given on a rectangular, uniform grid, although the
% spacing of the three dimensions need not be identical.

    fprintf('%s\n','Processing Decelerator Fields from COMSOL...');
    
    % COMSOL files usually have 9 header lines.
    data = importdata(['../Decels/' r.decel '.dat'],' ',9);
    
    % data is a struct with data and text header. We ignore the header and
    % take out the data.
    all = data.data;
    
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
    
    % not all decels will have b-field information.
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
    
    % z is assumed to run from $-90^\circ$ to $+90^\circ$, but no
    % assumptions are made about its actual coordinate range in COMSOL.
    % Thus it is shifted so that $+90^\circ$ phase is at zero.
    z = z - zs(end);
    zs = zs - zs(end);
    zstagel = -2*zs(1);
    
    % create lookup functions that tell you 'n' for a given x value such
    % that x is the nth x value.
    x2i = @(x) int16(1+x/xsp);
    y2i = @(y) int16(1+y/ysp);
    z2i = @(z) int16(length(zs)+z/zsp);
    
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
    
    % for each single-letter variable in varname, create a new variable
    % given by a double-letter, which is a 3D matrix of size fullsize, and
    % fill it by indexing into it with the locs linear index column.
    for varname={'b','e','t','m','x','y','z'}
        eval([varname{1} varname{1} '=zeros(fullsize);']);
        eval([varname{1} varname{1} '(locs) = ' varname{1} ';']);
    end
    mm = logical(mm);
    
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
    % vv2 = smooth3(vv,'gaussian',7,3);
    
    %% Bilateral Filtering to preserve fields near rods
    % It turns out that this will work better on the derivatives, since the
    % bilateral filter makes smoothing less effective on slopes.
    %{
    vv3 = vv; vv3(mm) = 5e-21;
    vv3 = shiftableBF3D(vv3,3,5e-23,0.0001,5e-21);
    
    %%
    
    figure('Position',[0 0 1500 800]);
    subplot(1,3,1)
    surf(cap(squeeze(vv(1,:,:)),1e-19))
    title('Plain')
    
    subplot(1,3,2)
    surf(cap(squeeze(vv2(1,:,:)),1e-19))
    title('Smoothed')
    
    subplot(1,3,3)
    surf(cap(squeeze(vv3(1,:,:)),1e-19))
    title('Bilateral Filtering')

    figure('Position',[0 0 1600 600])
    subplot(1,3,1)
    surf(cap(del2(squeeze(vv(1,:,:))),1e-25))
    title('Plain')
    
    subplot(1,3,2)
    surf(cap(del2(squeeze(vv2(1,:,:))),1e-25))
    title('Plain')
    
    subplot(1,3,3)
    surf(cap(del2(squeeze(vv3(1,:,:))),1e-25))
    title('Bilateral Filtering')
    
    vv1 = vv;
    vv = vv3;
    %}
    
   
    %% Take derivatives to get force fields
    
    % get symmetric derivative kernels so we can differentiate the
    % potential matrix via convolution
    xd = reshape([-.5 0 .5],3,1,1);
    yd = reshape(xd,1,3,1);
    zd = reshape(xd,1,1,3);
    
    % perform the derivatives. Convolution style differentiation requires
    % scaling by the matrix point spacing.
    dvdxu = convn(vv,xd,'same')/xsp;
    dvdyu = convn(vv,yd,'same')/ysp;
    dvdzu = convn(vv,zd,'same')/zsp;
    
    % zero the forces outside of the geometry mask. This ensures that
    % molecules aren't accidentally reflected off of pathologically large
    % field spikes that can occur near conductor or magnet surfaces in
    % COMSOL. Instead, molecules that hit geometry will continue through
    % until they fall out of the force field and are removed by the
    % molecule stepper which looks for this.
    mmx = convn(mm,abs(xd),'same') > 0;
    mmy = convn(mm,abs(yd),'same') > 0;
    mmz = convn(mm,abs(zd),'same') > 0;
    dvdxu(mmx)=0;
    dvdyu(mmy)=0;
    dvdzu(mmz)=0;
    dvdxu(1,:,:)=0; dvdxu(end,:,:) = dvdxu(end-1,:,:);
    dvdyu(:,1,:)=0; dvdyu(:,end,:) = dvdyu(:,end-1,:);
    dvdzu(:,:,[1 end])=0;    
    dvdxu(mmx)=2e-19;
    dvdyu(mmy)=2e-19;
    dvdzu(mmz)=2e-19;
    dvdxm = shiftableBF3D(dvdxu,2,2e-20,1e-4,2e-19);
    dvdym = shiftableBF3D(dvdyu,2,2e-20,1e-4,2e-19);
    dvdzm = shiftableBF3D(dvdzu,2,2e-20,1e-4,2e-19);
%     dvdxm = dvdxu;
%     dvdym = dvdyu;
%     dvdzm = dvdzu;
    dvdxm(mmx)=0;
    dvdym(mmy)=0;
    dvdzm(mmz)=0;

    
    % zero the x,y force along their respective lines of symmetry. This
    % should already be the case but convolution based derivatives can
    % behave strangely near borders due to zero padding assumptions.
    dvdxm(1,:,:)=0;
    dvdym(:,1,:)=0;
    dvdzm(:,:,[1 end])=0;
    
    %% Create interpolants
    % These convenient datatypes can be evaluated directly as functions and
    % carry with them all of the gridded data used to instantiate them.
    dvdxg  = griddedInterpolant(xx,yy,zz,dvdxm,'linear','none');
    dvdyg  = griddedInterpolant(xx,yy,zz,dvdym,'linear','none');
    dvdzg  = griddedInterpolant(xx,yy,zz,dvdzm,'linear','none');
    vfg = griddedInterpolant(xx,yy,zz,vv,'linear','none');
    % bf = griddedInterpolant(xx,yy,zz,bb);
    % ef = griddedInterpolant(xx,yy,zz,ee);
    % mf = griddedInterpolant(xx,yy,zz,mm);
    % tf = griddedInterpolant(xx,yy,zz,tt);
    
    %% Create helpful lookup functions
    % Anonymous functions created here will be saved along with the
    % relevant workspace when defined. This is an automatic matlab feature.
    % This allows the intricacies of stage parity and reflection symmetry
    % to be encapsulated here and not dealt with in the functions that use
    % these lookup tables.
    
    
    % This phase lookup function gives the phase angle as a function of the
    % z coordinate and the parity of the stage number n. It gives the phase
    % in the range -270 to 90, which is convenient since the decelerator
    % will be run with at a phase angle close to +90 degrees, and thus
    % during a single stage most well-decelerated molecules won't wrap
    % their phase angles as returned by this lookup.
    phase = @(z,n) mod(z/zstagel+n/2,1)*360 - 270;
    
    % wrap returns a z-coordinate within the force lookup table given a
    % general z-coordinate and the stage parity. It achieves this in two
    % steps. First z is converted to a coordinate between -zstagel and 0
    % which is also the range -270 to +90 in phase angle. Then, the mirror
    % symmetry between the forces in -270--90 and -90-+90 is exploited. The
    % side lookup indicates whether this symmetry is exploited so the
    % z-forces can be inverted, since molecules in the -270 to -90 range
    % are accelerated, not decelerated.
    wrapc = @(z,n) (phase(z,n)-90)/360 * zstagel;
    wrap = @(z,n) abs(wrapc(z,n)+zstagel/2)-zstagel/2;
    side = @(z,n) (wrapc(z,n) > -zstagel/2)*2 - 1;
    
    % This returns the energy removed per stage as a function of phase
    % angle. Its inverse enables quickly choosing the phase angle given a
    % final velocity. 
    renergy = @(phi) vfg(0,0,(phi-90)/360 * zstagel) - ...
        vfg(0,0,(-phi-90)/360 * zstagel);
    
    % This is the potential energy at a given phase angle, measured
    % relative to the potential energy at phi=-90 degrees. One could
    % subtract this from itself reversed to get renergy above.
    aenergy = @(phi) vfg(0,0,(phi-90)/360 * zstagel) - vfg(0,0,-zstagel/2);

    % These functions reference the gridded interpolants, but with the
    % coordinates appropriately wrapped.
    dvdxa = @(x,y,z,n) dvdxg(abs(x),abs(y),wrap(z,n)).*sign(x);
    dvdya = @(x,y,z,n) dvdyg(abs(x),abs(y),wrap(z,n)).*sign(y);
    dvdza = @(x,y,z,n) dvdzg(abs(x),abs(y),wrap(z,n)).*side(z,n);
    vfa = @(x,y,z,n) vfg(abs(x),abs(y),wrap(z,n));
    dvdx = @(xyz,n) dvdxa(xyz(:,1),xyz(:,2),xyz(:,3),n);
    dvdy = @(xyz,n) dvdya(xyz(:,1),xyz(:,2),xyz(:,3),n);
    dvdz = @(xyz,n) dvdza(xyz(:,1),xyz(:,2),xyz(:,3),n);
    vf = @(xyz,n) vfa(xyz(:,1),xyz(:,2),xyz(:,3),n);
    
    % Save in a file for loading and propagating during decel simulation.
    save(['../Decels/' r.decel '.mat'],'dvdx','dvdy','dvdz',...
        'vf','phase','zstagel','renergy','aenergy');

    %% Produce Output Figure for Debugging
    % There are many potential errors that could be made in the COMSOL
    % output, data input processing. 
    cc = 4e-20;
    
    figure('position',[50,50,1100,1100])
    subplot(2,2,1)
    surf(cap(squeeze(dvdzm(:,1,:)),cc));
    title(['Decelerator dvdz, X-Z plane, ' r.decel '.dat']);

    % You might think the labels are backwards, but they're not. surf uses
    % the second index (the column of the matrix) as the x-axis and the
    % first index (the row) as the y-axis. It makes sense if you think of
    % matrices as oriented with x going left-right and y going up-down, but
    % its crazy when working in 3D.
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    subplot(2,2,2)
    surf(cap(squeeze(dvdzm(1,:,:)),cc));
    title(['Decelerator dvdz, Y-Z plane, ' r.decel '.dat']); 
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    subplot(2,2,3)
    surf(cap(squeeze(dvdxm(:,1,:)),cc));
    title(['Decelerator dvdx, X-Z plane, ' r.decel '.dat']);
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    subplot(2,2,4)
    surf(cap(squeeze(dvdym(1,:,:)),cc));
    title(['Decelerator dvdy, Y-Z plane, ' r.decel '.dat']);
    xlabel(['Z axis (' num2str(zsp) ')']);
    ylabel(['X axis (' num2str(xsp) ')']);
    zlim([-cc cc])

    %{
    figure('position',[100,200,400,400])
    plot(abs(zs),vf([zeros(length(zs),2) zs(:)],1))
    title(['Decel Potential along Z-axis, ' r.trapname '.dat'])
    xlabel('Z axis')
    ylabel('Potential Energy (J)')
    %}
end

function x = cap(x,varargin)
    if length(varargin)==1
        c = varargin{1};
        x(x>c) = c;
        x(x<-c) = -c;
    else
        [cl, ch] = varargin{:};
        x(x>ch) = ch;
        x(x<cl) = cl;
    end
end
