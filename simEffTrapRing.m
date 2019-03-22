%% Simulate Effective Trap
% Here we simulate the effective trap by randomly initializing a
% homogeneous phase space distribution and evolving it for a few
% milliseconds.
function varargout = simEffTrapRing(varargin)
    assert(~isempty(varargin),'simEffTrap Requires an input potential.');
    pot = varargin{1};
    assert(~isempty(pot),'Can''t simulate an empty potential.');

    settings = struct();
    settings.maxVel = 10;
    settings.time = 3e-3;
    settings.step = 1e-6;
    settings.num = 1e3;
    settings.zwidth = 3; %half width of the initialized volume in z direction
    
    lv = length(varargin);
    assert(~~mod(lv,2),'There must be an odd number of inputs, first the potential, then name value pairs.');
    if lv>1
        properties = varargin(2:2:end);
        values = varargin(3:2:end);
        for i=1:length(properties)
            assert(ischar(properties{i}),['Input #' num2str(i*2+1) ' must be the name of a property, but is not even a string.'])
            assert(isfield(settings,properties{i}),['Property ' properties{i} ' is not a valid setting for simEffTrap.'])
            settings.(properties{i}) = values{i};
        end
    end
    
    mOH = 2.82328e-26; % Accounts for Oxygen binding energy
    
    pos = rand(settings.num,2)*2;
    pos(:,2) = (pos(:,2)-1)*settings.zwidth;
    pos = pos*1e-3; % convert to mm
    vel = rand(settings.num,2)*2-1;
    vel = vel * settings.maxVel;
    
    
    
    % get symmetric derivative kernels so we can differentiate the
    % potential matrix via convolution
    rd = [-.5; 0; .5];
    zd = rd';
    
    % perform the derivatives. Convolution style differentiation requires
    % scaling by the matrix point spacing.
    sp = 100e-6;
    dvdr = convn(pot,rd,'same')/sp;
    dvdz = convn(pot,zd,'same')/sp;
    
    % convolutions behave funny on the borders, zero these.
    dvdr([1 end],:) = 0;
    dvdz(:,[1 end]) = 0;
    
    tv = 0:sp:2e-3;
    lg = -3e-3:sp:3e-3;
    [rr,zz] = ndgrid(tv,lg);
    ar  = griddedInterpolant(rr,zz,dvdr/mOH,'linear','none');
    az  = griddedInterpolant(rr,zz,dvdz/mOH,'linear','none');    
    
    for time=0:settings.step:settings.time
        pos = pos + vel*settings.step/2;
        vel = vel + [ax(pos(:,1),pos(:,2),pos(:,3)), ay(pos(:,1),pos(:,2),pos(:,3)), az(pos(:,1),pos(:,2),pos(:,3))]*settings.step;
        pos = pos + vel*settings.step/2;
        if ~mod(time/settings.step,100) || time == settings.time
            lost = any(isnan(pos),2);
            pos = pos(~lost,:);
            vel = vel(~lost,:);
        end
    end
    
    num = size(pos,1);
    temp = mean(.5*mOH*sum(vel.^2,2))/1.38e-23;
    initPSV = (2*settings.maxVel)^3 * (2e-3 * settings.zwidth)^3;
    psv = num/settings.num * initPSV;
    varargout = {num,psv,temp};
    varargout = varargout(1:max(nargout,1));
end











