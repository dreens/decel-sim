%%
% Make use of the Min Depth Finder to study some plots!
%

phis = [-35:5:85, 89, 90, 95];
c = @(p) p + (p<=0)*180;

if ~exist('modes','var') || ~isstruct(modes)
    modes = struct();
end
if ~isfield(modes,'s1')
    modes.s1 = struct();
    modes.s1.decels = {'longdecel','longdecel'};
    modes.s1.phifunc = @(p) [p-180, p-90, p];
end
if ~isfield(modes,'s3')
    modes.s3 = struct();
    modes.s3.decels = {'longdecel','longdecel'};
    modes.s3.phifunc = @(p) [p-540, p-90, p];
end
if ~isfield(modes,'sf')
    modes.sf = struct();
    modes.sf.decels = {'singlerod','longdecel'};
    modes.sf.phifunc = @(p) [p-180, -p, p];
end
if ~isfield(modes,'vsf')
    modes.vsf = struct();
    modes.vsf.decels = {'ppgg','longdecel'};
    modes.vsf.phifunc = @(p) [p-180, -p, p];
end
if ~isfield(modes,'xsf')
    modes.xsf = struct();
    modes.xsf.decels = {'ppmm_2mm','longdecel'};
    modes.xsf.phifunc = @(p) [p-180, -p, p];
end
%% Vary Phi2
for p2=10:10:170
    thisM = ['vsfe' num2str(p2)];
    if ~isfield(modes,thisM)
        modes.(thisM) = struct();
        modes.(thisM).decels = {'ppgg','longdecel'};
        modes.(thisM).phifunc = @(p) [p-180, p-p2, p];
    end
    thisM = ['xsfe' num2str(p2)];
    if ~isfield(modes,thisM)
        modes.(thisM) = struct();
        modes.(thisM).decels = {'ppmm_2mm','longdecel'};
        modes.(thisM).phifunc = @(p) [p-180, p-p2, p];
    end

end
%% All of the needed potentials and min depths.
allmodes = fields(modes);
for i=1:length(allmodes)
    fprintf('Mode %s:\n',allmodes{i})
    thisM = modes.(allmodes{i});
    if length(thisM) < 90
        thisM(2:90) = thisM(1);
    end
    if length(thisM) < 180
        thisM(91:180) = thisM(6);
    end
    for p=phis
        fprintf(' Phi=%d\n',p)
        if ~isfield(thisM(c(p)),'pot') || isempty(thisM(c(p)).pot)
            [thisM(c(p)).pot, thisM(c(p)).acc] = ...
                efftrap3Dgen(thisM(c(p)).decels,thisM(c(p)).phifunc(p),[1 1]);
        end
        if ~isfield(thisM(c(p)),'dep') || isempty(thisM(c(p)).dep)
            [thisM(c(p)).dep, thisM(c(p)).typ, thisM(c(p)).vol, thisM(c(p)).psv] = ...
                effTrapMinDepth(thisM(c(p)).pot);
        end
    end
    modes.(allmodes{i}) = thisM;
end

%%
%figure;
%phis = [1 5:5:20 70:5:85 89];
%plot(phis,vsfp(phis),'DisplayName','VSF''')
%hold on
%plot(phis,xsfp(phis),'DisplayName','XSF''')
%plot(phis,sf(phis),'DisplayName','SF')
%plot(phis,vsf(phis),'DisplayName','VSF')
%plot(phis,xsf(phis),'DisplayName','XSF')
%legend('show')

%% Plot Lots of Things, Edit as Needed
figure; hold on
allModes = fields(modes);
for i=1:length(allModes)
    modeName = allModes{i};
    if length(modeName)>=6 && strcmp(modeName(1),'v')
        thisM = modes.(modeName);
        accs = [thisM([180 1:90]).acc]/1e3;
        depths = [thisM([180 1:90]).dep]/1.38e-23/1e-3;
        plot(accs,depths,'DisplayName',allModes{i})
    end
end

xlabel('Deceleration (km/s/s)')
ylabel('Trap Depth (mK)')
title('Effective Trap Depths')

%% Plot the Non-Optimized, Well-Defined Modes
figure; hold on
allModes = fields(modes);
for i=1:length(allModes)
    modeName = allModes{i};
    if length(modeName)<6
        thisM = modes.(modeName);
        accs = [thisM([180 1:90]).acc]/1e3;
        depths = [thisM([180 1:90]).dep]/1.38e-23/1e-3;
        plot(accs,depths,'DisplayName',allModes{i})
    end
end

xlabel('Deceleration (km/s/s)')
ylabel('Trap Depth (mK)')
title('Effective Trap Depths')

%% Optimized envelope for XSF mod
% Let's just grab the peaks.
%figure;
allModes = fields(modes);
locs = [];
peaks = [];
for i=1:length(allModes)
    modeName = allModes{i};
    if length(modeName)>=6 && strcmp(modeName(1:2),'xs')
        thisM = modes.(modeName);
        depths = [thisM.dep];
        accs = [thisM.acc];
        deptweak = accs*1e-7*1.38e-23;
        [~, loc] = max(depths+deptweak);
        peaks(end+1) = depths(loc);
        locs(end+1) = accs(loc);
    end
end

% above a certain cutoff they stop working well:
locs = locs(2:11);
peaks = peaks(2:11);

% next we need to add the leading results onto the end:
locs  = [locs  [modes.xsfe110([70:5:85 89]).acc]];
peaks = [peaks [modes.xsfe110([70:5:85 89]).dep]];

% change units
locs = locs * 1e-3;
peaks = peaks / 1.38e-23 / 1e-3;

% close to zero phase angle, nothing seems to work perfectly. I'm just
% going to hard code it, these peaks come by plotting full envelopes.
locs = [1.46 15.02 21.49 locs];
peaks = [197.4 217.2 228.8 peaks];

plot(locs,peaks,'DisplayName','XSF*')

%% Optimized envelope for VSF mod
% Let's just grab the peaks.
%figure;
allModes = fields(modes);
locs = [];
peaks = [];
for i=1:length(allModes)
    modeName = allModes{i};
    if length(modeName)>=6 && strcmp(modeName(1:2),'vs')
        thisM = modes.(modeName);
        depths = [thisM.dep];
        accs = [thisM.acc];
        deptweak = accs*5e-7*1.38e-23;
        [~, loc] = max(depths+deptweak);
        peaks(end+1) = depths(loc);
        locs(end+1) = accs(loc);
    end
end

% above a certain cutoff they stop working well:
locs = locs([2:10 12]);
peaks = peaks([2:10 12]);

% next we need to add the leading results onto the end:
locs  = [locs  [modes.vsfe130([70:5:85 89]).acc]];
peaks = [peaks [modes.vsfe130([70:5:85 89]).dep]];

% change units
locs = locs * 1e-3;
peaks = peaks / 1.38e-23 / 1e-3;

% close to zero phase angle, nothing seems to work perfectly. I'm just
% going to hard code it, these peaks come by plotting full envelopes.
%locs = [1.46 15.02 21.49 locs];
%peaks = [197.4 217.2 228.8 peaks];

plot(locs,peaks,'DisplayName','VSF*')

%% Delete Phi2 Related Modes
if false
    allModes = fields(modes);

    for i=1:length(allModes)
        if length(allModes{i})==7
            modes = rmfield(modes,allModes{i});
        end
    end

end


%% Add phase information


