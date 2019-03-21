%%
% Make use of the Min Depth Finder to study some plots!
%

phis = [-20:5:100];
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
if ~isfield(modes,'vsf0')
    modes.vsf0 = struct();
    modes.vsf0.decels = {'ppgg','None'};
    modes.vsf0.phifunc = @(p) [p-180, p-90, p];
end
if ~isfield(modes,'xsf0')
    modes.xsf0 = struct();
    modes.xsf0.decels = {'ppmm_2mm','None'};
    modes.xsf0.phifunc = @(p) [p-180, p-90, p];
end
if ~isfield(modes,'vsf01')
    modes.vsf01 = struct();
    modes.vsf01.decels = {'ppgg','None'};
    modes.vsf01.phifunc = @(p) [p-190, p-80, p-10];
end
if ~isfield(modes,'vsf02')
    modes.vsf02 = struct();
    modes.vsf02.decels = {'ppgg','None'};
    modes.vsf02.phifunc = @(p) [p-200, p-70, p-20];
end
if ~isfield(modes,'vsf03')
    modes.vsf03 = struct();
    modes.vsf03.decels = {'ppgg','None'};
    modes.vsf03.phifunc = @(p) [p-210, p-60, p-30];
end
if ~isfield(modes,'xsf01')
    modes.xsf01 = struct();
    modes.xsf01.decels = {'ppmm_2mm','None'};
    modes.xsf01.phifunc = @(p) [p-190, p-80, p-10];
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
if false
figure; hold on
allModes = fields(modes);
for i=1:length(allModes)
    modeName = allModes{i};
    if length(modeName)>=5 && strcmp(modeName(1),'x')
        thisM = modes.(modeName);
        accs = [thisM([105:180 1:100]).acc]/1e3;
        depths = [thisM([105:180 1:100]).dep]/1.38e-23/1e-3;
        plot(accs,depths,'DisplayName',allModes{i})
    end
end

xlabel('Deceleration (km/s/s)')
ylabel('Trap Depth (mK)')
title('Effective Trap Depths')
end
%% Plot the Non-Optimized, Well-Defined Modes
figure; hold on
allModes = fields(modes);
plotModes = {'s1','s3','sf'};
for i=1:length(allModes)
    modeName = allModes{i};
    if any(strcmp(modeName,plotModes))
        thisM = modes.(modeName);
        accs = [thisM([180 1:90]).acc]/1e3;
        depths = [thisM([180 1:90]).dep]/1.38e-23/1e-3;
        plot(accs,depths,'DisplayName',allModes{i},'LineWidth',2)
    end
end

xlabel('Deceleration (km/s/s)')
ylabel('Trap Depth (mK)')
title('Effective Trap Depths')
xlim([0 400])

%% Optimized envelope for VSF mod
% Let's just grab the peaks.
%figure;
allModes = fields(modes);
bestModesVSF = {};
bestPhisVSF = [];
for i=1:length(allModes)
    modeName = allModes{i};
    if length(modeName)>=5 && strcmp(modeName(1:2),'vs')
        thisM = modes.(modeName);
        phidef = [1:90 -89:0];
        thesePhis = phidef(~cellfun(@isempty,{thisM.dep}));
        depths = [thisM.dep];
        accs = [thisM.acc];
        deptweak = accs*5e-7*1.38e-23;
        [~, loc] = max(depths+deptweak);
        bestPhisVSF(end+1) = thesePhis(loc);
        bestModesVSF{end+1} = modeName;
    end
end

% above a certain cutoff they stop working well:
bestPhisVSF = bestPhisVSF([2:10 12]);
bestModesVSF = bestModesVSF([2:10 12]);

% next we need to add the leading results onto the end:
bestPhisVSF = [bestPhisVSF 70:5:90];
bestModesVSF(end+1:end+5) = {'vsfe130'};

% close to zero phase angle, nothing seems to work perfectly. I'm just
% going to hard code it, these peaks come by plotting full envelopes.
% Actually, let's grab peaks from the VSF modification where no
% conventional charge config is used (vsf0, vsf01, vsf02 above)
bestPhisVSF = [45 -5 bestPhisVSF];
bestModesVSF = {'vsf03' 'vsfe20' bestModesVSF{:}};

bestAccs = bestPhisVSF;
bestPeaks = bestPhisVSF;
for ii = 1:length(bestAccs)
    thisM = modes.(bestModesVSF{ii});
    bestAccs(ii) = thisM(c(bestPhisVSF(ii))).acc;
    bestPeaks(ii) = thisM(c(bestPhisVSF(ii))).dep;
end

% change units
bestAccs = bestAccs * 1e-3;
bestPeaks = bestPeaks / 1.38e-23 / 1e-3;

plot(bestAccs,bestPeaks,'DisplayName','vsf*','LineWidth',2)

%% Optimized envelope for XSF mod
% Let's just grab the peaks.
%figure;
allModes = fields(modes);
bestModesXSF = {};
bestPhisXSF = [];
for i=1:length(allModes)
    modeName = allModes{i};
    if length(modeName)>=5 && strcmp(modeName(1:2),'xs')
        thisM = modes.(modeName);
        phidef = [1:90 -89:0];
        thesePhis = phidef(~cellfun(@isempty,{thisM.dep}));
        depths = [thisM.dep];
        accs = [thisM.acc];
        deptweak = accs*5e-7*1.38e-23;
        [~, loc] = max(depths+deptweak);
        bestPhisXSF(end+1) = thesePhis(loc);
        bestModesXSF{end+1} = modeName;
    end
end

% above a certain cutoff they stop working well:
bestPhisXSF = bestPhisXSF([2:11]);
bestModesXSF = bestModesXSF([2:11]);

% next we need to add the leading results onto the end:
bestPhisXSF = [bestPhisXSF 70:5:90];
bestModesXSF(end+1:end+5) = {'xsfe110'};

% close to zero phase angle, nothing seems to work perfectly. I'm just
% going to hard code it, these peaks come by plotting full envelopes.
% Actually, let's grab peaks from the VSF modification where no
% conventional charge config is used (vsf0, vsf01, vsf02 above)
bestPhisXSF = [0 0 0 bestPhisXSF];
bestModesXSF = {'xsfe50' 'xsfe40' 'xsfe30' bestModesXSF{:}};

bestAccs = bestPhisXSF;
bestPeaks = bestPhisXSF;
for ii = 1:length(bestAccs)
    thisM = modes.(bestModesXSF{ii});
    bestAccs(ii) = thisM(c(bestPhisXSF(ii))).acc;
    bestPeaks(ii) = thisM(c(bestPhisXSF(ii))).dep;
end

% change units
bestAccs = bestAccs * 1e-3;
bestPeaks = bestPeaks / 1.38e-23 / 1e-3;

plot(bestAccs,bestPeaks,'DisplayName','xsf*','LineWidth',2)

%% Add in T-Wave Decel
ringprocess
xlabel('Deceleration (km/s/s)','FontSize',13)
ylabel('Worst Case Trap Depth (mK)','FontSize',13)
title('Trap Depth v Deceleration','FontSize',14)
set(gca,'FontSize',13)
hl = legend('show');
hl.FontSize = 13;

%% Delete Phi2 Related Modes
if false
    allModes = fields(modes);

    for i=1:length(allModes)
        if length(allModes{i})==7
            modes = rmfield(modes,allModes{i});
        end
    end

end

%% Actual Simulation?
% Here we make an attempt to actually load up these traps and see what
% happens!
%
% Plan is to load the traps with a randomly initialized but homogeneous
% phase space density that more than overlaps the trap, let them evolve for
% some typical time like 3 milliseconds, and then infer the trapped phase
% space volume from the remaining number.
%
% Did this mostly in simEffTrap. Here we just loop through and plot.
figure; hold on
allModes = fields(modes);
plotModes = {'s1','s3','sf'};
for i=1:length(allModes)
    modeName = allModes{i};
    if any(strcmp(modeName,plotModes))
        thisM = modes.(modeName);
        accs = [thisM([180 c(phis(2:end))]).acc]/1e3;
        psvs = [];
        for p=phis
            fprintf(' Phi=%d\n',p)
            if isfield(thisM(c(p)),'pot') && ~isempty(thisM(c(p)).pot)
                [num, psv] = simEffTrap(thisM(c(p)).pot,'num',1e2);
                psvs = [psvs psv];
            end
        end
 
        
        plot(accs,psvs,'DisplayName',allModes{i},'LineWidth',2)
    end
end

xlabel('Deceleration (km/s/s)')
ylabel('Phase Space Volume (m^6/s^3)')
title('Effective Trap Depths')
xlim([0 400])


%% Next add in the VSF and XSF stuff.
accs = bestPhisVSF;
psvs = bestPhisVSF;
for i=1:length(bestPhisVSF)
    thisM = modes.(bestModesVSF{i});
    thisM = thisM(c(bestPhisVSF(i)));
    accs(i) = thisM.acc;
    [num, psv] = simEffTrap(thisM.pot,'num',1e2);
    psvs(i) = psv;
end
plot(accs,psvs,'DisplayName','vsf*','LineWidth',2)

accs = bestPhisXSF;
psvs = bestPhisXSF;
for i=1:length(bestPhisXSF)
    thisM = modes.(bestModesXSF{i});
    thisM = thisM(c(bestPhisXSF(i)));
    accs(i) = thisM.acc;
    [num, psv] = simEffTrap(thisM.pot,'num',1e2);
    psvs(i) = psv;
end
plot(accs,psvs,'DisplayName','xsf*','LineWidth',2)
    
    
    
    
    
    
