%%
% Make use of the Min Depth Finder to study some plots!
%

phis = [0:5:35];
c = @(p) p + (p<=0)*180 + (p>90)*90;

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
% if ~isfield(modes,'xsfxb')
%     modes.xsfxb = struct();
%     modes.xsfxb.decels = {'longdecel','ppmm_2mm','longdecel'};
%     modes.xsfxb.phifunc = @(p) [p-180, p-135, p-45, p];
% end
% if ~isfield(modes,'xsfxb10')
%     modes.xsfxb10 = struct();
%     modes.xsfxb10.decels = {'longdecel','ppmm_2mm','longdecel'};
%     modes.xsfxb10.phifunc = @(p) [p-180, p-145, p-35, p];
% end
% if ~isfield(modes,'xsfxb20')
%     modes.xsfxb20 = struct();
%     modes.xsfxb20.decels = {'longdecel','ppmm_2mm','longdecel'};
%     modes.xsfxb20.phifunc = @(p) [p-180, p-155, p-25, p];
% end
if ~isfield(modes,'xsfxb30')
    modes.xsfxb30 = struct();
    modes.xsfxb30.decels = {'longdecel','ppmm_2mm','longdecel'};
    modes.xsfxb30.phifunc = @(p) [p-180, p-165, p-15, p];
end
% if ~isfield(modes,'xsfxbm10')
%     modes.xsfxbm10 = struct();
%     modes.xsfxbm10.decels = {'longdecel','ppmm_2mm','longdecel'};
%     modes.xsfxbm10.phifunc = @(p) [p-180, p-125, p-55, p];
% end

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
    if length(thisM) < 270
        thisM(181:270) = thisM(6);
    end
    for p=phis
        fprintf(' Phi=%d\n',p)
        if ~isfield(thisM(c(p)),'pot') || isempty(thisM(c(p)).pot)
            phisies = thisM(c(p)).phifunc(p);
            onesies = phisies(2:end)*0+1;
            [thisM(c(p)).pot, thisM(c(p)).acc] = ...
                efftrap3Dgen(thisM(c(p)).decels,phisies,onesies);
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
figure; subplot(2,1,1); hold on
allModes = fields(modes);
plotModes = {'s1','s3','sf'};
for i=1:length(plotModes)
    modeName = plotModes{i};
    thisM = modes.(modeName);
    accs = [thisM([180 1:90]).acc]/1e3;
    depths = [thisM([180 1:90]).dep]/1.38e-23/1e-3;
    if strcmp(modeName,'s3')
        accs(4) = [];
        depths(4) = [];
    end
    plot(accs,depths,'DisplayName',upper(modeName),'LineWidth',2)
end

xlabel('Deceleration (km/s/s)')
ylabel('Trap Depth (mK)')
title('Effective Trap Depths')
xlim([0 400])

%% Optimized envelope for VSF mod
% Let's just grab the peaks.
%figure;
grabNums = [20:10:100 120];
bestModesVSF = {};
for j=grabNums
    bestModesVSF{end+1} = ['vsfe' num2str(j)];
end
bestPhisVSF = grabNums*0;
for i=1:length(bestModesVSF)
    modeName = bestModesVSF{i};
    thisM = modes.(modeName);
    phidef = [1:90 -89:0];
    thesePhis = phidef(~cellfun(@isempty,{thisM.dep}));
    depths = [thisM.dep];
    accs = [thisM.acc];
    deptweak = accs*5e-7*1.38e-23;
    [~, loc] = max(depths+deptweak);
    bestPhisVSF(i) = thesePhis(loc);
end

% next we need to add the leading results onto the end:
bestPhisVSF = [bestPhisVSF 70:5:90];
bestModesVSF(end+1:end+5) = {'vsfe130'};

% close to zero phase angle, nothing seems to work perfectly. I'm just
% going to hard code it, these peaks come by plotting full envelopes.
% Actually, let's grab peaks from the VSF modification where no
% conventional charge config is used (vsf0, vsf01, vsf02 above)
bestPhisVSF = [45 -5 bestPhisVSF];
bestModesVSF = [{'vsf03' 'vsfe20'} bestModesVSF];

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

plot(bestAccs,bestPeaks,'DisplayName','VSF*','LineWidth',2)

%% Optimized envelope for XSF mod
% Let's just grab the peaks.
%figure;
grabModes = cell(1,length(20:10:110));
for j=20:10:110
    grabModes{j/10-1} = ['xsfe' num2str(j)];
end
bestModesXSF = grabModes;
bestPhisXSF = zeros(1,length(grabModes));
for i=1:length(grabModes)
    modeName = grabModes{i};
    thisM = modes.(modeName);
    phidef = [1:90 -89:0];
    thesePhis = phidef(~cellfun(@isempty,{thisM.dep}));
    depths = [thisM.dep];
    accs = [thisM.acc];
    deptweak = accs*5e-7*1.38e-23;
    [~, loc] = max(depths+deptweak);
    bestPhisXSF(i) = thesePhis(loc);
end

% next we need to add the leading results onto the end:
bestPhisXSF = [bestPhisXSF 70:5:90];
bestModesXSF(end+1:end+5) = {'xsfe110'};

% close to zero phase angle, nothing seems to work perfectly. I'm just
% going to hard code it, these peaks come by plotting full envelopes.
% Actually, let's grab peaks from the VSF modification where no
% conventional charge config is used (vsf0, vsf01, vsf02 above)
bestPhisXSF = [0 bestPhisXSF];
bestModesXSF = [{'xsfxb30'} bestModesXSF];

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

plot(bestAccs,bestPeaks,'DisplayName','XSF*','LineWidth',2)

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
simphis = 0:5:100;
NN = 1e4;
subplot(2,1,2); hold on
allModes = fields(modes);
plotModes = {'s1','s3','sf'};
for mode=plotModes
    modeName = mode{:};
    fprintf(['Mode: ' modeName '\n'])
    thisM = modes.(modeName);
    psvs = [];
    accs = [];
    for p=simphis
        fprintf(' Phi=%d:',p)
        if isfield(thisM(c(p)),'pot') && ~isempty(thisM(c(p)).pot)
            [num, psv] = simEffTrap(thisM(c(p)).pot,'num',NN);
            psvs = [psvs psv];
            accs = [accs thisM(c(p)).acc];
        elseif p>90
            if strcmp(modeName,'s3')
                a = 4*((p-90)/5) + 86;
            else
                a = 8*((p-90)/5) + 259;
            end
            a = a * 1e3; % convert to m/s/s
            pot = efftrap3Dgen(thisM(1).decels,thisM(1).phifunc(90),[1 1],a);
            [num, psv] = simEffTrap(pot,'num',NN);
            psvs = [psvs psv];
            accs = [accs a];
        end
        fprintf(' %d, %1.1f\n',num,psv)
    end
    accs = accs / 1e3; % back to km/s/s
    plot(accs,psvs,'DisplayName',upper(modeName),'LineWidth',2)
end

xlabel('Deceleration (km/s/s)','FontSize',13)
ylabel('Phase Space Volume (m^6/s^3)','FontSize',13)
title('Monte-Carlo in Effective Trap','FontSize',14)
xlim([0 400])
set(gca,'FontSize',13)
hl = legend('show');
hl.FontSize = 13;


%% Next add in the VSF
%figure;
%hold on
accs = bestPhisVSF;
psvs = bestPhisVSF;
for i=1:length(bestPhisVSF)
    fprintf(['Mode: ' bestModesVSF{i} ' Phi: ' num2str(bestPhisVSF(i))])
    thisM = modes.(bestModesVSF{i});
    thisM = thisM(c(bestPhisVSF(i)));
    accs(i) = thisM.acc;
    [num, psv] = simEffTrap(thisM.pot,'num',NN);
    fprintf(', N: %d, PSV: %f\n',num,psv)
    psvs(i) = psv;
end

% Add in some higher phase angles
for a=[300e3 310e3]
    mode = modes.(bestModesVSF{end});
    mode = mode(90);
    pot = efftrap3Dgen(mode.decels,mode.phifunc(90),[1 1],a);
    [num, psv] = simEffTrap(pot,'num',NN);
    psvs = [psvs psv];
    accs = [accs a];
end
            
accs = accs * 1e-3;
plot(accs,psvs,'DisplayName','VSF*','LineWidth',2)
%% And the XSF
accs = bestPhisXSF;
psvs = bestPhisXSF;
for i=1:length(bestPhisXSF)
    fprintf(['Mode: ' bestModesXSF{i} ' Phi: ' num2str(bestPhisXSF(i))])
    thisM = modes.(bestModesXSF{i});
    thisM = thisM(c(bestPhisXSF(i)));
    accs(i) = thisM.acc;
    [num, psv] = simEffTrap(thisM.pot,'num',NN,'maxVel',17);
    fprintf(', N: %d, PSV: %f\n',num,psv)
    psvs(i) = psv;
end    
    
% Add in some higher phase angles
for a=[370e3]
    mode = modes.(bestModesXSF{end});
    mode = mode(90);
    pot = efftrap3Dgen(mode.decels,mode.phifunc(90),[1 1],a);
    [num, psv] = simEffTrap(pot,'num',NN,'maxVel',17);
    psvs = [psvs psv];
    accs = [accs a];
end

accs = accs*1e-3;
plot(accs,psvs,'DisplayName','XSF*','LineWidth',2)

    
%% Now figure out how to get Traveling Wave in here.    
raccs
tpots
psvs = raccs;
i = 1;
for tp=tpots
    fprintf(' Acc=%d: ',raccs(i))
    [num, psv] = simEffTrapRing(tp{:},'num',NN,'maxVel',17);
    psvs(i) = psv; i = i+1;
    fprintf('N=%d, PSV=%f\n',num,psv)
end

plot(raccs,psvs/pi,'DisplayName','TW','LineWidth',2)





