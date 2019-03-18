%%
% Make use of the Min Depth Finder to study some plots!
%

phis = [-15:5:85, 89];
check = @(p) p + (p<=0)*180;

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
if ~isfield(modes,'vsfe20')
    modes.vsfe20 = struct();
    modes.vsfe20.decels = {'ppgg','longdecel'};
    modes.vsfe20.phifunc = @(p) [p-180, (p<45)*(p-20)+(p>45)*(p-160), p];
end
if ~isfield(modes,'xsfe20')
    modes.xsfe20 = struct();
    modes.xsfe20.decels = {'ppmm_2mm','longdecel'};
    modes.xsfe20.phifunc = @(p) [p-180, (p<45)*(p-20)+(p>45)*(p-160), p];
end
if ~isfield(modes,'vsfe40')
    modes.vsfe40 = struct();
    modes.vsfe40.decels = {'ppgg','longdecel'};
    modes.vsfe40.phifunc = @(p) [p-180, (p<45)*(p-40)+(p>45)*(p-140), p];
end
if ~isfield(modes,'xsfe40')
    modes.xsfe40 = struct();
    modes.xsfe40.decels = {'ppmm_2mm','longdecel'};
    modes.xsfe40.phifunc = @(p) [p-180, (p<45)*(p-40)+(p>45)*(p-140), p];
end
if ~isfield(modes,'vsfe60')
    modes.vsfe60 = struct();
    modes.vsfe60.decels = {'ppgg','longdecel'};
    modes.vsfe60.phifunc = @(p) [p-180, (p<45)*(p-60)+(p>45)*(p-120), p];
end
if ~isfield(modes,'xsfe60')
    modes.xsfe60 = struct();
    modes.xsfe60.decels = {'ppmm_2mm','longdecel'};
    modes.xsfe60.phifunc = @(p) [p-180, (p<45)*(p-60)+(p>45)*(p-120), p];
end
if ~isfield(modes,'vsfe80')
    modes.vsfe80 = struct();
    modes.vsfe80.decels = {'ppgg','longdecel'};
    modes.vsfe80.phifunc = @(p) [p-180, (p<45)*(p-80)+(p>45)*(p-100), p];
end
if ~isfield(modes,'xsfe80')
    modes.xsfe80 = struct();
    modes.xsfe80.decels = {'ppmm_2mm','longdecel'};
    modes.xsfe80.phifunc = @(p) [p-180, (p<45)*(p-80)+(p>45)*(p-100), p];
end
if ~isfield(modes,'vsfe10')
    modes.vsfe10 = struct();
    modes.vsfe10.decels = {'ppgg','longdecel'};
    modes.vsfe10.phifunc = @(p) [p-180, (p<45)*(p-10)+(p>45)*(p-170), p];
end
if ~isfield(modes,'xsfe10')
    modes.xsfe10 = struct();
    modes.xsfe10.decels = {'ppmm_2mm','longdecel'};
    modes.xsfe10.phifunc = @(p) [p-180, (p<45)*(p-10)+(p>45)*(p-170), p];
end
if ~isfield(modes,'vsfe30')
    modes.vsfe30 = struct();
    modes.vsfe30.decels = {'ppgg','longdecel'};
    modes.vsfe30.phifunc = @(p) [p-180, (p<45)*(p-30)+(p>45)*(p-150), p];
end
if ~isfield(modes,'xsfe30')
    modes.xsfe30 = struct();
    modes.xsfe30.decels = {'ppmm_2mm','longdecel'};
    modes.xsfe30.phifunc = @(p) [p-180, (p<45)*(p-30)+(p>45)*(p-150), p];
end
if ~isfield(modes,'vsfe50')
    modes.vsfe50 = struct();
    modes.vsfe50.decels = {'ppgg','longdecel'};
    modes.vsfe50.phifunc = @(p) [p-180, (p<45)*(p-50)+(p>45)*(p-130), p];
end
if ~isfield(modes,'xsfe50')
    modes.xsfe50 = struct();
    modes.xsfe50.decels = {'ppmm_2mm','longdecel'};
    modes.xsfe50.phifunc = @(p) [p-180, (p<45)*(p-50)+(p>45)*(p-130), p];
end
if ~isfield(modes,'vsfe70')
    modes.vsfe70 = struct();
    modes.vsfe70.decels = {'ppgg','longdecel'};
    modes.vsfe70.phifunc = @(p) [p-180, (p<45)*(p-70)+(p>45)*(p-110), p];
end
if ~isfield(modes,'xsfe70')
    modes.xsfe70 = struct();
    modes.xsfe70.decels = {'ppmm_2mm','longdecel'};
    modes.xsfe70.phifunc = @(p) [p-180, (p<45)*(p-70)+(p>45)*(p-110), p];
end
if ~isfield(modes,'vsfe90')
    modes.vsfe90 = struct();
    modes.vsfe90.decels = {'ppgg','longdecel'};
    modes.vsfe90.phifunc = @(p) [p-180, (p<45)*(p-90)+(p>45)*(p-90), p];
end
if ~isfield(modes,'xsfe90')
    modes.xsfe90 = struct();
    modes.xsfe90.decels = {'ppmm_2mm','longdecel'};
    modes.xsfe90.phifunc = @(p) [p-180, (p<45)*(p-90)+(p>45)*(p-90), p];
end

%% All of the needed potentials and min depths.
allmodes = fields(modes);
for i=1:length(allmodes)
    fprintf('Mode %s:\n',allmodes{i})
    thisM = modes.(allmodes{i});
    if length(thisM) < 180
        thisM(2:180) = thisM(1);
    end
    for p=phis
        fprintf(' Phi=%d\n',p)
        if ~isfield(thisM(p),'pot') || isempty(thisM(p).pot)
            [thisM(p).pot, thisM(p).acc] = ...
                efftrap3Dgen(thisM(p).decels,thisM(p).phifunc(p),[1 1]);
        end
        if ~isfield(thisM(p),'dep') || isempty(thisM(p).dep)
            [thisM(p).dep, thisM(p).typ, thisM(p).vol, thisM(p).psv] = ...
                effTrapMinDepth(thisM(p).pot);
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

%%
figure; hold on
allModes = fields(modes);
for i=1:length(allModes)
    modeName = allModes{i};
    if length(modeName)~=6 || ~strcmp(modeName(end-2),'e')
        thisM = modes.(modeName);
        plot([thisM(:).acc]/1e3,[thisM(:).dep]/1.38e-23/1e-3,'DisplayName',allModes{i})
    end
end

xlabel('Deceleration (km/s/s)')
ylabel('Trap Depth (mK)')
title('Effective Trap Depths')

%% Optimized envelope for XSF, VSF mods
% Let's just grab the peaks.
figure;
allModes = fields(modes);
locs = [];
peaks = [];
for i=1:length(allModes)
    modeName = allModes{i};
    if length(modeName)==6 && strcmp(modeName(1:2),'xs')
        thisM = modes.(modeName);
        depths = [thisM.dep];
        accs = [thisM.acc];
        deptweak = (1:length(depths))*1.38e-23*1e-3;
        [~, loc] = max(depths+deptweak);
        peaks(end+1) = depths(loc);
        locs(end+1) = accs(loc);
    end
end

figure;
plot(locs,peaks,'.')
    