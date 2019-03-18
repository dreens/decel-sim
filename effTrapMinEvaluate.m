%%
% Make use of the Min Depth Finder to study some plots!
%

phis = [1, 5:5:85, 89];

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
    modes.sf = struct();
    modes.sf.decels = {'ppgg','longdecel'};
    modes.sf.phifunc = @(p) [p-180, -p, p];
end
if ~isfield(modes,'xsf')
    modes.sf = struct();
    modes.sf.decels = {'ppmm_2mm','longdecel'};
    modes.sf.phifunc = @(p) [p-180, -p, p];
end
if ~isfield(modes,'vsfe')
    modes.sf = struct();
    modes.sf.decels = {'ppgg','longdecel'};
    modes.sf.phifunc = @(p) [p-180, p-60, p];
end
if ~isfield(modes,'xsfe')
    modes.sf = struct();
    modes.sf.decels = {'ppmm_2mm','longdecel'};
    modes.sf.phifunc = @(p) [p-180, p-60, p];
end

%% All of the needed potentials and min depths.
allmodes = fields(modes);
for i=1:length(allmodes)
    modes.(allmodes(i))
    for phi=phis
        [outs1.pot{phi}, outs1.acc(phi)] = efftrap3Dgen({'longdecel','longdecel'},[phi-180,0,phi],[1 1]);
        [outs1.dep outs1.typ outs1.vol outs1.psv] = effTrapMinDepth(outs1.pot);
    end
end
if ~exist('outs3')
    outs3 = struct();
    for phi=phis
        [outs3.pot,outs3.acc] = efftrap3Dgen({'longdecel','longdecel'},[phi-540,0,phi],[1 1]);
        [outs3.dep outs3.typ outs3.vol outs3.psv] = effTrapMinDepth(outs1.pot);
    end
end
if ~exist('outpggg')
    outpggg = struct();
    for phi=phis
        [outpggg.pot,outpggg.acc] = efftrap3Dgen({'singlerod','longdecel'},[phi-180,-phi,phi],[1 1]);
    end
end
if ~exist('outppgg')
    outppgg = struct();
    for phi=phis
        [outppgg.pot, outppgg.acc] = efftrap3Dgen({'ppgg','longdecel'},[phi-180,-phi,phi],[1 1]);
    end
end
if ~exist('outppmm')
    outppmm = struct();
    for phi=phis
        [outppmm.pot, outppmm.acc] = efftrap3Dgen({'ppmm_2mm','longdecel'},[phi-180,-phi,phi],[1 1]);
    end
end
if ~exist('outppggs')
    outppggs = struct();
    for phi=phis
        [outppggs.pot, outppggs.acc] = efftrap3Dgen({'ppgg','longdecel'},[phi-180,phi-60,phi],[1 1]);
    end
end
if ~exist('outppmms')
    outppmms = struct();
    for phi=phis
        [outppmms.pot, outppmms.acc] = efftrap3Dgen({'ppmm_2mm','longdecel'},[phi-180,phi-60,phi],[1 1]);
    end
end



for phi=phis
    phi
    midphi = phi-120;
    if phi>length(outacc) || isempty(outacc{phi})
        [~, outacc{phi}] = efftrap3Dgen({'longdecel','longdecel'},[phi-180,midphi,phi],[1 1]);
        %outs3{phi} = efftrap3Dgen({'longdecel','longdecel'},[phi-540,midphi,phi],[1 1]);
        %outpggg{phi} = efftrap3Dgen({'singlerod','longdecel'},[phi-180,midphi,phi],[1 1]);
        %[outppggs{phi}, outppggsacc{phi}] = efftrap3Dgen({'ppgg','longdecel'},[phi-180,midphi,phi],[1 1]);
        %[outppmms{phi}, outppmmsacc{phi}] = efftrap3Dgen({'ppmm_2mm','longdecel'},[phi-180,midphi,phi],[1 1]);
    end
    %s1(phi)  = effTrapMinDepth(out{phi});
    %s3(phi)  = effTrapMinDepth(outs3{phi});
    %sf(phi)  = effTrapMinDepth(outpggg{phi});
    %vsf(phi) = effTrapMinDepth(outppgg{phi});
    %xsf(phi) = effTrapMinDepth(outppmm{phi});
    %vsfp(phi) = effTrapMinDepth(outppggs{phi});
    %xsfp(phi) = effTrapMinDepth(outppmms{phi});
end
%%
figure;
phis = [1 5:5:20 70:5:85 89];
plot(phis,vsfp(phis),'DisplayName','VSF''')
hold on
plot(phis,xsfp(phis),'DisplayName','XSF''')
%plot(phis,sf(phis),'DisplayName','SF')
%plot(phis,vsf(phis),'DisplayName','VSF')
%plot(phis,xsf(phis),'DisplayName','XSF')
legend('show')

%%
figure
plot(3:3:87,[outacc{3:3:87}])
    
    
    
    