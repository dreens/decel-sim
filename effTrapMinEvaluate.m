%%
% Make use of the Min Depth Finder to study some plots!
%
phis = [1, 5:5:85, 89];

%% All of the needed potentials and min depths.
if ~exist('outs1')
    outs1 = {};
    ous1a = [];
    for phi=phis
        [outs1{phi}, outs1a(phi)] = efftrap3Dgen({'longdecel','longdecel'},[phi-180,0,phi],[1 1]);
    end
end
if ~exist('outs3')
    outs3 = {};
    outs3a = [];
    for phi=phis
        [outs3{phi},outs3a(phi)] = efftrap3Dgen({'longdecel','longdecel'},[phi-540,0,phi],[1 1]);
    end
end
if ~exist('outpggg')
    outpggg = {};
    outpggga = [];
    for phi=phis
        [outpggg{phi},outpggga(phi)] = efftrap3Dgen({'singlerod','longdecel'},[phi-180,-phi,phi],[1 1]);
    end
end
if ~exist('outppgg')
    outppgg = {};
    outppgga = [];
    for phi=phis
        [outppgg{phi}, outppgga(phi)] = efftrap3Dgen({'ppgg','longdecel'},[phi-180,-phi,phi],[1 1]);
    end
end
if ~exist('outppmm')
    outppmm = {};
    outppmma = [];
    for phi=phis
        [outppmm{phi}, outppmma(phi)] = efftrap3Dgen({'ppmm_2mm','longdecel'},[phi-180,-phi,phi],[1 1]);
    end
end
if ~exist('outppggs')
    outppggs = {};
    outppggsa = [];
    for phi=phis
        [outppggs{phi}, outppggsa(phi)] = efftrap3Dgen({'ppgg','longdecel'},[phi-180,phi-60,phi],[1 1]);
    end
end
if ~exist('outppmms')
    outppmms = {};
    outppmmsa = [];
    for phi=phis
        [outppmms{phi}, outppmmsa(phi)] = efftrap3Dgen({'ppmm_2mm','longdecel'},[phi-180,phi-60,phi],[1 1]);
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
    
    
    
    