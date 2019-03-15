%%
% Make use of the Min Depth Finder to study some plots!
%
phis = 5:5:80;
s1 = [];
sf = [];
vsf = [];
xsf = [];
for phi=phis
    phi
    if ~exist('out')
        out{phi} = efftrap3Dgen({'longdecel','longdecel'},[phi-180,-phi,phi],[1 1]);
        outpggg{phi} = efftrap3Dgen({'singlerod','longdecel'},[phi-180,-phi,phi],[1 1]);
        outppgg{phi} = efftrap3Dgen({'ppgg','longdecel'},[phi-180,-phi,phi],[1 1]);
        outppmm{phi} = efftrap3Dgen({'ppmm_2mm','longdecel'},[phi-180,-phi,phi],[1 1]);
    end
    s1 = [s1 effTrapMinDepth(out{phi})];
    sf = [sf effTrapMinDepth(outpggg{phi})];
    vsf = [vsf effTrapMinDepth(outppgg{phi})];
    xsf = [xsf effTrapMinDepth(outppmm{phi})];
end

figure;
plot(phis,s1,'DisplayName','S=1')
hold on
plot(phis,sf,'DisplayName','SF')
plot(phis,vsf,'DisplayName','VSF')
plot(phis,xsf,'DisplayName','XSF')
legend('show')
    