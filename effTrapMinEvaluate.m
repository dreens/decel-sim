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
    out = efftrap3Dgen({'longdecel','longdecel'},[phi-180,-phi,phi],[1 1]);
    s1 = [s1 effTrapMinDepth(out)];
    outpggg = efftrap3Dgen({'singlerod','longdecel'},[phi-180,-phi,phi],[1 1]);
    sf = [sf effTrapMinDepth(outpggg)];
    outppgg = efftrap3Dgen({'ppgg','longdecel'},[phi-180,-phi,phi],[1 1]);
    vsf = [vsf effTrapMinDepth(outppgg)];
    outppmm = efftrap3Dgen({'ppmm_2mm','longdecel'},[phi-180,-phi,phi],[1 1]);
    xsf = [xsf effTrapMinDepth(outppmm)];
end

figure;
plot(phis,s1,'DisplayName','S=1')
hold on
plot(phis,sf,'DisplayName','SF')
plot(phis,vsf,'DisplayName','VSF')
plot(phis,xsf,'DisplayName','XSF')
legend('show')
    