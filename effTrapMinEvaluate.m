%%
% Make use of the Min Depth Finder to study some plots!
%
phis = [3:3:87];
outs3 = {};
if false%~exist('out')
    out = {};
    outpggg = {};
    outppgg = {};
    outppmm = {};
    s1 = [];
    sf = [];
    vsf = [];
    xsf = [];
end
for phi=phis
    phi
    midphi = -phi;
    if phi>length(out) || isempty(out{phi})
        out{phi} = efftrap3Dgen({'longdecel','longdecel'},[phi-180,midphi,phi],[1 1]);
        outpggg{phi} = efftrap3Dgen({'singlerod','longdecel'},[phi-180,midphi,phi],[1 1]);
        outppgg{phi} = efftrap3Dgen({'ppgg','longdecel'},[phi-180,midphi,phi],[1 1]);
        outppmm{phi} = efftrap3Dgen({'ppmm_2mm','longdecel'},[phi-180,midphi,phi],[1 1]);
    end
    s1(phi)  = effTrapMinDepth(out{phi});
    sf(phi)  = effTrapMinDepth(outpggg{phi});
    vsf(phi) = effTrapMinDepth(outppgg{phi});
    xsf(phi) = effTrapMinDepth(outppmm{phi});
end

figure;
phis = 1:89;
plot(phis,s1(phis),'DisplayName','S=1')
hold on
plot(phis,sf(phis),'DisplayName','SF')
plot(phis,vsf(phis),'DisplayName','VSF')
plot(phis,xsf(phis),'DisplayName','XSF')
legend('show')
    