%%
% Make use of the Min Depth Finder to study some plots!
%
phis = [3:3:87];
%outs3 = {};
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
    midphi = phi-180;
    if phi>length(outs3) || isempty(outs3{phi})
        out{phi} = efftrap3Dgen({'longdecel','longdecel'},[phi-180,midphi,phi],[1 1]);
        %outs3{phi} = efftrap3Dgen({'longdecel','longdecel'},[phi-540,midphi,phi],[1 1]);
        outpggg{phi} = efftrap3Dgen({'singlerod','longdecel'},[phi-180,midphi,phi],[1 1]);
        outppgg{phi} = efftrap3Dgen({'ppgg','longdecel'},[phi-180,midphi,phi],[1 1]);
        outppmm{phi} = efftrap3Dgen({'ppmm_2mm','longdecel'},[phi-180,midphi,phi],[1 1]);
    end
    s1(phi)  = effTrapMinDepth(out{phi});
    %s3(phi)  = effTrapMinDepth(outs3{phi});
    sf(phi)  = effTrapMinDepth(outpggg{phi});
    vsf(phi) = effTrapMinDepth(outppgg{phi});
    xsf(phi) = effTrapMinDepth(outppmm{phi});
end

figure;
phis = 1:29;
plot(phis,s3(phis*3),'DisplayName','S=3')
hold on
%plot(phis,sf(phis),'DisplayName','SF')
%plot(phis,vsf(phis),'DisplayName','VSF')
%plot(phis,xsf(phis),'DisplayName','XSF')
legend('show')
