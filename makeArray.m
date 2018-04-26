phi = 45;
out = efftrap3D('longdecel','longdecel',phi-180,-phi,phi);
outpggg = efftrap3D('longdecel','singlerod',phi-180,-phi,phi);
outppgg = efftrap3D('longdecel','ppgg',phi-180,-phi,phi);
outppmm = efftrap3D('longdecel','ppmm_2mm',phi-180,-phi,phi);
outpmpm = efftrap3D('pmpm_2mm_no-sym','ppgg',-109,-19,71);
%plotTrapSlice({outppgg,outpggg,out},0,'contourArray')

%%

outs3 = efftrap3D('longdecel','longdecel',phi-540,-phi,phi);
%plotTrapSlice(outs3,39,'yz')
%%
plotTrapSlice({outppmm,outppgg,outpggg,out},39,'contourArray');

%% try generalized efftrap
outmod = efftrap3Dgen({'ppgg','longdecel','ppgg','pmpm_2mm_no-sym'},[-132 -50 68 160 228],[1 1 0 0]);
plotTrapSlice({outmod, outpmpm, outppgg,outpggg},39,'contourArray');
