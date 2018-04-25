phi = 45;
out = efftrap3D('longdecel','longdecel',phi-180,-phi,phi);
outpggg = efftrap3D('longdecel','singlerod',phi-180,-phi,phi);
outppgg = efftrap3D('longdecel','ppgg',phi-180,-phi,phi);
outpmpm = efftrap3D('pmpm_2mm_no-sym','ppgg',-110,-20,70);
%plotTrapSlice({outppgg,outpggg,out},0,'contourArray')

%%

outs3 = efftrap3D('longdecel','longdecel',phi-540,-phi,phi);
%plotTrapSlice(outs3,39,'yz')
%%
plotTrapSlice({outppgg,outpggg,out,out},39,'contourArray');

%% try generalized efftrap
outmod = efftrap3Dgen({'ppgg','longdecel','ppgg','pmpm_2mm_no-sym'},[-124 -46 56 160 236],[1 1 0 0]);
plotTrapSlice({outmod, outpmpm, outppgg,outpggg,outs3,out},39,'contourArray');
