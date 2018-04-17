phi = 45;
out = efftrap3D('longdecel','longdecel',phi-180,-phi,phi);
outpggg = efftrap3D('longdecel','singlerod',phi-180,-phi,phi);
outppgg = efftrap3D('longdecel','ppgg',phi-180,-phi,phi);
%plotTrapSlice({outppgg,outpggg,out},0,'contourArray')

%%

outs3 = efftrap3D('longdecel','longdecel',phi-540,-phi,phi);
%plotTrapSlice(outs3,39,'yz')
%%
plotTrapSlice({outppgg,outpggg,outs3,out},39,'contourArray')