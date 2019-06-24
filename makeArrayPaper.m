phi = 45;
out = efftrap3D('longdecel','longdecel',phi-180,-phi,phi);
outs3 = efftrap3D('longdecel','longdecel',phi-540,-phi,phi);
outpggg = efftrap3D('longdecel','singlerod',phi-180,-phi,phi);
outppgg = efftrap3D('longdecel','ppgg',phi-180,-phi,phi);
outppmm = efftrap3D('longdecel','ppmm_2mm',phi-180,-phi,phi);
out2p2p = efftrap3D('singlerod','ppmm_2mm',phi-180,-phi,phi);
%%
plotTrapSlice({outppmm,outppgg,outpggg,outs3,out},0,'contourArray')
