phi = 55;
out = efftrap3D('longdecel','longdecel',phi-180,-phi,phi);
outpggg = efftrap3D('longdecel','singlerod',phi-180,-phi,phi);
outppgg = efftrap3D('longdecel','ppgg',phi-180,-phi,phi);
outppmm = efftrap3D('longdecel','ppmm_2mm',phi-180,-phi,phi);
%outpmpm = efftrap3D('pmpm_2mm_no-sym','ppgg',-109,-19,71);
outppmm2 = efftrap3D('pmpm_2mm_no-sym','ppmm_2mm',-112.5,-20,67.5);
%plotTrapSlice({outppgg,outpggg,out},0,'contourArray')

%%

outs3 = efftrap3D('longdecel','longdecel',phi-540,-phi,phi);
%see if delay switching im S=3 matters. It doesnt. 10% effect maybe. 
outs3p = efftrap3D('longdecel','singlerod',phi-540,-360-phi,phi);

%plotTrapSlice(outs3,39,'yz')
%%
plotTrapSlice({outppmm,outppgg,outpggg,out},39,'contourArray');
plotTrapSlice({outs3p,outs3,outpggg,out},39,'contourArray');


%% try generalized efftrap
outmod = efftrap3Dgen({'ppgg','longdecel','ppgg','pmpm_2mm_no-sym'},[-132 -50 68 160 228],[1 1 0 0]);
plotTrapSlice({outppmm2, outppmm, outppgg,outpggg},39,'contourArray');

%% DC guiding.
guide = efftrap3D('ppgg','ppgg',-180,0,180);
guideppmm = efftrap3D('ppmm_2mm','ppmm_2mm',-180,0,180);
plotTrapSlice(guideppmm,59,'xy')