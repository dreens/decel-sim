phi = 45;
out = efftrap3D('longdecel','longdecel',phi-180,-phi,phi);
outs3 = efftrap3D('longdecel','longdecel',phi-540,-phi,phi);
outpggg = efftrap3D('longdecel','singlerod',phi-180,-phi,phi);
outppgg = efftrap3D('longdecel','ppgg',phi-180,-phi,phi);
outppmm = efftrap3D('longdecel','ppmm_2mm',phi-180,-phi,phi);
%outpmpm = efftrap3D('pmpm_2mm_no-sym','ppgg',-109,-19,71);
%outppmm2 = efftrap3D('pmpm_2mm_no-sym','ppmm_2mm',-112.5,-20,67.5);
plotTrapSlice({outppmm,outppgg,outpggg,outs3,out},0,'contourArray')
%plotTrapSlice({guide},0,'contourArray')
%%

outs3 = efftrap3D('longdecel','longdecel',phi-540,-phi,phi);
%see if delay switching im S=3 matters. It doesnt. 10% effect maybe. 
outs3p = efftrap3D('longdecel','singlerod',phi-540,-360-phi,phi);

plotTrapSlice(outs3,39,'yz')
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

%% Trap Strength
phi = 55;
outppgg1 = efftrap3D('longdecel','singlerod',phi-180,-phi,phi);
out = efftrap3D('longdecel','longdecel',phi-180,-phi,phi);


%% VSF Bunching
outvsfb10 = efftrap3D('None','ppgg',-125,-35,55);

%%
plotTrapSlice(outppgg,0,'1D')
hold on
x = (-2:.01:2) *1e-3;
v = 200 * abs(x);
v = 1e3 * v * 1.4 * 9e-24 / 1.38e-23;
plot(x,v/2,'DisplayName','1T/cm')
plot(x,v,'DisplayName','2T/cm')
legend('show')
xlim([-2e-3 2e-3])
%%
plotTrapSlice({outvsfb10,outvsfb,outppgg1,out},0,'contourArray')
%%
phi=0;
outmod = efftrap3Dgen({'longdecel','ppgg','longdecel','longdecel','ppgg','longdecel'},[-180,-162,-18,0,18,162,180],[0 0 0 1 1 1 ]);

plotTrapSlice({outmod},0,'contourArray')

