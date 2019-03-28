clear all_out
j=1;
for phi = 0.1:5:79
phi2=phi;
%[out,accel] = efftrap3D('longdecel','longdecel',phi-180,-phi,phi);
%[outpggg,accel] = efftrap3D('longdecel','singlerod',phi-180,-phi,phi);
%[outppgg,accel] = efftrap3D('longdecel','ppgg',phi-180,-phi,phi);
%[outmod2,accel] = efftrap3Dgen({'longdecel','ppgg','longdecel','longdecel','ppgg','longdecel'},[-180,-162+phi,-18,0,18+phi,162,180],[0 0 0 1 0 1 ]);
 %outmod = efftrap3Dgen({'longdecel','ppgg','longdecel','longdecel','ppgg','longdecel'},[-180,-180+phi,-18,0,phi,162,180],[0 0 1 1 1 0]);
[outs3,accel] = efftrap3D('longdecel','longdecel',phi-540,-phi,phi);
%[outvsfb10,accel] = efftrap3D('None','ppgg',-155,-25,25);
%[outmod4,accel] = efftrap3Dgen({'longdecel','ppgg','longdecel','longdecel','ppgg','longdecel'},[-180,-180+phi,phi2-90,0,phi,90+phi2,180],[0 0 1 1 1 0]);
all_out(j)=cal_vol2(outs3);
all_acel(j)=accel;
j=j+1;
end

figure(99);hold on;subplot(2,1,1);plot( 0.1:5:79,all_out,'DisplayName','putpggg')
subplot(2,1,2);hold on;plot( 0.1:5:79,all_acel,'DisplayName','outpggg')
%%
phi = 80;
phi2=80;
%out = efftrap3D('longdecel','longdecel',phi-180,-phi,phi);
outpggg = efftrap3D('longdecel','singlerod',phi-180,-phi,phi);
outppgg = efftrap3D('longdecel','ppgg',phi-180,-phi,phi);

 
outmod = efftrap3Dgen({'longdecel','ppgg','longdecel','longdecel','ppgg','longdecel'},[-180,-180+phi,phi2-90,0,phi,90+phi2,180],[0 0 1 1 1 0]);

plotTrapSlice({outmod ,outppgg,outpggg,out},0,'contourArray')