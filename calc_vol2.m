% with more serious integral calculation
function vol=cal_vol2(outmod)

all=outmod;
allz=squeeze(all(40,40,:));
allx=squeeze(all(40,40:end,121));

%figure;plot(allx)
x= (0:.025:.975)*1e-3;
y=x;
z=(-3:.025:3)*1e-3;
[xx,yy,zz]=meshgrid(x,y,z);
potential=@(x,y,z)interp3(xx,yy,zz,all(40:end,40:end,:),x,y,z);

% find x direction and z direction trap depth

[pksz,locsz] = findpeaks(allz);
[pksx,locsx] = findpeaks(allx);

depthz=min(pksz);
if isempty(pksx)
        depthx=max(allx);
else
    a
    depthx=min(pksx);
end

if isempty(pksz)
        depthz=max(allz);
       locsz=241;

end

depth=min([depthx,depthz]);

% p=isosurface(xx,yy,zz,all,depth);
% if max(size(locsz))>1
%     
%     z_bound=(p.vertices(:,3)<z(locsz(end)))&(p.vertices(:,3)>z(locsz(1)));
% else
%     z_bound=p.vertices(:,3)<z(locsz);
% end
%figure;plot3(p.vertices(z_bound,1),p.vertices(z_bound,2),p.vertices(z_bound,3),'.')
if depthx>=depthz
    x_bound=interp1(allx,x,depth);
    
    z_bound_l=interp1(allz(1:121),z(1:121),depthz);
    z_bound_h=z(locsz(end));
else
    x_bound=max(x);
    z_bound_l=interp1(allz(1:121),z(1:121),depth);
    z_bound_h=interp1(allz(121:locsz(end)),z(121:locsz(end)),depth);
end

x_new=0:x_bound/200:x_bound;
y_new=x_new;
z_new=z_bound_l:(z_bound_h-z_bound_l)/400:z_bound_h;

%calculation of trap volume
% for xxx=1:51
%     for yyy=1:51
%         for zzz=1:101
%             
%             all_new(xxx,yyy,zzz)=potential(x_new(xxx),y_new(yyy),z_new(zzz));
%             
%         end             
%     end       
% end

[xxx,yyy,zzz]=meshgrid(x_new,y_new,z_new);
all_new=interp3(xx,yy,zz,all(40:end,40:end,:),xxx,yyy,zzz);


volume=(all_new<=depth);
volume_all=(depth-all_new).^1.5.*volume;
vol=sum(sum(sum(volume_all)));