function vol=cal_vol(outmod)

all=outmod;
allz=squeeze(all(40,40,:));
allx=squeeze(all(40,:,121));
dEdz=gradient(allz);
dEdx=gradient(allx);
index_z=sign(gradient(allz));
index_x=sign(gradient(allx));

for i=1:max(size(index_z))-1
    total_z(i)=index_z(i)+index_z(i+1);
    
end
for i=1:max(size(index_x))-1
    total_x(i)=index_x(i)+index_x(i+1);
    
end

   allindex_z=find(total_z~=2&total_z~=-2);
   m=1;
   min_index_z=[];
for j=allindex_z
   
    if    sign(dEdz(j))>0&sign(dEdz(j+1))<0
        
        min_index_z(m)=j;m=m+1;
    end
    
end

if isempty(min_index_z)
     minz_po=min([allz(1),allz(end)]);
     min_index_z=241;
else
  minz_po=min(allz(min_index_z));
 
end

% x/y trap depth
  allindex_x=find(total_x~=2&total_x~=-2);
   m=1;
   min_index_x=[];
for j=allindex_x
   
    if    sign(dEdx(j))>0&sign(dEdx(j+1))<0
        
        min_index_x(m)=j;m=m+1;
    end
    
end
if isempty(min_index_x)
     minx_po=min([allx(1),allx(end)]);
else
  minx_po=min(allx(min_index_x));
 
end
min_po=min([minx_po,minz_po]);
%calculation of trap volume
for xx=1:79
    for yy=1:79
        for zz=1:241
            
            comstrain(xx,yy,zz)=zz;
            
        end             
    end       
end
if max(size(min_index_z))>1
    num_con=comstrain<min_index_z(end)  ;
    num_con=num_con&(comstrain>min_index_z(1));
else
    
    num_con=(comstrain<min_index_z);    
    
end

volume=(all<=min_po);
volume_all=(min_po-all).^1.5.*volume.*num_con;
vol=sum(sum(sum(volume_all)));