
function [nn, ind_max]=findNN(ga,conn,threshold, mask)

if conn==18
    th=1.6;
elseif conn==26
    th=3;
elseif conn==124
    th=5;
elseif conn==6
    th=1.2;
end

[nx ny nz]= size(ga);
ind_max=find(ga>threshold);

        
[p1,p2,p3]=ind2sub([nx ny nz],ind_max);

if conn==124
[x,y,z]=meshgrid(-2:1:2);
else
[x,y,z]=meshgrid(-1:1:1);
end
x=x(:);
y=y(:);
z=z(:);

check_xyz=sqrt(x(:).^2+y(:).^2+z(:).^2);
ind_neigh_to_exclude=find(check_xyz>th);


x(ind_neigh_to_exclude)=[];
y(ind_neigh_to_exclude)=[];
z(ind_neigh_to_exclude)=[];

x=repmat(x(:)',[ size(p1,1)  1]);
y=repmat(y(:)',[ size(p2,1)  1]);
z=repmat(z(:)',[ size(p3,1)  1]);

p1_repmat27=repmat(p1,[1 conn+1])-x;
p2_repmat27=repmat(p2,[1 conn+1])-y;
p3_repmat27=repmat(p3,[1 conn+1])-z;


p1_repmat27(p1_repmat27<1)=1;
p1_repmat27(p1_repmat27>nx)=nx;

p2_repmat27(p2_repmat27<1)=1;
p2_repmat27(p2_repmat27>ny)=ny;

p3_repmat27(p3_repmat27<1)=1;
p3_repmat27(p3_repmat27>nz)=nz;
% 
% [ind_max3(:,1), ind_max3(:,2), ind_max3(:,3)] = ind2sub([nx ny nz],ind_max); 
% 
% dist1=repmat(ind_max3(:,1),1,size(p1_repmat27,2))-p1_repmat27;
% dist2=repmat(ind_max3(:,2),1,size(p2_repmat27,2))-p2_repmat27;
% dist3=repmat(ind_max3(:,3),1,size(p3_repmat27,2))-p3_repmat27;
% 
% %intensityDiff = 
% dist=sqrt(dist1.*dist1+dist2.*dist2+dist3.*dist3);
% [sorted,ind_sorted]=sort(dist,2);

% p1_repmat27 = p1_repmat27(ind_sorted);
% p2_repmat27 = p2_repmat27(ind_sorted);
% p3_repmat27 = p3_repmat27(ind_sorted);



 for i=1:size(p1_repmat27,1)
    flag_this_voxel_is_corected=0;
    j=1;
    while j<=size(p1_repmat27,2) && ~flag_this_voxel_is_corected
        current_ind_nearby=sub2ind([nx ny nz], p1_repmat27(i,j) ,p2_repmat27(i,j) ,p3_repmat27(i,j) );
        if ga(current_ind_nearby)<threshold && isempty(find(ind_max==current_ind_nearby)) && mask(current_ind_nearby)
          nn(i,1)=p1_repmat27(i,j);
          nn(i,2)=p2_repmat27(i,j);
          nn(i,3)=p3_repmat27(i,j);
          flag_this_voxel_is_corected=1;
        end
        j=j+1;
    end
    
    if ~flag_this_voxel_is_corected
        nn(i,1)=0;
        nn(i,2)=0;
        nn(i,3)=0;
    end
   
 end
        
