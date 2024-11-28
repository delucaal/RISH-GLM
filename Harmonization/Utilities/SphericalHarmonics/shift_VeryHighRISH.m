
function nnn=shift_VeryHighRISH(img,conn,mask,threshold)

if conn==18
    th=1.6;
elseif conn==26
    th=3;
elseif conn==6
    th=1.2;
end

[nx ny nz]= size(img);
nnn=reshape(img,nx,ny, nz);
ind_max=find(nnn>threshold);

        
[p1,p2,p3]=ind2sub([nx ny nz],ind_max);


[x,y,z]=meshgrid(-1:1:1);
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

 for i=1:size(p1_repmat27,1)
    flag_this_voxel_is_corected=0;
    j=1;
    while j<=size(p1_repmat27,2) && ~flag_this_voxel_is_corected
        current_ind_nearby=sub2ind([nx ny nz], p1_repmat27(i,j) ,p2_repmat27(i,j) ,p3_repmat27(i,j) );
        if nnn(current_ind_nearby)<threshold && isempty(find(ind_max==current_ind_nearby)) && mask(current_ind_nearby)==1
          nnn(ind_max(i))=nnn(current_ind_nearby);
          flag_this_voxel_is_corected=1;
        end
        j=j+1;
    end
    if ~flag_this_voxel_is_corected
        nnn(ind_max(i))=threshold;
    end

 end
        
