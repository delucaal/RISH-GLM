
function [corrected_cs,ind_max]=threshold_SH_Coeffs(IMG,option)
threshold = option.threshold;
if threshold < 10
    mask = option.skull;
    conn=6;
else
    mask = option.mask - option.skull;
    conn=6;
end


ind_max=[];
if conn==18
    th=1.6;
elseif conn==26
    th=3;
elseif conn==6
    th=1.2;
end
tem= (sum((IMG.Cs(: ,:).^2),1));
tem=tem(:);

corrected_cs=IMG.Cs.*double(repmat(mask(:)',size(IMG.Cs,1),1));
nnn=reshape(tem,[IMG.nx IMG.ny IMG.nz]).*double(mask);%double(IMG.mask);%nz ny nx
ind_max=find(nnn>threshold);
if  strcmp(option.process,'after') && threshold >10
     
        
        [p1,p2,p3]=ind2sub([IMG.nx IMG.ny IMG.nz],ind_max);
        
        
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
        p1_repmat27(p1_repmat27>IMG.nx)=IMG.nx;
        
        p2_repmat27(p2_repmat27<1)=1;
        p2_repmat27(p2_repmat27>IMG.ny)=IMG.ny;
        
        p3_repmat27(p3_repmat27<1)=1;
        p3_repmat27(p3_repmat27>IMG.nz)=IMG.nz;
      
         for i=1:size(p1_repmat27,1)
            flag_this_voxel_is_corected=0;
            j=1;
            while j<=size(p1_repmat27,2) && ~flag_this_voxel_is_corected
                current_ind_nearby=sub2ind([IMG.nx IMG.ny IMG.nz], p1_repmat27(i,j) ,p2_repmat27(i,j) ,p3_repmat27(i,j) );
                if nnn(current_ind_nearby)<threshold && isempty(find(ind_max==current_ind_nearby)) && mask(current_ind_nearby)==1
                   corrected_cs(:,ind_max(i))=IMG.Cs(:,current_ind_nearby);
                  flag_this_voxel_is_corected=1;
                  %  break;
%                 elseif  isempty(find(ind_max==current_ind_nearby)) && mask(current_ind_nearby)==0
%                     corrected_cs(:,ind_max(i))=0;
%                     flag_this_voxel_is_corected=1;
%                   %  break;
                end
                j=j+1;
            end
            if ~flag_this_voxel_is_corected
                corrected_cs(:,ind_max(i))=0;
            end
                
         end
        
else 
corrected_cs(:,ind_max)=0; %=IMG.Cs;

end
