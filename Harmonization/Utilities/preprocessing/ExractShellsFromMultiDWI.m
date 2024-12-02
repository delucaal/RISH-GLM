function ExractShellsFromMultiDWI(inputdwi, outputdwi, bvecs, bvals, binc)

% read data
nii = load_untouch_nii(inputdwi);
S = nii.img;
fid = fopen(bvecs,'r'); u = fscanf(fid, '%f',[inf ]); fclose(fid);
indT = strfind(inputdwi, '.');
type = inputdwi(indT(1):end);

%check if vectors lie on the column format or row format
[status,output]=system(['wc -l < ', bvecs]);
if status==0 & str2num(output)==3
    usz=length(u)/3;
    uhat=zeros(usz,3);
    for i=1:usz
        uhat(i,1)=u(i);
        uhat(i,2)=u(i+usz);
        uhat(i,3)=u(i+2*usz);
    end
    u=uhat;
else
    usz=length(u)/3;
    uhat=zeros(usz,3);
    ccc=1;
    for i=1:usz
        uhat(i,1)=u(ccc,1);
        uhat(i,2)=u(ccc+1,1);
        uhat(i,3)=u(ccc+2,1);
        ccc=ccc+3;
    end
    u=uhat;
end

% read b vals
fid = fopen(bvals,'r'); b = fscanf(fid, '%f'); fclose(fid);

% extract baseline
u_sum = (sum(u.^2,2)); 
ind_0 =  u_sum <= 0.05 | b<500;
s0 = S(:,:,:,ind_0);
disp (['b0 is extracted ...']);

% We assume the shell increment by binc and smallest shell is binc
bThreshold = binc;


% extract first shell
ind_s =  b>=bThreshold-binc/2 & b<bThreshold+binc/2;
Ssingle = S(:,:,:,ind_s);
while ~isempty(Ssingle)
    
    % extract single shell
    uSingle = u(ind_s, :);
    bSingle = cat(1, zeros(size(s0,4), 1), repmat(bThreshold, size(uSingle,1),1));
    uSingle = cat(1, repmat([0 0 0], size(s0,4), 1), uSingle);
   
    
    % save single shell dwi
    nii.img = cat(4, s0, Ssingle);
    nii.hdr.dime.dim(5) = size(nii.img,4);
    save_untouch_nii(nii, [outputdwi '_b', num2str(bThreshold), type]);

    % save single shell bvec
    fid = fopen([outputdwi '_b', num2str(bThreshold), '.bvec'],'w'); 
    for i=1:length(uSingle)
     fprintf(fid, '%f %f %f\n',uSingle(i,1), uSingle(i,2), uSingle(i,3));
    end
    fclose(fid);
    
    % save single bval
    fid = fopen([outputdwi '_b', num2str(bThreshold), '.bval'],'w'); fprintf(fid, '%f\n', bSingle);
    fclose(fid);

    disp (['b' num2str(bThreshold), ' is extracted ...']);

    bThreshold = bThreshold+binc;
    ind_s =  b>=bThreshold-binc/2 & b<bThreshold+binc/2;
    Ssingle = S(:,:,:,ind_s);

end

end