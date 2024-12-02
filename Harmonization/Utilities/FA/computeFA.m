function [fa3d,md3d, All_eigs]=computeFA(uu, bb, S,mask, option)
[nx ny nz d] = size(S);
S_double =  reshape(S, nx*ny*nz, d);
uu2=uu(1:size(uu,1)/2,:);
mask_d = reshape(mask, nx*ny*nz, 1);
[dd, DD] = direct_1T(uu2, bb, S_double');

DD(isnan(DD))=0;
dd(isnan(dd))=0;

[fa, T_]=tensor2fa_vec(dd);
fa(isnan(fa))=0;


All_eigs=[];
MD = zeros(size(S_double,1 ),1);
if option.compute_eigvec
for iii=1:size(dd,2)
    if mask_d(iii)
    a=dd(1,iii)';
    b=dd(2,iii)';
    c=dd(3,iii)';
    d=dd(4,iii)';
    e=dd(5,iii)';
    f=dd(6,iii)';
    
    %syms a b c d e f
    T=[a b c;b d e; c e f];
    
    [VV,DD]=eig(T);
    
    MD(iii) = mean(diag(DD));
    %eig_vec=VV(:,3);
    %All_eigs=[All_eigs; eig_vec'];
    end
end

end

fa3d=single(reshape(fa,[nx,ny,nz]));
md3d=single(reshape(MD,[nx,ny,nz]));



end
    