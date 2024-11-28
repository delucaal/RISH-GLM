% S = 4D diffusion signal in a subject
% grad = gradient directions
% mask = 3D brain mask
% Single shell
function [IMG]=SHCoefficients(S, grad, mask, SH_level_docomposition)

IMG.shs_same_level=[1 1; 2 6; 7 15 ; 16 28; 29 45];
S = single(S).*single(repmat(mask,[1,1,1,size(S,4)])); 
[IMG.nx, IMG.ny, IMG.nz, d] = size(S);
L = IMG.nx * IMG.ny * IMG.nz;
S = reshape(S,[L d]);


lambda=10^-6;
[IMG.Y, B, R, P] = MultiResParam(0, grad, SH_level_docomposition);

%get the spherical harmonics coefficients
IMG.Cs=double(IMG.Y'*IMG.Y+lambda.*diag(B))\(IMG.Y'*S');
IMG.mask=mask;

% taking too much space 
%IMG.Estimated_S_using_Cs = ((real(single(IMG.Y)*single(IMG.Cs(:,1:lenght(IMG.Cs)))');

IMG.rishImgs = zeros(IMG.nx, IMG.ny, IMG.nz, 5);
for i=1:SH_level_docomposition/2+1
    
    temp_rish1D= IMG.Cs(IMG.shs_same_level(i,1):IMG.shs_same_level(i,2) ,:);
    IMG.rishImgs(:,:,:,i)=reshape(temp_rish1D(:).*double(mask(:)), IMG.nx, IMG.ny, IMG.nz);
end

end
  
