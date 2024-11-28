% ODF (2D): [nx*ny*nz,d]
% gFA (3D): [nx, ny, nz ] 
function gFA = generalizedFA(ODF, sz) 
nx  = sz(1) ; ny = sz(2); nz = sz(3); d = sz(4);
mu = mean(ODF,2);
gFA = real(sqrt((d*sum((ODF - repmat(mu, 1, d)).^2,2))./((d-1)*sum(ODF.^2,2))));
gFA =  reshape(gFA, nx, ny, nz);
end
    
