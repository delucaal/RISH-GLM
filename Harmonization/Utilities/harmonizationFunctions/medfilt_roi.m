function ScaleLn = medfilt_roi(scaleL, id)

[nx, ny, nz] = size(scaleL);

[xi, yi, zi] = ind2sub([nx ny nz], id);

ScaleLn = scaleL;
ws = 2;
iz  = ones(length(zi),1)*ws;
iz (zi>=nz-1) = 0;
iz (zi<=2 ) =0;

iy  = ones(length(yi),1)*ws;
iy (yi>=ny-1) = 0;
iy (yi<=2) = 0;

ix  = ones(length(xi),1)*ws;
ix (xi>=nx-1 ) =0;
ix (xi<=2) = 0;

for l=1:length(id)
    A=(scaleL(xi(l)-ix(l): xi(l)+ix(l), yi(l)-iy(l): yi(l)+iy(l), zi(l)-iz(l):zi(l)+iz(l)));
    
    szA = size(A,1)*size(A,2)*size(A,3);
    ScaleLn(xi(l), yi(l), zi(l)) = squeeze(median(reshape(A, szA, 1),1));
end
end