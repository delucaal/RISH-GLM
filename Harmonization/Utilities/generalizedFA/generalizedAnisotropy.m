function ga3D = generalizedAnisotropy(S)
[nx ny nz d ] = size(S);
rms = sqrt(mean(S.^2, ndims(S)));
ga = std(S, 0, ndims(S)) ./ (rms + eps);
ga3D = reshape(ga, nx, ny, nz);
end
