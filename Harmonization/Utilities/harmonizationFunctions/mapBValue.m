% mapBValue:
% S: diffusion signal  (S)
% b_o: old b value
% b_n: new b value
% Sn: output
function Sn=mapBValue(S, s0, b_o, b_n)
[sx, sy, sz, d] = size(S);
S = S./repmat(s0, [1, 1, 1, d]);
Sr = reshape(S, sx*sy*sz,d);
Sr(Sr<=0) = 0;
Sr(Sr>1) = 1;
D = real(-log(single(Sr))./repmat(b_o', sx*sy*sz,1));
s0(s0<eps) = eps; 
Sn = reshape(exp(-b_n*D), sx, sy, sz, d).*repmat(s0, [1, 1, 1, d]);
end


