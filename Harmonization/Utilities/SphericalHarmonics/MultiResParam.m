function [Y B R P u conn] = MultiResParam(res, grads, L)

if(nargin == 2)
    L = 6;
end
conn = [];
R = 0;
if(res > 0)
    u = icosahedron(res);
    fcs=convhulln(u); % convhull is a matlab function
    %u = sortrows(u);
    Y  = makespharms(u, 0);
    for l = 2:2:L
      Y = [Y makespharms(u,l)];
      R = [R; repmat(l,[2*l+1 1])];
    end
    B=(R.*(R+1)).^2;
    P=R(:);
    P(1)=1;
    for k=2:length(P),
        P(k)=((-1)^(P(k)/2))*(prod(1:2:P(k)-1)/prod(2:2:P(k)));
    end
    if(nargout == 6)
     conn = FindConnectivity(fcs,size(u,1));
    end
else
    u = grads;
    Y  = makespharms(u, 0);
    for l = 2:2:L
      Y = [Y makespharms(u,l)];
      R = [R; repmat(l,[2*l+1 1])];
    end
    B=(R.*(R+1)).^2;
    P=R(:);
    P(1)=1;
    for k=2:length(P),
        P(k)=((-1)^(P(k)/2))*(prod(1:2:P(k)-1)/prod(2:2:P(k)));
    end
    
end
