function conn = FindConnectivity(fcs,N)

conn(N).elem = [];
for i = 1:N
    a1 = find(fcs(:,1) == i);
    a2 = find(fcs(:,2) == i);
    a3 = find(fcs(:,3) == i);
 
    u = [flat(fcs(a1,2:3))' flat(fcs(a2,1))' flat(fcs(a2,3))' flat(fcs(a3,1:2))'];
    conn(i).elem = unique(u);
end

end

function s = flat(S)
   s = S(:);
end
