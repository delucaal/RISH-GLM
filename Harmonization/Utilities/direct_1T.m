function [d, D] = direct_1T(u, b, s)
  ux = u(:,1); uy = u(:,2); uz = u(:,3);
  if numel(b) == 1
    B = -b * [ux.^2  2*ux.*uy  2*ux.*uz  uy.^2  2*uy.*uz  uz.^2];
  else
    B = repmat(-b,[1 6]) .* [ux.^2  2*ux.*uy  2*ux.*uz  uy.^2  2*uy.*uz  uz.^2];
  end

  d = real(B \ log(s)); % ensure real since unconstrained
  
  if nargout == 2
    D = reshape(d([1 2 3 2 4 5 3 5 6],:), 3, 3, []);
  end
end
