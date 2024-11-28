function [Y] = makespharms(u,L,flag)
% Compute Spherical Harmonics
%   MAKESPHARMS(u,L), where u is an N-by-3 matrix of N samples of the unit
% sphere, computes the spherical harmonics of degree L at the points given
% by u. The output is an N-by-2*L+1 matrix of the spherical harmonics.
%
%                           [Y] = makespharms(u,L)
% Input:
%              u - samples of the unit sphere
%              L - degree of spherical harmonics
% Output:
%              Y - matrix of spherical harmonics
%
% See also SPHERIHARMS, ICOSAHEDRON
%
% written by Oleg Michailovich, October 2007

n=size(u,1);
theta=acos(u(:,3));                 % deflection from z-axis: [0 pi]
varphi=atan2(u(:,2),u(:,1));
varphi=varphi+(varphi<0)*(2*pi);    % rotation in xy-plane: [0 2*pi]

if (L==0)
    Y=repmat(1/(sqrt(4*pi)),[n 1]);
else
    P=legendre(L,cos(theta))';
    P=[P(:,L+1:-1:2)*diag((-1).^(L:-1:1)) P];

    E=exp(1i*varphi);
    E=repmat(E,[1 2*L+1]).^repmat((-L:L),[n 1]); % repmat(-2 -1 0 1 2, [n 1])

    m=abs(-L:L);
    C=sqrt(((2*L+1)/(4*pi))*(factorial(L-m)./factorial(L+m)));
    C=repmat(C,[n 1]);
    Y=C.*P.*E;
    p=1:L;
    
    if nargin == 2
    % restrict basis: real, symmetric  (Jimi)
    Y_ = sqrt(2) * real(Y(:,1:L));  % m < 0 
    Y_ = [Y_ Y(:,L+1)];             % m = 0
    Y_ = [Y_  repmat((-1).^p,[n 1]) .* (sqrt(2) * imag(Y(:,L+2:end)))]; % m > 0
    Y = Y_;
    end
end