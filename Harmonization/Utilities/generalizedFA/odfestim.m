function [Cs,S,Cw,W,Y,Q] = odfestim(S,u,L,lambda)

% Estimate ADC and ODF signals
%    ODFESTIM(S,u,L,...) approximates the ADC signal S measured at the 
% (unit) directions defined by u (N-by-3 matrix of samples of the unit
% sphere) by a linear combination of even-order spherical harmonics up
% to the order L inclusive. The function also computes the related ODF
% as well as its representation coefficients by the harmonics.
%
%               [Cs,S,Cw,W,Y,Q] = odfestim(S,u,L,lambda)
% Input:
%                 S - data ADC  
%                 u - N-by-3 matrix of unit directions
%                 L - maximum order of spherical harmonics
%            lambda - regularization parameter
% Output:
%                Cs - representation coefficient of the ADC
%                Se - estimated ADC
%                Cw - representation coefficient of the ODF
%                We - estimated ODF
%                 Y - model matrix of spherical harmonics
%                 Q - projection matrix, e.g. Cw=Q*S.
%
% written by Oleg Michailovich, November 2007

if (nargin<4),
    lambda=0.006;
end

Y=makespharms(u,0);
R=0;
for l=2:2:L,
    Y=[Y makespharms(u,l)];
    R=[R; repmat(l,[2*l+1 1])];
end
B=(R.*(R+1)).^2;

Cs=(Y'*Y+lambda.*diag(B))\(Y'*S);
S=single(real(Y*Cs));

if nargout > 2
n = size(Cs,2);    
P=R(:);
P(1)=1;
for k=2:length(P),
    %P(k)=((-1)^(P(k)/2))*(prod(1:2:P(k)-1)/prod(2:2:P(k)));
    P(k)=((-1)^(P(k)/2))*(prod(2:2:P(k))/prod(1:2:P(k)-1)); %inverse FRT
end


    Cw=2*pi*repmat(P,[1 n]).*Cs;
    W=real(Y*Cw);
    %W=W/sum(W);

    if (nargout>5),
        Q=diag(2*pi*P)*(inv(Y'*Y+lambda.*diag(B))*Y');
    end
end