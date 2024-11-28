%rawS = 4D diffusion signal in a subject
%grad = gradient directions
%mask = 3D brain mask
function [ODFall]=odfReconstruction(DWI, grad, mask)  


% - Loading the reconstruction scheme (mesh where the ODF will be reconstructed)
%   ---------------------------------
% unit vectors on the sphere
disp(' ')
V = load('/home2/pnl/Downloads/hardi_tools_1.2/sampling_and_reconstruction_schemes/On_the_sphere/724_shell.txt');
% facets
display(['The reconstruction scheme contains ' num2str(size(V,1)) ' directions']);
F = load('/home2/pnl/Downloads/hardi_tools_1.2/sampling_and_reconstruction_schemes/On_the_sphere/724_sphere_facets.txt');
disp(' ');


%  --- Regularization parameter
Lambda = 0.001;


% --- real spherical harmonic reconstruction: parameters definition --- %
[Lmax Nmin] = obtain_Lmax(grad);
if Lmax >= 8
    Lmax = 8;
end

[basisG, thetaG, phiG] = construct_SH_basis (Lmax, grad, 2, 'real');
[basisV, thetaV, phiV] = construct_SH_basis (Lmax, V, 2, 'real');
K = []; Laplac2 = [];
for L=0:2:Lmax
    for m=-L:L
        % factor1 = ((-1)^(L/2))*doublefactorial(L+1)./( (L+1)*doublefactorial(L) );
        Pnm = legendre(L,0); factor1 = Pnm(1);
        K = [K; factor1];
        Laplac2 = [Laplac2; (L^2)*(L + 1)^2];
    end
end
Laplac = diag(Laplac2);

% Creating the kernel for the reconstruction
A = recon_matrix(basisG,Laplac,Lambda);
% where A = inv(Ymatrix'*Ymatrix + Lambda*Laplacian)*Ymatrix';
% and Ymatrix = basisG.
% The matrix, A, is the same for all the voxels


ODFall=zeros(size(DWI, 1),size(DWI, 2), size(DWI,3), length(V),'single');
for xx=1:size(DWI, 1)
    for yy=1:size(DWI, 2)
        for zz=1:size(DWI,3)
            if(mask(xx,yy,zz))
                S=double(squeeze(DWI(xx,yy,zz,:)));
                coeff = A*S; 
                ss = coeff.*K;
                ODF = basisV*ss;
                ODF = ODF - min(ODF);
                ODFall(xx,yy,zz,:) = ODF/sum(ODF); % normalization
            end
        end
    end
   xx
end   
  

% set(gcf, 'color', 'white');
% screen_size = get(0, 'ScreenSize');
% f1 = figure(1); hold on;
% set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
% cameramenu;
% % -----
% Origen = [0 0 0.25]; % Voxel coordinates [x y z] where the ODF will be displayed
% angle = 0:.01:2*pi;
% R = ones(1,length(angle));
% R1 = R.*f2/f1;
% subplot(1,2,1)
% polarm(angle,R,'k'); hold on; polarm(angle, R1,'k');
% % the ODF is scaled only for the display
% plot_ODF(ODF./max(ODF), V, F, Origen); 
% title('\fontsize{14} QBI from the signal without noise', 'FontWeight','bold')
% axis equal; axis off;
% 
% subplot(1,2,2)

end
