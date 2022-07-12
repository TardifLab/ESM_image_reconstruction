function H=cg_conc_func(x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates spherical harmonics basis functions.
%Input:
%   x: meshgridded values of x- [Nx*Nx,1]
%   y: meshgridded values of y- [Ny*Ny,1]
%   z: meshgridded values of z- [Nz*Nz,1]
% Output:
%   H: concomitant gradient basis functions- [4,Nimg^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h(1,:,:) = z.*z;
h(2,:,:) = x.*x + y.*y;
h(3,:,:) = x.*z;
h(4,:,:) = y.*z;

for k=1:4
    tmp=squeeze(h(k,:,:));
    H(k,:)=tmp(:);
end