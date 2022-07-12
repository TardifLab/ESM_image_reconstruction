% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-cg_maps_interpolate-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description:
% -----------
% 
% 3D interpolation of coil sensitivity, B0 nonuniformity and mask to match
% positions of final the scan.
%
% Inputs:
% ------
%
%    sens: coil sensitivity maps [Nx,Ny,Nz,Ncoils]
% 
%    b0: B0 map [Nx,Ny,Nz]
% 
%    mask: mask for image recon [Nx,Ny,Nz]
% 
%    header: ISMRMRD header with position information of the scan and maps
% 
% Outputs:
% -------
% 
%    sens_new: interpolated coil sensitivity map [Nx,Ny,Nz,Ncoils]
% 
%    b0_new:  interpolated B0 nonuniformity map [Nx,Ny,Nz]
% 
%    mask_new: interpolated mask for image recon [Nx,Ny,Nz]
% 
% Article: Feizollah and Tardif (2022)
% -------
%
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function [sens_new,b0_new,mask_new]=cg_maps_interpolate(sens,b0,mask,header)

[Ny,Nx,~]=size(b0);
fovx=header.maps.hdr.encoding.reconSpace.fieldOfView_mm.x*1e-3;
fovy=header.maps.hdr.encoding.reconSpace.fieldOfView_mm.y*1e-3;
if(length(header.maps.position)>size(b0,3))
    header.maps.position=header.maps.position(:,1:size(b0,3));
end

[z,I]=sort(single(header.maps.position(3,:)),'ascend');
y=single(header.maps.position(2,1));
x=single(header.maps.position(1,1));
x=linspace(x-fovx/2,x+fovx/2-fovx/Nx,Nx);
y=linspace(y-fovy/2,y+fovy/2-fovy/Ny,Ny);
[X,Y,Z]=meshgrid(x,y,z);

zv=single(header.position(3,:));
[Xv,Yv,Zv]=meshgrid(x,y,zv);
b0_new=interp3(X,Y,Z,b0(:,:,I),Xv,Yv,Zv,'spline');

[Ny,Nx,~,~]=size(sens);
y=single(header.maps.position(2,1));
x=single(header.maps.position(1,1));
x=linspace(x-fovx/2,x+fovx/2-fovx/Nx,Nx);
y=linspace(y-fovy/2,y+fovy/2-fovy/Ny,Ny);
[X,Y,Z]=meshgrid(x,y,z);

for k=1:size(sens,4)
    sens_new(:,:,:,k)=interp3(X,Y,Z,sens(:,:,I,k),Xv,Yv,Zv,'makima');
end

[Nx,Ny,~,~]=size(mask);
y=single(header.maps.position(2,1));
x=single(header.maps.position(1,1));
x=linspace(x-fovx/2,x+fovx/2-fovx/Nx,Nx);
y=linspace(y-fovy/2,y+fovy/2-fovy/Ny,Ny);
[X,Y,Z]=meshgrid(x,y,z);

mask_new=interp3(X,Y,Z,mask(:,:,I),Xv,Yv,Zv,'spline');
mask_new(mask_new~=0)=1;
