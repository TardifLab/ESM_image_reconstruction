% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-cg_prepare_recon_data_gpu-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description:
% -----------
% 
% Prepares data for image reconstruction using Conjugate Gradient method.
%
% Inputs:
% ------
%
%    data:    structure containing data for image reconstruction:
%
%       data.kloc:  trajectory coefficients for spherical harmonic basis function (rad/m) [Nk,Norder,Ncontrast,Nslice]
% 
%       data.kdata: measured k-space samples [Nk,Ncoil,Ncontrast,Nslice]
% 
%       data.header:  header file of measured data in ISMRMRD format that includes:
%                   
%                 data.header.position: position of slice center
%                   
%                 data.header.hdr.encoding.reconSpace.matrixSize.x: matrix size
% 
%                 data.header.hdr.encoding.reconSpace.fieldOfView_mm.x: Field of view (m)
% 
%                 data.header.hdr.encoding.trajectory: type of trajectory
% 
%       data.sens:  coil sensitivity map [Nx,Ny,Nz,Ncoil]
% 
%       data.b0:    B0 nonuniformity map [Nx,Ny,Nz]
% 
%       data.mask:  mask for image reconstruction [Nx,Ny,Nz]
%
%   n:  number of slice to be reconstructed
% 
% Outputs:
% -------
% 
%    recon_data: structure with parameters for image recon.
%   
%       recon_data.h: spherical harmonics basis function [Norder,N*N]
%       recon_data.kloc: trajectory data [Nk,Norder,Ncontrast]
%       recon_data.kdatak-space raw data [Nk,Ncoil,Ncontrast]
%       recon_data.Np: total number of pixels (N*N)
%       recon_data.sens: coil sensitivity map [N*N,Ncoils]
%       recon_data.b0: B0 nonuniformity map [1,N*N]
%       recon_data.t: time vector of samples [Nk,1]
%       recon_data.Nc: number of coils (Ncoil)
%       recon_data.Nk: number of k-space samples (Nk)
%       recon_data.Nimg: matrix size (N)
%       recon_data.Naq: number of contrasts (Ncontrast)
% 
% Article: Feizollah and Tardif (2022)
% -------
%
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function recon_data=cg_prepare_recon_data_gpu(data,n)

[~,~,Naq,~]=size(data.kdata(:,:,:,n));
pos=data.header.position(:,n);
N=data.header.hdr.encoding.reconSpace.matrixSize.x;
fov=data.header.hdr.encoding.reconSpace.fieldOfView_mm.x*1e-3;

mask=data.mask(:,:,n);
mask=mask(:);
recon_data.indx_im=find(mask);

b0=data.b0(:,:,n);
b0=b0(:);
b0=b0(recon_data.indx_im);

sens=data.sens(:,:,n,:);
sens=reshape(sens,[N*N,size(sens,4)]);
sens=sens(recon_data.indx_im,:);

[X,Y,Z]=meshgrid((-.5:1/N:.5-1/N)*fov-pos(1),(-.5:1/N:.5-1/N)*fov+pos(2),-pos(3));
hh=sph_harmonics(X,-Y,Z);
conc=conc_func(X,-Y,Z)*0;
h=cat(1,hh(1:size(data.kloc,2)-4,:),conc);

recon_data.h=gpuArray(single(h(:,recon_data.indx_im)));
recon_data.kloc=gpuArray(complex(single(1i*data.kloc(:,:,:,n))));
recon_data.kdata=gpuArray(complex(single(data.kdata(:,:,:,n))));
recon_data.Np=gpuArray(size(b0,1));
recon_data.sens=gpuArray(complex(single(sens)));
recon_data.b0=gpuArray(single(b0(:).'));
recon_data.t=gpuArray(complex(single(1i*data.header.time)));
recon_data.Nc=gpuArray(size(data.kdata,2));
recon_data.Nk=gpuArray(size(recon_data.kloc,1));
recon_data.Nimg=gpuArray(N);
recon_data.fov=gpuArray(fov);
recon_data.Naq=gpuArray(Naq);

if(data.header.hdr.encoding.trajectory=="spiral"||data.header.hdr.encoding.trajectory=="sp")
    recon_data.dcf=gpuArray(cg_dcf_spiral(imag(recon_data.kloc))/N/N);
else
    recon_data.dcf=gpuArray(ones(length(recon_data.kloc),1)/N/N);
end
