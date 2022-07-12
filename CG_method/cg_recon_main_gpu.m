% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-cg_recon_main_gpu-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description:
% -----------
% 
% Implementation of expanded signal model image reconstruction (Wilm et al., 2011).
% Image is reconstructed using k-space data, trajectory information, coil 
% sensitivity map, and B0 nonuniformity map. The forward model is then
% solved using Conjugate Gradient algorithm.
% Current implementation works for single-band excitation sequences (NO SMS).
% 
%   **GPU is required for image reconstruction**
%
% Inputs:
% ------
%
%    data:    structure with data for image reconstruction:
%
%       data.kloc:  coefficients for spherical harmonic basis function of trajectory (rad/m) [Nk,Norder,Ncontrast,Nslice]
% 
%       data.kdata: measured k-space samples for each coil and each acquisition [Nk,Ncoil,Ncontrast,Nslice]
% 
%       data.header:  header file of measured data in ISMRMRD format that includes:
%                   
%                 data.header.position: slice center
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
%    nIter: number of itterations for Conjugate Gradient (CG)
% 
%    vis:   turn on/off showing information about every iteration in CG
% 
% Outputs:
% -------
% 
%    im: reconstructed image [Nx,Ny,Nz,Naquisition,Niter]
% 
% Article: Feizollah and Tardif (2022)
% -------
%
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function im=cg_recon_main_gpu(data,nIter,vis)

N=data.header.hdr.encoding.reconSpace.matrixSize.x;
Slice=double(max(data.header.idx.slice)+1);

if(~vis)
    tic
end

im=zeros(N,N,Slice,nIter,size(data.kdata,3));
for n=1:Slice
    if(vis)
        tic
        fprintf("Slices: "+num2str(n)+"/"+num2str(Slice)+"...");
    end
    
    param=cg_prepare_recon_data_gpu(data,n);
    I=cg_solver(param,nIter,vis);  
    img=zeros(N*N,nIter+1);
    for k=1:param.Naq
        img(param.indx_im,:)=squeeze(gather(I(:,k,:)));
        im(:,:,n,:,k)=reshape(img(:,2:end),[N,N,nIter]);
    end
    if(vis)
        toc
    end
end
if(~vis)
    toc
end