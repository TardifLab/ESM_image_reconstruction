% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-cg_solver-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description:
% -----------
% 
% Conjugate Gradient method to reconstruct images iteratively.
% 
%   **Requires GPU for calculations**
%
% Inputs:
% ------
%
%    recon_data: structure with data for image recon.
%   
%       recon_data.h: spherical harmonics basis function [Norder,N]
%       recon_data.kloc: trajectory data [Nk,Norder,Ncontrast]
%       recon_data.kdatak-space raw data [Nk,Ncoil,Ncontrast]
%       recon_data.Np: total number of pixels (N*N)
%       recon_data.sens: coil sensitivity map [N*N,Ncoils]
%       recon_data.b0: B0 nonuniformity map [N*N,1]
%       recon_data.t: time vector of samples [Nk,1]
%       recon_data.Nc: number of coils (Ncoil)
%       recon_data.Nk: number of k-space samples (Nk)
%       recon_data.Nimg: matrix size (N)
%       recon_data.Naq: number of contrasts (Ncontrast)
% 
%    nIter: number of iterations
% 
%    vis: turn on/off to show information for every step (true/false)
% 
% Outputs:
% -------
% 
%    b: reconstructed image for every iteration [N*N,Ncontrast,NIter]
% 
% Article: Feizollah and Tardif (2022)
% -------
%
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function b=cg_solver(recon_data,nIter,vis)
if(vis)
    fprintf('Reconstructing...')
end

alpha=0;
sen=recon_data.sens;
Csen=conj(recon_data.sens);

[L_block,N_block]=cg_encoding_loops_gpu(recon_data);

a=gpuArray(zeros(recon_data.Np,recon_data.Naq));
for k=1:N_block
    indx=((k-1)*L_block(k)+1:(k-1)*L_block(k)+L_block(k+1));
    for j=1:recon_data.Naq
        Evar=exp(recon_data.kloc(indx,:,j)*recon_data.h+recon_data.t(indx,:)*recon_data.b0);
        a(:,j)=a(:,j)+sum(Csen.*(Evar'*(recon_data.kdata(indx,:,j).*recon_data.dcf(indx))),2);
    end
end
p=a;
r(:,1)=p(:);

b=gpuArray(zeros(recon_data.Np,recon_data.Naq,nIter+1));
q=gpuArray(zeros(recon_data.Np,recon_data.Naq));

for i = 1:nIter

    delta=(r(:,i)'*r(:,i))/(a(:)'*a(:));
    q=0.*q;
    
    for k=1:N_block
        indx=((k-1)*L_block(k)+1:(k-1)*L_block(k)+L_block(k+1));
        for j=1:recon_data.Naq
            Evar=exp(recon_data.kloc(indx,:,j)*recon_data.h+recon_data.t(indx,:)*recon_data.b0);
            q(:,j)=q(:,j)+sum(Csen.*(Evar'*((Evar*(p(:,j).*sen)).*recon_data.dcf(indx))),2);
        end
    end
    q=q+alpha*p;
    
    b(:,:,i+1)=b(:,:,i)+((r(:,i)'*r(:,i))/(p(:)'*q(:)))*p;

    r(:,i+1)=r(:,i)-((r(:,i)'*r(:,i))/(p(:)'*q(:)))*q(:);
    p=reshape(r(:,i+1),[recon_data.Np,recon_data.Naq])+((r(:,i+1)'*r(:,i+1))/(r(:,i)'*r(:,i)))*p;
end
if(vis)
    fprintf("Itr "+num2str(i)+", "+"resid="+num2str(abs(delta))+", ");
end
