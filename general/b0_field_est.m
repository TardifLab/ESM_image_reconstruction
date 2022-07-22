% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-b0_field_est-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% Estimates B0 field map from phase images of a multi-echo GRE scan. It uses
% phase unrwraping over the echos per voxel, then does a linear fit. It is
% designed to estimate Bo map for every channel, but it can be set to 1 for
% a combined image.
%
% Inputs:
% -----
% 
%     pha: phase images of every channel from mGRE scan in (rad) [Nx,Ny,Nz,Necho,Ncoil]
% 
%     echo: echo time in (ms) [Necho,1]
% 
% Outputs:
% -------
%       
%     B0_map: estimated B0 map in (rad/s) [Nx,Ny,Nz,Ncoil]
% 
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function B0_map=b0_field_est(pha,echo)

[Nx,Ny,Nz,~,Nc]=size(pha);
dEcho=(echo-echo(1))*1e-3;
dEcho=[ones(length(dEcho),1) dEcho];

warning('off','all')
parfor c=1:Nc
    for k=1:Nz
        for j=1:Ny
            for i=1:Nx
%                 coef=[6,sum(dEcho);sum(dEcho),sum(dEcho.^2)]\...
%                     [sum(tmp_unw(i,j,k,:,c));sum(dEcho.*squeeze(tmp_unw(i,j,k,:,c)))];
%                 B0_map(i,j,k,c)=coef(2);
                b=dEcho\squeeze(unwrap(squeeze(pha(i,j,k,:,c))));
                B0_map(i,j,k,c)=b(2);
            end
        end
    end
end
warning('on','all')