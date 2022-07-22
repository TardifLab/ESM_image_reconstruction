% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-cg_b0_field_calc-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
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
%     echo: echo time vector in (ms) [Necho,1]
% 
% Outputs:
% -------
%       
%     B0_map: estimated B0 map in (rad/s) [Nx,Ny,Nz,Ncoil]
% 
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function B0_map=b0_field_calc(pha,echo)

[Nx,Ny,Nz,~,Nc]=size(pha);
dEcho=(echo-echo(1))*1e-3;
dEcho=[ones(length(dEcho),1) dEcho];

parfor c=1:Nc
    for k=1:Nz
        for j=1:Ny
            for i=1:Nx
                b=dEcho\squeeze(unwrap(squeeze(pha(i,j,k,:,c))));
                B0_map(i,j,k,c)=b(2);
            end
        end
    end
end
