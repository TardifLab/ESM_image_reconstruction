% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-cg_load_maps-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description:
% -----------
% 
% Loads coil sensitivity maps, B0 nonuniformity, and mask from a .mat file.
% The mat file should contain: 'sens_map', 'bField', 'mask'.
%
% Inputs:
% ------
%
%    map_adrs: location of ISMRMRD file that contains raw data, trajectory,
%               noise data and header
% 
% Outputs:
% -------
% 
%    sens_map: coil sensitivity map [Nx,Ny,Nz,Ncoils]
% 
%    bField:  B0 nonuniformity map [Nx,Ny,Nz]
% 
%    mask: mask for image recon [Nx,Ny,Nz]
% 
%    header: header information of maps containing position information
% 
% Article: Feizollah and Tardif (2022)
% -------
%
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function [sens_map,bField,mask,header]=cg_load_maps(map_adrs)

data=load(map_adrs);

names=fieldnames(rmfield(data,'header'));
if(names{1}=="bField")
    Nz=size(data.bField,3);
elseif(names{1}=="sens_map")
    Nz=size(data.sens_map,3);
elseif(names{1}=="mask")
    Nz=size(data.mask,3);
else
    Nz=1;
end
if(isfield(data,'bField'))
    bField=data.bField;
else
    warning('No B0 field found, set to zero!')
    bField=zeros(128,128,Nz);
end

if(isfield(data,'sens_map'))
    sens_map=data.sens_map;
else
    warning('No coil sensitivity found, set to one!')
    sens_map=ones(128,128,Nz);
end

if(isfield(data,'mask'))
    mask=data.mask;
else
    warning('No mask found, set to one!')
    mask=ones(128,128,Nz);
end

if(isfield(data,'header'))
    header=data.header;
end
