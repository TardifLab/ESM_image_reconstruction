% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-example_cg_recon-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description:
% -----------
% 
% Reads raw data of a scan in ISMRMRD format, coil sensitivity maps, B0
% nonuniformity, and mask, Then reconstructs an image using expanded 
% signal model.
%
% Inputs:
% ------
% 
% Outputs:
% -------
% 
% Article: Feizollah and Tardif (2022)
% -------
%
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

data_adrs='';
map_adrs='';

fprintf('Loading data...')
[data.kdata,kloc,header,data.noise]=cg_ismrmrd_sort(data_adrs);
fprintf('Finished\n')

N=header.hdr.encoding.reconSpace.matrixSize.x;

[b0_init,sens_init,mask_init,header.maps]=cg_load_maps(map_adrs);

[sens,b0,mask]=cg_maps_interpolate(b0_init,sens_init,mask_init,header);

mask=imresize(mask,[N,N]);
sens=imresize(sens,[N,N]);
b0=imresize(b0,[N,N]);

data.kloc=kloc;
data.header=header;
data.sens=sens;
data.b0=b0;
data.mask=mask;

im=cg_recon_main_gpu(data,12,false);
