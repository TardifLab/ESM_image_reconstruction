% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-cg_encoding_loops_gpu-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description:
% -----------
% 
% Determines number of blocks and their length for faster calculations and 
% to avoid memory issues.
% 
% **IMPORTANT: max_numel is selected based on GPU memory to
%   balance between maximum number of matrix size and the load on GPU. It is
%   different for every GPU, so please modify to achieve maximum
%   performance. It is when the maximum power is drawn by the GPU. The default value is
%   tested for RTX 3090 and Tesla V100.**
% 
%
% Inputs:
% ------
%
%    params.Np: number of all pixels (N*N)
% 
%    params.Nk: number of samples in trajectory
% 
% Outputs:
% -------
% 
%    L_blocks: length of every block for calculation
% 
%    N_block: total number of blocks
% 
% Article: Feizollah and Tardif (2022)
% -------
%
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function [L_blocks,N_block]=cg_encoding_loops_gpu(params)

max_numel=0.15e9;

Np=params.Np;
Nk=params.Nk;

if(Np*Nk<=max_numel)
    L_blocks=[0,Nk];
    N_block=1;
else
    L_block=round(max_numel/Np);
    N_block=floor(Nk/L_block)+1;
    L_blocks=[0,L_block*ones(1,N_block-1),rem(Nk,L_block)];
end
