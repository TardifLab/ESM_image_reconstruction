% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-cg_conc_func-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% Calculates concomitant field basis functions 
%
% Inputs:
% ------
%
%    x: meshgridded values of x- [Nx*Nx,1]
%    y: meshgridded values of y- [Ny*Ny,1]
%    z: meshgridded values of z- [Nz*Nz,1]
%
% Outputs:
% -------
% 
%    H: concomitant gradient basis functions [4,N*N]
%
% Article: Feizollah and Tardif (2022)
% -------
%
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function H=cg_conc_func(x,y,z)

h(1,:,:) = z.*z;
h(2,:,:) = x.*x + y.*y;
h(3,:,:) = x.*z;
h(4,:,:) = y.*z;

for k=1:4
    tmp=squeeze(h(k,:,:));
    H(k,:)=tmp(:);
end
