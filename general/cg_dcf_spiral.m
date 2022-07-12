% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-cg_dcf_spiral-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% Quick estimation for density compensation factor (DCF) for spiral trajectories.
%
% Inputs:
% ------
%
%    kloc: k-space trajectory points in (rad/m) [Nk,2]
%
% Outputs:
% -------
% 
%    dcf: estimated DCF [Nk,1]
%
% Article: Feizollah and Tardif (2022)
% -------
%
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function dcf=cg_dcf_spiral(kloc)

if size(kloc,2)>3
    kloc=kloc(:,2:3);
end

theta=atan2(kloc(:,2),kloc(:,1));
theta=unwrap(theta);
dcf=(theta(2:end)-theta(1:end-1));
dcf=[dcf;dcf(end)].*sqrt(kloc(:,1).^2+kloc(:,2).^2);
