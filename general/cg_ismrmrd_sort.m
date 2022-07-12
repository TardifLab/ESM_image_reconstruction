% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-cg_ismrmrd_sort-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description:
% -----------
% 
% Reads raw data of a scan in ISMRMRD format containing k-space data,
% trajectory, header, and noise, then sorts the data for image recon.
%
% Inputs:
% ------
%
%    data_adrs: location of ISMRMRD file that contains raw data, trajectory,
%               noise data and header
% 
% Outputs:
% -------
% 
%    kdata: raw k-space data [Nk,Ncoils,Ncontrast,Nslice]
% 
%    kloc:  trajectory data. Coefficients of spherical harmonics basis
%           functions [Nk,Norder,Ncontrast,Nslice]
% 
%    header_new: header information of ISMRMRD format
% 
%    nosie: acquired noise data
% 
% Article: Feizollah and Tardif (2022)
% -------
%
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function [kdata,kloc,header_new,noise]=cg_ismrmrd_sort(data_adrs)


if(nargin<1)
    disp('Select ISMRMRD file.....')
    [baseName, folder] = uigetfile('./*.h5','Data File');
    data_adrs = fullfile(folder, baseName);
end

[rawdata,traj,header,noise]=load_ismrmrd(data_adrs);
Nseg=header.hdr.encoding.echoTrainLength;
Norder=size(traj,1);
Ncoil=double(size(rawdata,3));
Nsample=double(size(rawdata,1));
Nz=double(max(header.idx.slice+1));
Naq=double(size(rawdata,2)/Nseg/Nz);

rawdata_new=single(zeros(Nsample*Nseg,Ncoil,Naq*Nz));
traj_new=single(zeros(Norder,Nsample*Nseg,Naq*Nz));
header_new=header;
header_new.position=[];

traj=traj(1:Norder,:,:);
for k=1:Naq*Nz
    rawdata_new(:,:,k)=reshape(rawdata(:,(k-1)*Nseg+1:Nseg*k,:),[Nsample*Nseg,Ncoil]);
    traj_new(:,:,k)=reshape(traj(1:Norder,:,(k-1)*Nseg+1:Nseg*k),[Norder,Nsample*Nseg]);
end
header.position=header.position(:,1:Nseg:end);
header.idx.slice=header.idx.slice(1:Nseg:end);

clear rawdata traj
kdata=single(zeros(Nsample*Nseg,Ncoil,Naq,Nz));
kloc=single(zeros(Norder,Nsample*Nseg,Naq,Nz));

aq_idx=repmat(1:Naq,[Nz,1]);
aq_idx=aq_idx(:);

for k=1:Naq*Nz
    kdata(:,:,aq_idx(k),header.idx.slice(k)+1)=rawdata_new(:,:,k);
    kloc(:,:,aq_idx(k),header.idx.slice(k)+1)=traj_new(:,:,k);
    header_new.position(:,header.idx.slice(k)+1,aq_idx(k))=header.position(:,k);

end

kloc=permute(kloc,[2,1,3,4]);
header_new.position=header_new.position(:,:,1);

if(header.hdr.encoding.trajectory=="epi")
    [~,indx]=min(abs(kloc(:,3,1,1)));
    header_new.seg_center=round(indx/Nsample);
    dt=double(header.sample_time_us(1));
    t_echo_spacing=double(header.hdr.sequenceParameters.echo_spacing*1000);
    time=repmat(linspace(0,Nsample*dt-dt,Nsample)',[1,Nseg])+...
        repmat(0:t_echo_spacing:t_echo_spacing*Nseg-t_echo_spacing,[Nsample,1]);
    time=time-((header_new.seg_center-1)*Nsample*dt+Nsample*dt/2);
    header_new.time=time(:)*1e-6;
    header_new.symetric_k_indx=Nsample*(header_new.seg_center*2-1);
    
elseif(header.hdr.encoding.trajectory=="spiral")
    header_new.time=1e-6*linspace(0,size(kdata,1)*header.sample_time_us(1)-header.sample_time_us(1),size(kdata,1))';
end


[~,I]=sort(header_new.position(3,:),'descend');
header_new.position=header_new.position(:,I);
kdata=kdata(:,:,:,I);
kloc=kloc(:,:,:,I);

