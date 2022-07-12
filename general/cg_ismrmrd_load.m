% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-cg_ismrmrd_load-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description:
% -----------
% 
% Reads raw data of a scan in ISMRMRD format containing k-space data,
% trajectory, header, and noise.
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
%    kdata: raw k-space data [Nk,Nsegment,Ncoils]
% 
%    trajectory:  trajectory data. Coefficients of spherical harmonics basis
%           functions [Norder,Nk,Nsegment]
% 
%    header: header information of ISMRMRD format
% 
%    nosie: acquired noise data [Nnoise,Ncoil]
% 
% Article: Feizollah and Tardif (2022)
% -------
%
% Sajjad Feizollah, July 2022
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function [kdata,trajectory,header,noise]=cg_ismrmrd_load(data_adrs)


if(nargin<1)
    disp('Select ISMRMRD file.....')
    [baseName, folder] = uigetfile('./*.h5','Data File');
    data_adrs = fullfile(folder, baseName);
end

mrd = ismrmrd.Dataset(data_adrs);
header.hdr=ismrmrd.xml.deserialize(mrd.readxml);
acq = mrd.readAcquisition();
L= mrd.getNumberOfAcquisitions();

noise_idx=find(acq.head.flags==262144|acq.head.flags==262208);
noise=[];
if(~isempty(noise_idx))
    for k=length(noise_idx):length(noise_idx)
        noise=acq.data{1,noise_idx(k)};
    end
end

data_start_idx=max(noise_idx);
if(isempty(noise))
    data_start_idx=0;
end
kdata=single(zeros(size(acq.data{1,data_start_idx+1},1),size(acq.data,2)-data_start_idx,32));
trajectory=single(zeros(size(acq.traj{1,data_start_idx+1},1),size(acq.data{1,data_start_idx+1},1),size(acq.data,2)-data_start_idx));

for k=1:L-data_start_idx
    kdata(:,k,:)=acq.data{1,k+data_start_idx};
    trajectory(:,:,k)=acq.traj{1,k+data_start_idx};
end

header.version = (acq.head.version(data_start_idx+1:end));
header.flags = (acq.head.flags(data_start_idx+1:end));
header.measurement_uid = (acq.head.measurement_uid(data_start_idx+1:end));
header.scan_counter = (acq.head.scan_counter(data_start_idx+1:end));
header.acquisition_time_stamp = (acq.head.acquisition_time_stamp(data_start_idx+1:end));
header.physiology_time_stamp = (acq.head.physiology_time_stamp(:,data_start_idx+1:end));
header.number_of_samples = (acq.head.number_of_samples(data_start_idx+1:end));
header.available_channels = (acq.head.available_channels(data_start_idx+1:end));
header.active_channels = (acq.head.active_channels(data_start_idx+1:end));
header.channel_mask = (acq.head.channel_mask(:,data_start_idx+1:end));
header.discard_pre = (acq.head.discard_pre(data_start_idx+1:end));
header.discard_post = (acq.head.discard_post(data_start_idx+1:end));
header.center_sample = (acq.head.center_sample(data_start_idx+1:end));
header.encoding_space_ref = (acq.head.encoding_space_ref(data_start_idx+1:end));
header.trajectory_dimensions = (acq.head.trajectory_dimensions(data_start_idx+1:end));
header.sample_time_us = (acq.head.sample_time_us(data_start_idx+1:end));
header.position = (acq.head.position(:,data_start_idx+1:end));
header.read_dir = (acq.head.read_dir(:,data_start_idx+1:end));
header.phase_dir = (acq.head.phase_dir(:,data_start_idx+1:end));
header.slice_dir = (acq.head.slice_dir(:,data_start_idx+1:end));
header.patient_table_position = (acq.head.patient_table_position(:,data_start_idx+1:end));
header.idx.kspace_encode_step_1 = (acq.head.idx.kspace_encode_step_1(data_start_idx+1:end));
header.idx.kspace_encode_step_2 = (acq.head.idx.kspace_encode_step_2(data_start_idx+1:end));
header.idx.average = (acq.head.idx.average(data_start_idx+1:end));
header.idx.slice = (acq.head.idx.slice(data_start_idx+1:end));
header.idx.contrast = (acq.head.idx.contrast(data_start_idx+1:end));
header.idx.phase = (acq.head.idx.phase(data_start_idx+1:end));
header.idx.repetition = (acq.head.idx.repetition(data_start_idx+1:end));
header.idx.set = (acq.head.idx.set(data_start_idx+1:end));
header.idx.segment =  (acq.head.idx.segment(data_start_idx+1:end));
header.idx.user = (acq.head.idx.user(:,data_start_idx+1:end));
header.user_int = (acq.head.user_int(:,data_start_idx+1:end));
header.user_float = (acq.head.user_float(:,data_start_idx+1:end));
