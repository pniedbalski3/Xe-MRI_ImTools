function [gas_k,dis_k,Gas_Traj_Array,Dis_Traj_Array] = get_raw_dixon(filename,delay_cor,nii_write)

if nargin < 2
    delay_cor = 0;
    nii_write = false;
end
if nargin < 3
    nii_write = false;
end

dset = ismrmrd.Dataset(filename,'dataset');
hdr = ismrmrd.xml.deserialize(dset.readxml);

ImSize = [hdr.encoding.reconSpace.matrixSize.z hdr.encoding.reconSpace.matrixSize.z hdr.encoding.reconSpace.matrixSize.z];

%% Read in all data
D = dset.readAcquisition();

%% Ignore noise scans
% TODO add a pre-whitening example
% Find the first non-noise scan
% This is how to check if a flag is set in the acquisition header
isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
firstScan = find(isNoise==0,1,'first');
if firstScan > 1
    noise = D.select(1:firstScan-1);
elseramp
    noise = [];
end
meas  = D.select(firstScan:D.getNumber);
clear D;

%% Get Gas Data
if nnz(meas.head.idx.set) > 0
    Gas_ind = zeros(max(meas.head.idx.set),length(meas.head.idx.set));
    for i = 1:max(meas.head.idx.set)
        Gas_ind(i,:) = meas.head.idx.contrast == 1 & meas.head.measurement_uid == 0 & meas.head.idx.set == i;
    end
else
    Gas_ind = meas.head.idx.contrast == 1 & meas.head.measurement_uid == 0;
end
for i = 1:size(Gas_ind,1)
    if nnz(Gas_ind(i,:)) ~=0
        Gas_FID(i,:) = meas.data(logical(Gas_ind(i,:)));
        Gas_Traj(i,:) = meas.traj(logical(Gas_ind(i,:)));
    end
end

%Need to extract data from cell and put into an array
Gas_FID_Array = zeros(length(Gas_FID{1}),length(Gas_FID),max([max(meas.head.idx.set) 1]));
Gas_Traj_Array = zeros(3,length(Gas_FID{1}),length(Gas_FID),max([max(meas.head.idx.set) 1]));

for j = 1:size(Gas_FID,1)
    for i = 1:length(Gas_FID)
        Gas_FID_Array(:,i,j) = Gas_FID{j,i};
        Gas_Traj_Array(:,:,i,j) = Gas_Traj{j,i};
    end
end

% %Kill the first 20 points to get to steady state (and avoid any first
% %projection weirdness)
% Gas_FID_Array(:,1:20) = [];
% Gas_Traj_Array(:,:,1:20) = [];

gas_k = Gas_FID_Array;

%% Get Dissolved Data

if nnz(meas.head.idx.set) > 0
    Dis_ind = zeros(max(meas.head.idx.set),length(meas.head.idx.set));
    for i = 1:max(meas.head.idx.set)
        Dis_ind(i,:) = meas.head.idx.contrast == 2 & meas.head.measurement_uid == 0 & meas.head.idx.set == i;
    end
else
    Dis_ind = meas.head.idx.contrast == 2 & meas.head.measurement_uid == 0;
end
for i = 1:size(Dis_ind,1)
    Dis_FID(i,:) = meas.data(logical(Dis_ind(i,:)));
    Dis_Traj(i,:) = meas.traj(logical(Dis_ind(i,:)));
end

%Need to extract data from cell and put into an array
Dis_FID_Array = zeros(length(Dis_FID{1}),length(Dis_FID),max([max(meas.head.idx.set) 1]));
Dis_Traj_Array = zeros(3,length(Dis_FID{1}),length(Dis_FID),max([max(meas.head.idx.set) 1]));

for j = 1:size(Dis_FID,1)
    for i = 1:length(Dis_FID)
        Dis_FID_Array(:,i,j) = Dis_FID{j,i};
        Dis_Traj_Array(:,:,i,j) = Dis_Traj{j,i};
    end
end

%Kill the first 20 points to get to steady state (and avoid any first
%projection weirdness)
% Dis_FID_Array(:,1:20) = [];
% Dis_Traj_Array(:,:,1:20) = [];
dis_k = Dis_FID_Array;

% if contains(hdr.acquisitionSystemInformation.institutionName,'Iowa')
%     holdDis_FID_Array = Dis_FID_Array;
%     holdDis_Traj_Array = Dis_Traj_Array;
%     Dis_FID_Array = (Gas_FID_Array);
%     Dis_Traj_Array = Gas_Traj_Array;
%     Gas_FID_Array = holdDis_FID_Array;
%     Gas_Traj_Array = holdDis_Traj_Array;
% end


%% Reconstruct Images
%Reshape to column vectors

% for i = 1:size(Gas_FID_Array,3)
%     Gas_FIDr = reshape(Gas_FID_Array(:,:,i),1,[])';
%     if nnz(Gas_FIDr) == 0
%         continue
%     end
%     Gas_Traj_Array(:,:,:,i) = Reconstruct.traj_delay_correction(Gas_Traj_Array(:,:,:,i),delay_cor);
%     Gas_Trajr = [reshape(Gas_Traj_Array(1,:,:,i),1,[])' reshape(Gas_Traj_Array(2,:,:,i),1,[])' reshape(Gas_Traj_Array(3,:,:,i),1,[])'];
% 
%     gas_hires(:,:,:,i) = Reconstruct.sharp_kern_recon(ImSize,Gas_FIDr,Gas_Trajr);
%     gas_lores(:,:,:,i) = Reconstruct.broad_kern_recon(ImSize,Gas_FIDr,Gas_Trajr);
% end
% 
% for i = 1:size(Dis_FID_Array,3)
%     Dis_FIDr = reshape(Dis_FID_Array(:,:,i),1,[])';
%     if nnz(Dis_FIDr) == 0
%         continue
%     end
%     Dis_Traj_Array(:,:,:,i) = Reconstruct.traj_delay_correction(Dis_Traj_Array(:,:,:,i),delay_cor);
%     Dis_Trajr = [reshape(Dis_Traj_Array(1,:,:,i),1,[])' reshape(Dis_Traj_Array(2,:,:,i),1,[])' reshape(Dis_Traj_Array(3,:,:,i),1,[])'];
% 
%     dis(:,:,:,i) = Reconstruct.broad_kern_recon(ImSize,Dis_FIDr,Dis_Trajr);
% end
% 
% figure('Name','Check Trajectory Delay')
% montage(abs(gas_hires(:,:,:,1)/max(abs(gas_hires(:)))));
% 
% %% Write Images to file if desired
% if nii_write
%     right_path = 1;
%     [path1,~,~] = fileparts(filename);
%     while right_path 
%         [path2,~,~] = fileparts(path1);
%         check_path = dir(path2);
%         check_path = struct2cell(check_path);
%         if sum(find(strcmp(check_path(1,:),'QC')))
%             participant_folder = path1;
%             right_path = 0;
%         end
%         path1 = path2;
%     end
% 
%     folders = dir(participant_folder);
%     folders = struct2cell(folders);
%     getnames = folders(1,:);
%     myfolderind = find(contains(getnames,'sub-'));
%     bidsfolder = getnames{myfolderind};
% 
%     Subj_ID = bidsfolder(5:end);
% 
%     nii_name1 = ['sub-' Subj_ID '_sgas'];
%     nii_name2 = ['sub-' Subj_ID '_bgas'];
%     nii_name3 = ['sub-' Subj_ID '_dis'];
% 
%     writesgas = ReadData.mat2canon(abs(gas_hires(:,:,:,1)));
%     writebgas = ReadData.mat2canon(abs(gas_lores(:,:,:,1)));
%     writedis = ReadData.mat2canon(abs(dis(:,:,:,1)));
% 
%     niftiwrite(writesgas,fullfile(participant_folder,bidsfolder,'xegx',nii_name1),'Compressed',true);
%     niftiwrite(writebgas,fullfile(participant_folder,bidsfolder,'xegx',nii_name2),'Compressed',true);
%     niftiwrite(writedis,fullfile(participant_folder,bidsfolder,'xegx',nii_name3),'Compressed',true);
% end
% 
