function [anat,anat_k] = gxanat_recon(filename,nii_write)

if nargin < 2
    nii_write = true;
end

dset = ismrmrd.Dataset(filename,'dataset');
hdr = ismrmrd.xml.deserialize(dset.readxml);

ImSize = [hdr.encoding.reconSpace.matrixSize.x hdr.encoding.reconSpace.matrixSize.y hdr.encoding.reconSpace.matrixSize.z];

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
else
    noise = [];
end
meas  = D.select(firstScan:D.getNumber);
clear D;

%% All data points should contribute to image

%Need to extract data from cell and put into an array
FID_Array = zeros(length(meas.data{1}),length(meas.data),size(meas.data{1},2));
Traj_Array = zeros(3,length(meas.data{1}),length(meas.data));

for i = 1:length(meas.data)
    FID_Array(:,i,:) = meas.data{i};
    Traj_Array(:,:,i) = meas.traj{i};
end

anat_k = FID_Array;
Traj_Array = Reconstruct.traj_delay_correction(Traj_Array,2);
%% Reconstruct Images
%Reshape to column vectors
Trajr = [reshape(Traj_Array(1,:,:),1,[])' reshape(Traj_Array(2,:,:),1,[])' reshape(Traj_Array(3,:,:),1,[])'];
Img = zeros([ImSize size(FID_Array,3)]);
for i = 1:size(FID_Array,3)
    tmp = FID_Array(:,:,i);
    tmpr = reshape(tmp,1,[])';
    Img(:,:,:,i) = Reconstruct.h1_recon(ImSize,tmpr,Trajr);
end
    
anat = sqrt(sum(Img.^2,4));
%%
if nii_write
    right_path = 1;
    [path1,~,~] = fileparts(filename);
    while right_path 
        [path2,~,~] = fileparts(path1);
        check_path = dir(path2);
        check_path = struct2cell(check_path);
        if sum(find(strcmp(check_path(1,:),'QC')))
            participant_folder = path1;
            right_path = 0;
        end
        path1 = path2;
    end
    
    folders = dir(participant_folder);
    folders = struct2cell(folders);
    getnames = folders(1,:);
    myfolderind = find(contains(getnames,'sub-'));
    bidsfolder = getnames{myfolderind};
    
    Subj_ID = bidsfolder(5:end);
    
    nii_name1 = ['sub-' Subj_ID '_anat'];
    
    writeanat = ReadData.mat2canon(abs(anat));
    niftiwrite(writeanat,fullfile(participant_folder,bidsfolder,'xegx',nii_name1),'Compressed',true);

end


