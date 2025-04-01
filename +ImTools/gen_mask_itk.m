function mask = gen_mask_itk(image,start_mask,other_image)

if nargin < 2
    start_mask = zeros(size(image));
end
if nargin < 3
    other_image = [];
end

itk_path = ImTools.get_itk_path();
ITKSNAP_Path = ['"C:\Program Files\' itk_path '\bin\ITK-SNAP.exe"'];

tmp_file = 'C:\Users\pniedbalski\OneDrive - University of Kansas Medical Center\Documents\local_tmp';

niftiwrite(image,fullfile(tmp_file,'tmp_im_4_segment'),'Compressed',true);

niftiwrite(start_mask,fullfile(tmp_file,'tmp_mask_4_segment'),'Compressed',true);

vent_full_path = fullfile(tmp_file,'tmp_im_4_segment.nii.gz');
maskpath = fullfile(tmp_file,'tmp_mask_4_segment.nii.gz');

if ~isempty(other_image)
    niftiwrite(other_image,fullfile(tmp_file,'tmp_other_im_4_segment'),'Compressed',true);
    other_path = fullfile(tmp_file,'tmp_other_im_4_segment.nii.gz');
    mycommand = [ITKSNAP_Path ' -g "' vent_full_path '" -o "' other_path '" -s "' maskpath '"'];
else
    mycommand = [ITKSNAP_Path ' -g "' vent_full_path '" -s "' maskpath '"'];
end

system(mycommand);

mask = niftiread(maskpath);