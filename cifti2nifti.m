function cifti2nifti(image)
% Beta version to convert cifti images into nifti. this function should be
% extended to convert different types of information in the cifti structure
% ATM only converts parcelations stored in the dlabel file
% 
% Input: CIFTI image (dlabel)
%
% Output: NIFTI image in MNI space(2mm)
%_________________________________________________________________________
% Carlos Hernandez 2018

[dir,n]=spm_fileparts(image);
name=strsplit(n,'.');

% Read CIFTI
cii=ft_read_cifti(image);

% Conversion from mm to vox
vox=(cii.pos/2)+repmat([45 63 36],size(cii.pos,1),1);

% Write Labels
DATA=zeros(91,109,91);
for i=1:31870
    DATA(vox(i,1),vox(i,2),vox(i,3))=cii.indexmax(i); % MK changed 'brainstructure' to 'indexmax'
end

% MNI 2mm data structure
V = struct('fname',[dir 'N2C_' name{1} '.nii'],...
            'dim',[91,109,91],...
            'dt',[4 0],...
            'mat',[-2 0 0 92;0 2 0 -128;0 0 2 -74;0 0 0 1],...
            'pinfo',[1;0;352],...
            'descrip','MNI_2mm');
V = spm_create_vol(V);
spm_write_vol(V,DATA);