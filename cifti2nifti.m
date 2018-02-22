function cifti2nifti(image)
% Beta version to covert cifti image into nifti
% This function should be extended to convert different cifti file types
% ATM only converts parcellations stored in the dlabel.nii file
% 
% Input: CIFTI image (dlabel)
%
% Output: NIFTI image in MNI space(2mm)
%_________________________________________________________________________
% Carlos Hernandez 2018
% Amended by Maedbh King (22/02/18). Implemented formula v=inv(A)*w to go from world to voxel
% coordinates. Also - now codes labels and not brainstructure

[dir,n]=spm_fileparts(image);
name=strsplit(n,'.');

% Read CIFTI
cii=ft_read_cifti(image);

numVox=size(cii.indexmax,1); 

% Conversion from mm to vox
%vox=(cii.pos/2)+repmat([47 63 36],size(cii.pos,1),1);

% this should give the same answer as above solution ? 
for v=1:numVox,
    tmp=inv(cii.transform)*[cii.pos(v,:),1]'; 
    vox(v,:)=tmp'; 
end

% Write Labels
DATA=zeros(cii.dim);
for i=1:numVox,
    DATA(vox(i,1),vox(i,2),vox(i,3))=cii.indexmax(i); % MK changed this to always code labels (not brainstructure)
end

% MNI 2mm data structure
V = struct('fname',[dir 'N2C_' name{1} '.nii'],... % this should be C2N ?
            'dim',cii.dim,...
            'dt',[4 0],...
            'mat',[-2 0 0 92;0 2 0 -128;0 0 2 -74;0 0 0 1],...
            'pinfo',[1;0;352],...
            'descrip','MNI_2mm');
V = spm_create_vol(V);
spm_write_vol(V,DATA);