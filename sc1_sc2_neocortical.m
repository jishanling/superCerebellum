function varargout=sc1_sc2_neocortical(what,varargin)
% Analysis of cortical data 
% This is modernized from the analysis in sc1_sc2_imana, in that it uses
% the new FS_LR template

baseDir    = '/Volumes/MotorControl/data/super_cerebellum_new';
wbDir      = fullfile(baseDir,'sc1','surfaceWB'); 
fsDir      = fullfile(baseDir,'sc1','surfaceFreesurfer'); 
atlasDir   = '~/Data/Atlas_templates/standard_mesh';       
anatomicalDir = fullfile(baseDir,'sc1','anatomicals'); 

subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11',...
    's12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24',...
    's25','s26','s27','s28','s29','s30','s31'};
returnSubjs=[2,3,4,6,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31];

switch(what) 
    case 'SURF:recon_all'                    % STEP 2.2: Calls recon_all
        % STUDY 1 ONLY
        % Calls recon-all, which performs, all of the
        % FreeSurfer cortical reconstruction process
        % example: sc1_imana('surf_freesurfer',1)
        sn=varargin{1}; % subjNum
        
        for i=sn
            freesurfer_reconall(fullfile(baseDir,'sc1','surfaceFreesurfer'),subj_name{i},fullfile(anatomicalDir,subj_name{i},['anatomical.nii']));
        end
    case 'SURF:wb_resample' 
        % This reslice from the individual surface into the the fs_lr standard mesh - This replaces        
        % calls to freesurfer_registerXhem, freesurfer_mapicosahedron_xhem, & caret_importfreesurfer. 
        % It requires connectome wb to be installed. 
        sn=returnSubjs; 
        vararginoptions(varargin,{'sn'}); 
        for i=sn
            fprintf('reslicing %d\n',i); 
            surf_resliceFS2WB(subj_name{i},fsDir,atlasDir,wbDir); 
        end; 
        
        
  case 'PREP:cortex:surface_betas'         % STEP 11.5: Map betas and ResMS (.nii) onto surface (.metric)
        % Run FREESURFER before this step!
        % map volume images to metric file and save them in individual
        % surface folder
        % example: sc1_sc2_imana('map_con_surf',2,1,4,'betas')
        sn   = varargin{1}; % subjNum
        study=varargin{2};
        glm  = varargin{3}; % glmNum
        contrast = varargin{4}; % 'beta'or 'ResMS'
        
        glmDir =[studyDir{study} sprintf('/GLM_firstlevel_%d',glm)];dircheck(glmDir);
        
        subjs=length(sn);
        
        vararginoptions({varargin{5:end}},{'atlas','regSide'});
        
        for s=1:subjs,
            glmSubjDir=fullfile(glmDir, subj_name{sn(s)});
            for h=regSide,
                caretSDir = fullfile(studyDir{study}, caretDir,[atlasA,subj_name{sn(s)}],hemName{h}); dircheck(caretSDir);
                white=fullfile(studyDir{1},caretDir,[atlasA,subj_name{sn(s)}],hemName{h},[hem{h} '.WHITE.coord']);
                pial=fullfile(studyDir{1},caretDir,[atlasA,subj_name{sn(s)}],hemName{h},[hem{h} '.PIAL.coord']);
                
                C1=caret_load(white);
                C2=caret_load(pial);
                fileList = [];
                outfile  = [];
                
                filenames = dir(fullfile(glmSubjDir,sprintf('*%s*',contrast)));
                outfile = sprintf('%s_glm%d_%s_cortex_%s.metric',subj_name{sn(s)},glm,contrast,hem{h});
                
                for t = 1:length(filenames),
                    fileList{t} = {filenames(t).name};
                end
                
                for f=1:length(fileList),
                    images(f)=spm_vol(fullfile(glmSubjDir,fileList{f}));
                end;
                M=caret_vol2surf_own(C1.data,C2.data,images);
                caret_save(fullfile(caretSDir,outfile),M);
                
                fprintf('%s map to surface for %s:%s \n',contrast,subj_name{sn(s)},hemName{h});
            end;
        end
end; 