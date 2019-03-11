function varargout=sc1_sc2_neocortical(what,varargin)
% Analysis of cortical data 
% This is modernized from the analysis in sc1_sc2_imana, in that it uses
% the new FS_LR template

baseDir    = '/Volumes/MotorControl/data/super_cerebellum_new';
wbDir      = fullfile(baseDir,'sc1','surfaceWB'); 
fsDir      = fullfile(baseDir,'sc1','surfaceFreesurfer'); 
atlasDir   = '~/Data/Atlas_templates/standard_mesh';       
anatomicalDir = fullfile(baseDir,'sc1','anatomicals'); 
studyDir  = {'sc1','sc2'}; 
Hem       = {'L','R'}; 
hemname   = {'CortexLeft','CortexRight'}; 
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
        
        
    case 'SURF:map_con'         % STEP 11.5: Map con / ResMS (.nii) onto surface (.gifti)
        sn    = returnSubjs;     % subjNum
        study = [1 2];
        glm  = 4;     % glmNum
        hemis = [1 2];
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt')); 
        vararginoptions(varargin,{'sn','study','glm','what','hemis'}); 
        
        for s=sn 
            for h=hemis
                surfDir = fullfile(wbDir, subj_name{s}); 
                white=fullfile(surfDir,sprintf('%s.%s.white.164k.surf.gii',subj_name{s},Hem{h})); 
                pial=fullfile(surfDir,sprintf('%s.%s.pial.164k.surf.gii',subj_name{s},Hem{h})); 
                C1=gifti(white);
                C2=gifti(pial);
                
                for st = study
                    glmDir =fullfile(baseDir,studyDir{st},sprintf('GLM_firstlevel_%d',glm),subj_name{s});
                    T=getrow(D,D.StudyNum==st); 
                    filenames={}; 
                    for i=1:length(T.condNames); 
                        filenames{i} = fullfile(glmDir,sprintf('con_%s-rest.nii',T.condNames{i})); 
                    end; 
                    outfile = fullfile(surfDir,sprintf('%s.%s.%s.con.func.gii',subj_name{s},Hem{h},studyDir{st}));
                
                    G=surf_vol2surf(C1.vertices,C2.vertices,filenames,'column_names',T.condNames,'anatomicalStruct',hemname{h});
                    save(G,outfile); 
                
                    fprintf('mapped %s %s %s \n',subj_name{s},Hem{h},studyDir{st});
                end; 
            end;
        end
end; 