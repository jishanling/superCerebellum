function varargout=sc1_sc2_HCP(what,varargin)

% Directories
baseDir_orig          = '/Users/maedbhking/Documents/Cerebellum_Cognition';
baseDir          = '/Volumes/Seagate Backup Plus Drive';
% baseDir            = '/Volumes/MotorControl/data/super_cerebellum_new';
% baseDir          = '/Users/jdiedrichsen/Data/super_cerebellum_new';

% atlasDir='/Users/maedbhking/Documents/Atlas_templates/';

studyDir{1}     =fullfile(baseDir_orig,'sc1');
studyDir{2}     =fullfile(baseDir_orig,'sc2');
HCPDir          =fullfile(baseDir,'hcp');
suitDir         ='/suit';
caretDir        ='/surfaceCaret';
encodeDir       ='/encoding';
contrastDir     ='/contrasts';
anatDir         ='/anatomicals';

HCP_badSubjs=[122620];

% ijk coordinates
loc_AC = {[-94,-131,-114],... %s100307
    [-97,-131,-111],... %s100408
    [-93,-133,-114],...  %s101107
    [-93,-128,-114],... %s101309
    [-99,-136,-121],...%s101915
    [-94,-131,-112],... %s103111
    [-100,-133,-115],...%s103414
    [-97,-142,-119],... %s103818
    [-94,-127,-109],... %s105014
    [-94,-149,-132],... %s105115
    [-94,-141,-115],... %s106016
    [-101,-143,-126],... %s108828
    [-99,-131,-113],... %s110411
    [-99,-140,-124],... %s111312
    [-90,-127,-104],... %s111716
    [-95,-125,-105],... %s113619
    [-94,-126,-109],... %s113922
    [-92,-141,-124],... %s114419
    [-100,-147,-129],... %s115320
    [-100,-137,-115],... %s116524
    [-92,-143,-124],... %s117122
    [-94,-133,-114],... %s118528
    [-99,-129,-111],... %s118730
    [-99,-136,-107],...%s118932
    [-98,-132,-115],... %s120111
    [-94,-130,-106],... %s122317
    [-94,-128,-111],... %s122620
    [-91,-137,-115],... %s123117
    [-100,-130,-113],... %s123925
    [-94,-153,-127],... %s124422
    [-96,-137,-119],... %s125525
    [-99,-130,-108],... %s126325
    [-94,-137,-113],... %s127630
    [-91,-139,-123],... %s127933
    [-98,-140,-122],... %s128127
    [-100,-136,-114],... %s128632
    [-98,-133,-107],...%s129028
    [-92,-140,-125],... %s130013
    [-97,-133,-111],... %s130316
    [-102,-131,-115],...%s131217
    [-96,-145,-120],... %s131722
    [-100,-135,-120],... %s133019
    [-92,-130,-112],... %s133928
    [-105,-140,-123],...%s135225
    [-94,-126,-110],... %s135932
    [-103,-159,-137],... %s136833
    [-96,-139,-120],... %s138534
    [-98,-136,-115],... %s139637
    [-95,-147,-124],... %s140925
    [-94,-130,-104],... %s144832
    [-97,-131,-109],... %s146432
    [-104,-142,-116],...%s147737
    [-101,-132,-116],... %s148335
    [-101,-141,-113],...%s148840
    [-97,-133,-116],... %s149337
    [-91,-128,-108],... %s149539
    [-96,-144,-122],... %s149741
    [-96,-135,-117],... %s151223
    [-100,-149,-127],... %s151526
    [-98,-131,-111],... %s151627
    [-100,-133,-114],...%s153025
    [-95,-134,-112],... %s154734
    [-99,-143,-127],... %s156637
    [-96,-134,-116],... %s159340
    [-99,-135,-111],... %s160123
    [-92,-135,-117],... %s161731
    [-97,-145,-124],... %s162733
    [-101,-136,-125],... %s163129
    [-102,-132,-109],...%s176542
    [-95,-144,-127],... %s178950
    [-95,-138,-113],... %s188347
    [-93,-143,-122],... %s189450
    [-99,-129,-114],... %s190031
    [-95,-138,-116],... %s192540
    [-101,-128,-110],... %s196750
    [-100,-132,-114],... %s198451
    [-101,-131,-116],...%s199655
    [-97,-133,-111],... %s201111
    [-94,-129,-109],... %s208226
    [-95,-142,-117],... %s211417
    [-94,-136,-122],... %s211720
    [-103,-134,-112],...%s212318
    [-93,-140,-123],... %s214423
    [-99,-140,-120],... %s221319
    [-97,-127,-108],... %s239944
    [-96,-139,-118],... %s245333
    [-93,-137,-110],... %s280739
    [-100,-135,-113],... %s298051
    [-95,-130,-109],... %s366446
    [-95,-138,-121],... %s397760
    [-98,-133,-117],... %s414229
    [-98,-133,-112],... %s499566
    [-101,-140,-123],... %s654754
    [-95,-131,-108],... %s672756
    [-92,-125,-108],... %s751348
    [-91,-135,-111],... %s756055
    [-100,-138,-122],... %s792564
    [-102,-130,-113],... %s856766
    [-100,-127,-111],... %s857263
    [-99,-130,-108],... %s899885
    };

taskNames_new = {'EM_FACES','EM_SHAPES','GA_PUNISH','GA_REWARD','LA_MATH','LA_STORY','MO_CUE',...
    'MO_LF','MO_LH','MO_RF','MO_RH','MO_T','WM_TOOL','WM_PLACE','WM_0BK','WM_2BK',...
    'WM_BODY','WM_FACE','RE_MATCH','RE_REL','SO_RANDOM','SO_TOM','rest'};

taskNames={'EM_FACES','EM_SHAPES','GA_PUNISH','GA_REWARD','LA_MATH','LA_STORY','MO_CUE',...
    'MO_LF','MO_LH','MO_RF','MO_RH','MO_T','WM_TOOL','WM_PLACE','WM_0BK','WM_2BK',...
    'WM_BODY','WM_FACE','RE_MATCH','RE_REL','SO_RANDOM','SO_TOM'};


switch what
    
    case 'ANAT:load_files'
        anatRawDir=fullfile(HCPDir,'anatomicals_raw'); % where orig. HCP data are stored
        anats=dir(fullfile(anatRawDir,sprintf('*.zip*')));
        
        cd(anatRawDir);
        for i=1:length(anats),
            
            HCP_subjs=str2double(anats(i).name(1:6)); % subjNum
            
            if exist(fullfile(anatRawDir,sprintf('%s',anats(i).name)))
                % unzip
                unzip(fullfile(anatRawDir,sprintf('%s',anats(i).name)))
                outName1=sprintf('%d_3T_T1w_MPR1.nii.gz',HCP_subjs); opt1=fullfile(anatRawDir,'scans','103_T1w','NIFTI',outName1);
                opt2=fullfile(anatRawDir,'scans','105_T1w','NIFTI',outName1);
                outName2=sprintf('%d_3T_T1w_MPR2.nii.gz',HCP_subjs); opt3=fullfile(anatRawDir,'scans','105_T1w','NIFTI',outName2);
                opt4=fullfile(anatRawDir,'scans','107_T1w','NIFTI',outName1);
                opt5=fullfile(anatRawDir,'scans','106_T1w','NIFTI',outName1);
                opt6=fullfile(anatRawDir,'scans','505_T1w','NIFTI',outName1);
                opt7=fullfile(anatRawDir,'scans','104_T1w','NIFTI',outName1);
                opt8=fullfile(anatRawDir,'scans','503_T1w','NIFTI',outName1);
                opt9=fullfile(anatRawDir,'scans','110_T1w','NIFTI',outName1);
                
                tmp=strfind(outName1,'.gz');
                % get .nii file and move
                if exist(opt1),
                    gunzip(opt1)
                    niiFile=outName1(1:tmp-1);
                    mkdir(fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs)));
                    movefile(fullfile(anatRawDir,'scans','103_T1w','NIFTI',niiFile),fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs),'anatomical_raw.nii'))
                elseif exist(opt2)
                    gunzip(opt2)
                    niiFile=outName1(1:tmp-1);
                    mkdir(fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs)));
                    movefile(fullfile(anatRawDir,'scans','105_T1w','NIFTI',niiFile),fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs),'anatomical_raw.nii'))
                elseif exist(opt3)
                    gunzip(opt3)
                    niiFile=outName2(1:tmp-1);
                    mkdir(fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs)));
                    movefile(fullfile(anatRawDir,'scans','105_T1w','NIFTI',niiFile),fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs),'anatomical_raw.nii'))
                elseif exist(opt4),
                    gunzip(opt4)
                    niiFile=outName1(1:tmp-1);
                    mkdir(fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs)));
                    movefile(fullfile(anatRawDir,'scans','107_T1w','NIFTI',niiFile),fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs),'anatomical_raw.nii'))
                elseif exist(opt5),
                    gunzip(opt5)
                    niiFile=outName1(1:tmp-1);
                    mkdir(fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs)));
                    movefile(fullfile(anatRawDir,'scans','106_T1w','NIFTI',niiFile),fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs),'anatomical_raw.nii'))
                elseif exist(opt6),
                    gunzip(opt6)
                    niiFile=outName1(1:tmp-1);
                    mkdir(fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs)));
                    movefile(fullfile(anatRawDir,'scans','505_T1w','NIFTI',niiFile),fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs),'anatomical_raw.nii'))
                elseif exist(opt7),
                    gunzip(opt7)
                    niiFile=outName1(1:tmp-1);
                    mkdir(fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs)));
                    movefile(fullfile(anatRawDir,'scans','104_T1w','NIFTI',niiFile),fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs),'anatomical_raw.nii'))
                elseif exist(opt8),
                    gunzip(opt8)
                    niiFile=outName1(1:tmp-1);
                    mkdir(fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs)));
                    movefile(fullfile(anatRawDir,'scans','503_T1w','NIFTI',niiFile),fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs),'anatomical_raw.nii'))
                elseif exist(opt9),
                    gunzip(opt9)
                    niiFile=outName1(1:tmp-1);
                    mkdir(fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs)));
                    movefile(fullfile(anatRawDir,'scans','110_T1w','NIFTI',niiFile),fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs),'anatomical_raw.nii'))
                else
                    fprintf('subj %d not found \n',HCP_subjs)
                end
                if exist(fullfile(anatRawDir,'scans')),
                    rmdir(fullfile(anatRawDir,'scans'),'s')
                    fprintf('subj %d done \n',HCP_subjs)
                end
            else
                fprintf('subj %d not found \n',HCP_subjs)
            end
        end
    case 'ANAT:resample'
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        subjs=length(HCP_subjs);
        
        for s=1:subjs,
            
            % Downsample anat images to 1 mm ^3
            source    = fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs(s)),['anatomical_raw','.nii']);
            newFile   = fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs(s)),['anatomical_resample_raw','.nii']);
            spmj_resample(source,newFile,.7);
        end
    case 'ANAT:reslice_LPI'                  % STEP 1.2: Reslice anatomical image within LPI coordinate systems
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        subjs=length(HCP_subjs);
        
        for s=1:subjs,
            
            % (1) Reslice anatomical image to set it within LPI co-ordinate frames
            source  = fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs(s)),['anatomical_resample_raw','.nii']);
            dest    = fullfile(HCPDir,'anatomicals',sprintf('s%d',HCP_subjs(s)),['anatomical','.nii']);
            spmj_reslice_LPI(source,'name', dest);
            
            % (2) In the resliced image, set translation to zero
            V               = spm_vol(dest);
            dat             = spm_read_vols(V);
            V.mat(1:3,4)    = [0 0 0];
            spm_write_vol(V,dat);
            display 'Manually retrieve the location of the anterior commissure (x,y,z) before continuing'
        end
    case 'ANAT:centre_AC'                    % STEP 1.3: Re-centre AC
        % Set origin of anatomical to anterior commissure (must provide
        % coordinates in section (4)).
        % example: sc1_imana('ANAT:centre_AC',1)
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        subjs=length(HCP_subjs);
        
        for s=1:subjs,
            img    = fullfile(HCPDir,anatDir,sprintf('s%d',HCP_subjs(s)),['anatomical','.nii']);
            V               = spm_vol(img);
            dat             = spm_read_vols(V);
            V.mat(1:3,4)    = loc_AC{s};
            spm_write_vol(V,dat);
            fprintf('Done for s%d \n',HCP_subjs(s))
        end
        
    case 'FUNC:load_files'
        % unzip files
        fileNames=dir(fullfile(HCPDir,'*zip'));
        
        for i=1:length(fileNames),
            unzip(fullfile(HCPDir,fileNames(i).name))
            fprintf('%s unzipped \n',fileNames(i).name)
        end
        
    case 'SUIT:run_all'
        sn=varargin{1};
        
        sc1_sc2_ICB('SUIT:isolate_segment',sn)
        sc1_sc2_ICB('SUIT:make_maskImage',sn)
        sc1_sc2_ICB('SUIT:corr_cereb_cortex_mask',sn)
        sc1_sc2_ICB('SUIT:normalise_dartel',sn,'grey')
        sc1_sc2_ICB('SUIT:make_mask',sn,'grey')
        %         sc1_sc2_ICB('SUIT:reslice',sn)
    case 'SUIT:isolate_segment'              % STEP 9.2:Segment cerebellum into grey and white matter
        %        spm fmri
        
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        for s=1:length(HCP_subjs),
            suitSubjDir = fullfile(HCPDir,suitDir,'anatomicals',sprintf('s%d',HCP_subjs(s)));dircheck(suitSubjDir);
            source=fullfile(HCPDir,anatDir,sprintf('s%d',HCP_subjs(s)),'anatomical.nii');
            dest=fullfile(suitSubjDir,'anatomical.nii');
            copyfile(source,dest);
            cd(fullfile(suitSubjDir));
            suit_isolate_seg({fullfile(suitSubjDir,'anatomical.nii')},'keeptempfiles',1);
        end
    case 'SUIT:make_maskImage'               % STEP 3.7:Make mask images (noskull and grey_only)
        % Make maskImage meanepi
        % example: sc1_sc2_imana('FUNC:make_maskImage',1)
        
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        for s=1:length(HCP_subjs),
            
            % get example func image and mask
            nam{1}  = fullfile(HCPDir,'contrasts',sprintf('mask_EM_%d.nii',HCP_subjs(s))); % contrast image
            spm_imcalc(nam, fullfile(HCPDir,'contrasts',sprintf('mask_gray_%d.nii',HCP_subjs(s))), 'i1~=0')
            
            fprintf('subj %d done \n',HCP_subjs(s));
        end
    case 'SUIT:corr_cereb_cortex_mask'       % STEP 9.4:
        suitAnatDir=fullfile(HCPDir,'suit','anatomicals');
        
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        for s=1:length(HCP_subjs),
            
            cortexGrey= fullfile(suitAnatDir,sprintf('s%d',HCP_subjs(s)),'c3anatomical.nii');
            cerebGrey = fullfile(suitAnatDir,sprintf('s%d',HCP_subjs(s)),'c1anatomical.nii');
            bufferVox = fullfile(suitAnatDir,sprintf('s%d',HCP_subjs(s)),'buffer_voxels.nii');
            
            % isolate overlapping voxels
            spm_imcalc({cortexGrey,cerebGrey},bufferVox,'(i1.*i2)')
            
            % mask buffer
            spm_imcalc({bufferVox},bufferVox,'i1>0')
            
            cerebGrey2 = fullfile(suitAnatDir,sprintf('s%d',HCP_subjs(s)),'cereb_prob_corr_grey.nii');
            cortexGrey2= fullfile(suitAnatDir,sprintf('s%d',HCP_subjs(s)),'cortical_mask_grey_corr.nii');
            
            % remove buffer from cerebellum
            spm_imcalc({cerebGrey,bufferVox},cerebGrey2,'i1-i2')
            
            % remove buffer from cortex
            spm_imcalc({cortexGrey,bufferVox},cortexGrey2,'i1-i2')
        end
    case 'SUIT:normalise_dartel'             % STEP 9.5: Normalise the cerebellum into the SUIT template.
        % Normalise an individual cerebellum into the SUIT atlas template
        % Dartel normalises the tissue segmentation maps produced by suit_isolate
        % to the SUIT template
        % !! Make sure that you're choosing the correct isolation mask
        % (corr OR corr1 OR corr2 etc)!!
        % if you are running multiple subjs - change to 'job.subjND(s)."'
        % example: sc1_sc2_imana('SUIT:normalise_dartel',1,'grey')
        type=varargin{1}; % 'grey' or 'whole' cerebellar mask
        
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        for s=1:length(HCP_subjs),
            cd(fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s))));
            job.subjND.gray      = {'c_anatomical_seg1.nii'};
            job.subjND.white     = {'c_anatomical_seg2.nii'};
            switch type,
                case 'grey'
                    job.subjND.isolation= {'cereb_prob_corr_grey.nii'};
                case 'whole'
                    job.subjND.isolation= {'cereb_prob_corr.nii'};
            end
            suit_normalize_dartel(job);
        end
        
        % 'spm_dartel_warp' code was changed to look in the working
        % directory for 'u_a_anatomical_segment1.nii' file - previously it
        % was giving a 'file2mat' error because it mistakenly believed that
        % this file had been created
    case 'SUIT:make_mask'                    % STEP 9.7: Make cerebellar mask using SUIT
        type=varargin{1}; % 'grey' or 'whole'
        
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        for s=1:length(HCP_subjs),
            mask = fullfile(HCPDir,'contrasts',sprintf('mask_EM_%d.nii',HCP_subjs(s))); % mask for functional image
            switch type
                case 'grey'
                    suit  = fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s)),'cereb_prob_corr_grey.nii'); % cerebellar mask grey (corrected)
                    omask = fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s)),'maskbrainSUITGrey.nii'); % output mask image - grey matter
                case 'whole'
                    suit  = fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s)),'cereb_prob_corr.nii'); % cerebellar mask (corrected)
                    omask = fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s)),'maskbrainSUIT.nii'); % output mask image
            end
            cd(fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s))));
            spm_imcalc({mask,suit},omask,'i1>0 & i2>0.7',{});
        end
    case 'SUIT:reslice'                      % STEP 9.8: Reslice the contrast images from first-level GLM
        % Reslices the functional data (betas, contrast images or ResMS)
        % from the first-level GLM using deformation from
        % 'suit_normalise_dartel'.
        % example: sc1_sc2_imana('SUIT:reslice',1,1,4,'betas','cereb_prob_corr_grey')
        % make sure that you reslice into 2mm^3 resolution
        
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        T=[];
        for s=1:length(HCP_subjs),
            for t=1:length(taskNames),
                source=fullfile(HCPDir,'contrasts',sprintf('%s_%d.nii',taskNames{t},HCP_subjs(s))); % images to be resliced
                if exist(source),
                    job.subj.affineTr = {fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s)),'Affine_c_anatomical_seg1.mat')};
                    job.subj.flowfield= {fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s)),'u_a_c_anatomical_seg1.nii')};
                    job.subj.mask     = {fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s)),'cereb_prob_corr_grey.nii')};
                    job.vox           = [2 2 2];
                    cd(fullfile(HCPDir,'contrasts'))
                    job.subj.resample = {source};
                    dircheck(fullfile(HCPDir,'suit','contrasts'));
                    %             cd(fullfile(HCPDir,'suit','contrasts',HCP_subjs{sn(s)}));
                    suit_reslice_dartel(job);
                    temp=dir('*wd*');
                    for f=1:length(temp),
                        movefile(fullfile(HCPDir,'contrasts',temp(f).name),fullfile(HCPDir,'suit','contrasts',temp(f).name));
                    end
                else
                    fprintf('%s does not exist for %d',taskNames{t},sprintf('s%d',HCP_subjs(s)));
                    T.subj=HCP_subjs(s);
                    T.missingCon=taskNames{t};
                end
                
                fprintf('%s contrast has been resliced into suit space for %s \n\n',taskNames{t},sprintf('s%d',HCP_subjs(s)))
            end
        end
        
        varargout={T};
        keyboard;
        
    case 'PREP:HCP_subjs'
        % get all possible task conditions
        idx=1;
        conName=dir(fullfile(HCPDir,'anatomicals'));
        for c=4:length(conName),
            HCP_subjs(idx,1)=str2double(conName(c).name(2:end)); % assuming that the subjID is always 6 digits
            idx=idx+1;
        end
        
        varargout={HCP_subjs};
    case 'PREP:HCP_info'
        sn=varargin{1};
        % organise across subjs
        S=[];
        for s=1:length(sn),
            conName=dir(fullfile(HCPDir,'contrasts',HCP_subjs{sn(s)},'*sess*'));
            for c=1:length(conName),
                condName=conName(c).name(8:end-4);
                T.condNum(c,1)=D.condNum(strcmp(D.condNames,condName));
                T.condNames{c,1}=condName;
                tmp=str2double(conName(c).name(5:6));
                T.sessNum(c,1)=tmp;
            end
            T.SN=repmat(sn(s),length(T.condNames),1);
            S=addstruct(S,T);
            clear T conName tmp
        end
        
        varargout={S};
    case 'PREP:avrgMask_cereb'               % STEP 11.3:
        step=varargin{1}; % 'reslice' or 'mask'
        % don't include s01 and s04 in the 'mask' step
        
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        switch step,
            case 'reslice'
                for s=1:length(HCP_subjs),
                    cd(fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s))))
                    % normalise cerebellar grey into suit
                    job.subj.affineTr = {fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s)),'Affine_c_anatomical_seg1.mat')};
                    job.subj.flowfield= {fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s)),'u_a_c_anatomical_seg1.nii')};
                    job.subj.mask     = {fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s)),'cereb_prob_corr_grey.nii')};
                    job.vox           = [2 2 2];
                    job.subj.resample = {'c1anatomical.nii'};
                    suit_reslice_dartel(job);
                end
            case 'mask'
                for s=1:length(HCP_subjs),
                    nam{s}=fullfile(HCPDir,'suit','anatomicals',sprintf('s%d',HCP_subjs(s)),'wdc1anatomical.nii');
                end
                opt.dmtx = 1;
                cd(fullfile(HCPDir,'suit','anatomicals'));
                spm_imcalc(nam,'cerebellarGreySUIT.nii','mean(X)',opt);
                fprintf('averaged cerebellar grey mask in SUIT space has been computed \n')
        end
    case 'PREP:cereb:voxels'                 % STEP 11.6: Get UW cerebellar data (voxels)
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        % Load over all grey matter mask
        V=spm_vol(fullfile(HCPDir,'suit','anatomicals','cerebellarGreySUIT.nii'));
        
        X=spm_read_vols(V);
        grey_threshold = 0.1; % grey matter threshold
        volIndx=find(X>grey_threshold);
        [i,j,k]= ind2sub(size(X),volIndx');
        
        P=size(volIndx,1); % number of voxels
        
        T=[];
        for s=1:length(HCP_subjs),
            
            for c=1:length(taskNames),
                
                nam{1}=fullfile(HCPDir,'suit','contrasts',sprintf('wd%s_%d.nii',taskNames{c},HCP_subjs(s)));
                
                if exist(nam{1}),
                    Vi=spm_vol(nam{1});
                    C1=spm_sample_vol(Vi,i,j,k,0);
                    
                    S.data(c,:)=C1;
                    S.condName{c,1}=taskNames{c};
                    S.condNum(c,1)=c;
                    S.SN(c,1)=HCP_subjs(s);
                end
                fprintf('subj %d done for %s contrast \n',HCP_subjs(s),taskNames{c})
            end
            % add intercept (rest)
            S.data=[S.data; zeros(1,P)];
            S.condName=[S.condName;'rest'];
            S.condNum=[S.condNum;c+1];
            S.SN=[S.SN;HCP_subjs(s)];
            T=addstruct(T,S);
            clear S
        end
        outName=fullfile(HCPDir,'suit','contrasts','cereb_avrgDataStruct.mat');
        save(outName,'T','volIndx','V');
        
    case 'ACTIVITY:SNR'
        
        load(fullfile(HCPDir,'suit','contrasts','cereb_avrgDataStruct.mat'));
        
        outDir=fullfile(studyDir{2},caretDir,'suit_flat','glm4');
        
        CN=unique(T.condNum(T.condNum~=0)); % task conditions
        condNames=unique(T.condName(T.condNum~=0));
        P=size(T.data,2); % number of voxels
        
        % Tet up volume info
        Yy=zeros(length(CN),V.dim(1)*V.dim(2)*V.dim(3));
        C{1}.dim=V.dim;
        C{1}.mat=V.mat;
        
        % loop over tasks
        % calculate SSQ, variance, sd
        for c=1:length(CN),
            A=getrow(T,T.condNum==CN(c));
            A_avrg=nanmean(A.data,1);
            SSQ(c,:)=bsxfun(@minus,A_avrg,nanmean(A_avrg));
        end
        
        Yy(:,volIndx)=SSQ;
        outName='HCP_SSQ';
        
        % map vol2surf
        indices=reshape(Yy,[size(Yy,1) V.dim(1),V.dim(2),V.dim(3)]);
        for i=1:size(indices,1),
            data=reshape(indices(i,:,:,:),[C{1}.dim]);
            C{i}.dat=data;
        end
        P=caret_suit_map2surf(C,'space','SUIT','stats','nanmean','column_names',condNames');  % MK created caret_suit_map2surf to allow for output to be used as input to caret_save
        
        % save out metric
        caret_save(fullfile(outDir,sprintf('%s.metric',outName)),P);
    case 'ACTIVITY:patterns'
        outDir=fullfile(studyDir{2},caretDir,'suit_flat','glm4');
        
        load(fullfile(HCPDir,'suit','contrasts','cereb_avrgDataStruct.mat'));
        
        CN=unique(T.condNum(T.condNum~=0));
        SN=unique(T.SN(T.condNum~=0));
        
        % Set up volume info
        Yy=zeros(length(CN),V.dim(1)*V.dim(2)*V.dim(3));
        C{1}.dim=V.dim;
        C{1}.mat=V.mat;
        
        % do normalisation (subtract avrg taskCond)
        for c=1:length(CN),
            A=getrow(T,T.condNum==CN(c));
            B(c,:,:)=A.data;
        end
        % subtract baseline (avrg of all taskConds from each task)
        B_avrg=nanmean(B,1);
        Bb=bsxfun(@minus,B,B_avrg);
        
        % get group avrg
        B=permute(Bb,[1,3,2]);
        B_group=nanmean(B,3);
        
        Yy(:,volIndx)=B_group;
        
        % map vol2surf
        indices=reshape(Yy,[size(Yy,1) V.dim(1),V.dim(2),V.dim(3)]);
        for i=1:size(indices,1),
            data=reshape(indices(i,:,:,:),[C{1}.dim]);
            C{i}.dat=data;
        end
        
        P=caret_suit_map2surf(C,'space','SUIT','stats','nanmean','column_names',taskNames_new(CN));  % MK created caret_suit_map2surf to allow for output to be used as input to caret_save
        
        % save out metric
        outName='HCP_contrasts';
        caret_save(fullfile(outDir,sprintf('%s.metric',outName)),P);
    case 'ACTIVITY:vol2surf'
        taskContrast=varargin{1};
        
        % this function takes any labelled volume (already in SUIT space)
        % and plots to the surface
        HCP_subjs=sc1_sc2_HCP('PREP:HCP_subjs');
        
        someSubj=randi(100,1);
        
        mapDir=fullfile(HCPDir,'suit','contrasts',sprintf('wd%s_%d.nii',taskContrast,HCP_subjs(someSubj)));
        
        Vo=spm_vol(fullfile(mapDir));
        Vi=spm_read_vols(Vo);
        Vv{1}.dat=Vi;
        Vv{1}.dim=Vo.dim;
        Vv{1}.mat=Vo.mat;
        
        M=caret_suit_map2surf(Vv,'space','SUIT');
        
        suit_plotflatmap(M.data)
    case 'ACTIVITY:indiv'
        
        outDir=fullfile(studyDir{2},caretDir,'suit_flat','glm4');
        
        load(fullfile(HCPDir,'suit','contrasts','cereb_avrgDataStruct.mat'));
        
        CN=unique(T.condNum(T.condNum~=0));
        SN=unique(T.SN(T.SN~=0));
        P=size(T.data,2);
        B=zeros(length(SN),P,length(CN));
        
        % Set up volume info
        Yy=zeros(length(CN),V.dim(1)*V.dim(2)*V.dim(3),length(SN));
        C{1}.dim=V.dim;
        C{1}.mat=V.mat;
        
        % fix baseline here
        for c=1:length(CN),
            A=getrow(T,T.condNum==CN(c));
            A_avrg=nanmean(A.data,1);
            B(:,:,c)=bsxfun(@minus,A.data,A_avrg);
        end
        
        B=permute(B,[3 2 1]);
        Yy(:,volIndx,:)=B;
        
        % map vol2surf
        idx=1;
        indices=reshape(Yy,[size(Yy,1) V.dim(1),V.dim(2),V.dim(3) size(Yy,3)]);
        for i=1:size(indices,1), % loop over taskConds
            for ii=1:size(indices,5), % loop over subjs
                data=reshape(indices(i,:,:,:,ii),[C{1}.dim]);
                C{idx}.dat=data;
                colNames{idx}=sprintf('%s-subj%d',taskNames_new{i},ii);
                idx=idx+1;
                fprintf('%s-subj%d done \n',taskNames_new{CN(i)},ii)
            end
        end
        
        P=caret_suit_map2surf(C,'space','SUIT','stats','nanmean','column_names',colNames);  % MK created caret_suit_map2surf to allow for output to be used as input to caret_save
        
        % save out metric
        outName='HCP_contrasts_allSubjs';
        caret_save(fullfile(outDir,sprintf('%s.metric',outName)),P);
        
    case 'EVAL:HCP'       % evaluate group map on individual HCP data
        mapType=varargin{1}; % options are 'lob10','lob26','Buckner_17Networks','Buckner_7Networks', 'Cole_10Networks'
        condType=varargin{2}; % which subset of tasks are we choosing 'subset1', 'subset2' ...
        groupSize=varargin{3}; % how many subjs are we grouping together ? 10,25,50,100 etc
        
        % load in func data to test 
        load(fullfile(HCPDir,'suit','contrasts','cereb_avrgDataStruct.mat'));
        
        CN=unique(T.condNum(T.condNum~=0));
        condNames=T.condName(CN);
        SN=unique(T.SN(T.SN~=0));
        P=size(T.data,2);
        
        switch condType,
            case 'subset1' % includes all taskConds
                idx=CN';
            case 'subset2' % doesn't include cue, lf, rf, tongue, la-story
                idx=[1:5,9,11,13,14,16:23];
            case 'subset3' % doesn't include cue, lf, rf, tongue
                idx=[1:6,9,11,13,14,16:23];
        end
        
        % do normalisation (subtract avrg taskCond)
        for c=1:length(idx),
            A=getrow(T,T.condNum==idx(c));
            B(c,:,:)=A.data;
        end
        % subtract baseline (avrg of all taskConds from each task)
        B_avrg=nanmean(B,1);
        Bb=bsxfun(@minus,B,B_avrg);
        
        % group subjs into 10
        groups=[0:groupSize:100];
        for s=1:length(groups)-1,
            tmp=Bb(:,groups(s)+1:groups(s+1),:);
            Bb_group(:,s,:)=nanmean(tmp,2);
        end
        
        % load in group map
        mapName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'map.nii');
        outName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc_HCP_%s_groupSize%d.mat',condType,groupSize));
        
        volIndx=volIndx';
        
        % Now get the parcellation sampled into the same space
        [i,j,k]=ind2sub(V.dim,volIndx);
        [x,y,z]=spmj_affine_transform(i,j,k,V.mat);
        VA= spm_vol(mapName);
        [i1,j1,k1]=spmj_affine_transform(x,y,z,inv(VA.mat));
        Parcel = spm_sample_vol(VA,i1,j1,k1,0);
        % Divide the voxel pairs into all the spatial bins that we want
        fprintf('parcels\n');
        voxIn = Parcel>0;
        XYZ= [x;y;z];
        RR=[];
        [BIN,R]=mva_spatialCorrBin(XYZ(:,voxIn),'Parcel',Parcel(1,voxIn));
        clear XYZ i k l x y z i1 j1 k1 VA Parcel; % Free memory
        % Now calculate the estimation of the correlation for each subject
        for s=1:length(groups)-1,
            D=Bb_group(:,s,voxIn);
            D=permute(D,[1,3,2]); % conditions x voxels
            fprintf('%d cross\n',s);
            R.SN = ones(length(R.N),1)*s;
            R.corr=mva_spatialCorr(D,BIN);
            R.crossval = zeros(length(R.corr),1);
            RR = addstruct(RR,R);
        end;
        save(outName,'-struct','RR');
    case 'EVAL:PLOT:CURVES'
        mapType=varargin{1}; % options are 'lob10','lob26','bucknerRest','SC<studyNum>_<num>cluster', or 'SC<studyNum>_POV<num>'
        condType=varargin{2}; % 'subset1', 'subset2' etc
        
        vararginoptions({varargin{3:end}},{'CAT','sn'}); % option if doing individual map analysis
        
        T=load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc_HCP_%s.mat',condType)));
        
        % distances are diff across evals so need to get dist per bin:
        for b=1:length(unique(T.bin)),
            dist=mode(round(T.dist(T.bin==b)));
            idx=find(T.bin==b);
            T.dist(idx,1)=dist;
        end
        
        if exist('CAT'),
            xyplot(T.dist,T.corr,T.dist,'split',T.bwParcel,'subset',T.crossval==0 & T.dist<=35,'CAT',CAT,'leg',{'within','between'},'leglocation','SouthEast');
        else
            xyplot(T.dist,T.corr,T.dist,'split',T.bwParcel,'subset',T.crossval==0 & T.dist<=35,'leg',{'within','between'},'leglocation','SouthEast');
        end
    case 'EVAL:STATS:CURVES'
        mapType=varargin{1}; % options are 'lob10','lob26','bucknerRest','SC<studyNum>_<num>cluster', or 'SC<studyNum>_POV<num>'
        condType=varargin{2};
        crossval=0;
        
        T=load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc_HCP_%s.mat',condType)));
        
        % do stats (over all bins) for group only
        C=getrow(T,T.crossval==crossval & T.dist<=35); % only crossval and dist<35
        S=tapply(C,{'bwParcel','SN'},{'corr'});
        fprintf('overall \n')
        ttest(S.corr(S.bwParcel==0), S.corr(S.bwParcel==1),2,'paired');
        
        % calculate effect size
        Group1=S.corr(S.bwParcel==0);
        Group2=S.corr(S.bwParcel==1);
        
        num=((Group1-1)*std(Group1)^2 + (Group2-1)*std(Group2)^2);
        denom=Group1+Group2-2;
        
        pooledSTD= sqrt(mean(num)/mean(denom));
        
        ES_pooled=(mean(Group1)-mean(Group2))/pooledSTD;
        
        fprintf('Effect size for within and between for %s is %2.2f when denom is pooled std  \n',mapType,ES_pooled);
        
        % summary stats
        x1=nanmean(S.corr(S.bwParcel==0));x2=nanmean(S.corr(S.bwParcel==1));
        SEM1=std(S.corr(S.bwParcel==0))/sqrt(length(T.SN));SEM2=std(S.bwParcel==1)/sqrt(length(T.SN));
        fprintf('average within corr is %2.2f; CI:%2.2f-%2.2f \n average between corr is %2.2f; CI:%2.2f-%2.2f \n',...
            nanmean(S.corr(S.bwParcel==0)),x1-(1.96*SEM1),x1+(1.96*SEM1),nanmean(S.corr(S.bwParcel==1)),...
            x2-(1.96*SEM2),x2+(1.96*SEM2));
    case 'EVAL:PLOT:DIFF'
        mapType=varargin{1}; % options are 'lob10','lob26','bucknerRest','SC<studyNum>_<num>cluster', or 'SC<studyNum>_POV<num>'
        condType=varargin{2}; % 'subset1', 'subset2' etc
        
        vararginoptions({varargin{3:end}},{'CAT','sn'}); % option if plotting individual map analysis
        
        T=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc_HCP_%s.mat',condType)));  %sprintf('spatialBoundfunc_HCP_%s.mat',condType)
        
        % distances are diff across evals so need to get dist per bin:
        for b=1:length(unique(T.bin)),
            dist=mode(round(T.dist(T.bin==b)));
            idx=find(T.bin==b);
            T.dist(idx,1)=dist;
        end
        
        % plot boxplot of different clusters
        W=getrow(T,T.bwParcel==0); % within
        B=getrow(T,T.bwParcel==1); % between
        W.diff=W.corr-B.corr;
        W=rmfield(W,{'bwParcel','crossval','corr'});
        
        if exist('CAT'),
            lineplot(W.dist,W.diff,'subset',W.dist<=35,'CAT',CAT,'leg','auto');
        else
            lineplot(W.dist,W.diff,'subset',W.dist<=35);
        end
    case 'EVAL:STATS:DIFF'
        mapType=varargin{1}; % options are 'lob10','lob26','bucknerRest','SC<studyNum>_<num>cluster', or 'SC<studyNum>_POV<num>'
        condType=varargin{2}; % 'subset3_groupSize25'
        crossval=0;
        
        % do stats
        P=[];
        for m=1:length(mapType),
            T=load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType{m}),sprintf('spatialBoundfunc_HCP_%s.mat',condType)));
            
            A=getrow(T,T.crossval==crossval);
            A.type=repmat({sprintf('%d.%s',m,mapType{m})},length(A.bin),1);
            A.m=repmat(m,length(A.bin),1);
            P=addstruct(P,A);
            clear A
        end
        
        W=getrow(P,P.bwParcel==0); % within
        B=getrow(P,P.bwParcel==1); % between
        W.diff=W.corr-B.corr;
        % do stats (integrate over spatial bins)
        W=rmfield(W,{'bwParcel','crossval','corr'});
        C=getrow(W,W.dist<=35);
        S=tapply(C,{'m','SN','type'},{'diff'});
        
        % do F test (or t test if just two groups)
        if length(unique(S.m))>2,
            X=[S.diff(S.m==1),S.diff(S.m==2),S.diff(S.m==3)];
            results=anovaMixed(S.diff,S.SN,'between',S.m,'betweenNames',{'a','b','c'});
            
        else
            ttest(S.diff(S.m==1), S.diff(S.m==2),2,'paired');
            
            % calculate effect size
            Group1=S.diff(S.m==1);
            Group2=S.diff(S.m==2);
            
            ES_group1=(mean(Group1)-mean(Group2))/std(Group1); % uses the std of one of the groups
            ES_group2=(mean(Group1)-mean(Group2))/std(Group2); % uses the std of one of the groups
            % this is biased as the effect size changes depending on which
            % group you choose. Therefore, pooled estimate is better.
            
            num=((Group1-1)*std(Group1)^2 + (Group2-1)*std(Group2)^2);
            denom=Group1+Group2-2;
            
            pooledSTD= sqrt(mean(num)/mean(denom));
            
            ES_pooled=(mean(Group1)-mean(Group2))/pooledSTD;
            
            fprintf('Effect size between %s and %s is %2.2f when denom is std(Group1) \n',mapType{1},mapType{2},ES_group1);
            fprintf('Effect size between %s and %s is %2.2f when denom is std(Group2) \n',mapType{1},mapType{2},ES_group2);
            fprintf('Effect size between %s and %s is %2.2f when denom is pooled std  \n',mapType{1},mapType{2},ES_pooled);
            
        end

    case 'STRENGTH:eval_bound'
        mapType = varargin{1}; % 'SC12_10cluster','Buckner_7Networks'
        condType = varargin{2}; % 'subset1' etc
        groupSize=varargin{3}; % 1,10,25 etc
        
        spatialBins = [0:3:35];

        EvalDir = fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType));
        
        % load boundaries
        load(fullfile(EvalDir,'boundaries.mat'));
        numBins = length(spatialBins)-1;
        
        % Get the condition numbers
        P=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        CN=unique(P.condNum(P.condNum~=0));
        
        % load in func data to test
        load(fullfile(HCPDir,'suit','contrasts','cereb_avrgDataStruct.mat'));
        
        switch condType,
            case 'subset1' % includes all taskConds
                idx=CN';
            case 'subset2' % doesn't include cue, lf, rf, tongue, la-story
                idx=[1:5,9,11,13,14,16:23];
            case 'subset3' % doesn't include cue, lf, rf, tongue
                idx=[1:6,9,11,13,14,16:23];
        end
        
        % do normalisation (subtract avrg taskCond)
        for c=1:length(idx),
            A=getrow(T,T.condNum==idx(c));
            B(c,:,:)=A.data;
        end
        % subtract baseline (avrg of all taskConds from each task)
        B_avrg=nanmean(B,1);
        Bb=bsxfun(@minus,B,B_avrg);
        
        % group subjs into 10
        groups=[0:groupSize:100];
        for s=1:length(groups)-1,
            tmp=Bb(:,groups(s)+1:groups(s+1),:);
            Bb_group(:,s,:)=nanmean(tmp,2);
        end
        
        RR=[];
        % Now build a structure all boundaries
        for i=1:size(Edge,1)
            indx = Cluster==Edge(i,1) | Cluster==Edge(i,2);  % Get all the involved voxels
            fprintf('Edge %d %d\n',Edge(i,1),Edge(i,2));
            
            % Determine the spatial bins for this pair of regions
            [BIN,R]=mva_spatialCorrBin(coords(:,indx),'Parcel',Cluster(indx)','spatialBins',spatialBins);
            N=length(R.N);
            
            % Now determine the correlation for each subject
            for s=1:length(groups)-1;
                D=Bb_group(:,s,indx);
                D=permute(D,[1,3,2]); % conditions x voxels
                fprintf('%d cross\n',s);
                R.Edge = repmat(Edge(i,:),N,1);
                R.SN = ones(N,1)*s;
                R.corr=mva_spatialCorr(D,BIN,'numBins',N);
                R.crossval = zeros(N,1);
                if (length(R.corr)~=N)
                    keyboard;
                end;
                RR=addstruct(RR,R);
            end;
        end;
        save(fullfile(EvalDir,sprintf('BoundariesFunc3_HCP_%s_groupSize%d.mat',condType,groupSize)),'-struct','RR');
    case 'STRENGTH:visualise_bound'
        mapType = varargin{1}; % 'SC12_10cluster','Buckner_7Networks'
        condType = varargin{2}; % 'subset1' etc
        groupSize=varargin{3}; % 1,10,25 etc

        bcolor ='k';
        opacacy = 0.5;
        bscale = 600;
        bmax=20;
        
        vararginoptions(varargin(4:end),{'bcolor','opacacy','bscale','bmax'});
        
        EvalDir = fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType));
        SurfDir = fullfile(studyDir{1},'surfaceCaret','suit_flat');
        
        % Get boundaries for the map
        load(fullfile(EvalDir,'boundaries.mat'));
        
        % Evaluate the strength of each border
        T=load(fullfile(EvalDir,sprintf('BoundariesFunc3_HCP_%s_groupSize%d.mat',condType,groupSize)));
        
        % Map the clusters
        V.dat=zeros([V.dim(1) V.dim(2) V.dim(3)]);
        V.dat(volIndx)=Cluster;
        C{1}=V;
        C{1}=rmfield(C{1},{'fname','dt','pinfo','n','descrip','private'});
        Mcl=caret_suit_map2surf(C,'space','SUIT','stats','mode','stats',@mode);
        Mcl=Mcl.data;
        
        % Map the parcel
        C{1}.dat(volIndx)=Parcel;
        Mpa=caret_suit_map2surf(C,'space','SUIT','stats','mode','stats',@mode);
        Mpa=Mpa.data;
        
        % Determine the border points
        COORD=gifti(fullfile('FLAT.coord.gii'));
        TOPO=gifti(fullfile('CUT.topo.gii'));
        
        % Make matrix of all the unique edges of the flatmap
        Tedges=[TOPO.faces(:,[1 2]);TOPO.faces(:,[2 3]);TOPO.faces(:,[1 3])]; % Take 3 edges from the faces
        Tedges= [min(Tedges,[],2) max(Tedges,[],2)];     % Sort in ascending order
        Tedges = unique(Tedges,'rows');                  % Only retain each edge ones
        EdgeCl=Mcl(Tedges);                              % Which cluster does each node belong to?
        EdgeCl= [min(EdgeCl,[],2) max(EdgeCl,[],2)];     % Sort in ascending order
        
        % Assemble the edges that lie on the boundary between clusters
        for  i=1:size(Edge,1)
            indxEdge = find(EdgeCl(:,1)==Edge(i,1) & EdgeCl(:,2)==Edge(i,2));
            Border(i).numpoints=length(indxEdge);
            for e=1:length(indxEdge)
                % find the boundary point: In the middle of the edge
                Border(i).data(e,:)=(COORD.vertices(Tedges(indxEdge(e),1),:)+...
                    COORD.vertices(Tedges(indxEdge(e),2),:))/2; % Average of coordinates
            end;
        end;

        for  i=1:size(Edge,1)
            % Make sure that the bin is calculated both for within and
            % between
            A=pivottable(T.bin,T.bwParcel,T.corr,'nanmean','subset',T.crossval==0 & ...
                T.Edge(:,1)==Edge(i,1) & T.Edge(:,2)==Edge(i,2));
            EdgeWeight(i,1)=nanmean(A(:,1)-A(:,2));  % Difference within - between
        end;
        
        % check if a color map is provided
        if exist(fullfile(EvalDir,'colourMap.txt'),'file'),
            cmap=load(fullfile(EvalDir,'colourMap.txt'));
            cmap=cmap(:,2:4);
            cmapF=cmap/255;
        else
            cmapF=colorcube(max(Parcel));
        end
        %         cmapF = bsxfun(@plus,cmapF*opacacy,[1 1 1]*(1-opacacy));
        
        % Make the plot
        suit_plotflatmap(Mpa,'type','label','border',[],'cmap',cmapF);
        whitebg
        hold on;
        LineWeight=EdgeWeight*bscale;
        LineWeight(LineWeight>bmax)=bmax;
        for  b=1:length(Border)
            if (Border(b).numpoints>0 & EdgeWeight(b)>0),
                p=plot(Border(b).data(:,1),Border(b).data(:,2),'k.');
                set(p,'MarkerSize',LineWeight(b));
                weights(b)=EdgeWeight(b);
                % plot DCBC of each functional boundary ?
%                                 p=text(double(Border(b).data(1,1)),double(Border(b).data(1,2)),sprintf('%2.3f',weights(b)));
%                                 set(p,'FontSize',20);
            end;
        end
        hold off;
        tmp=weights(weights~=0);
        
        fprintf('min diff is %2.5f and max diff is %2.2f \n', min(tmp),max(tmp));
        
    case 'AXES:group_curves' % make separate graphs for 'lob10','Buckner_7Networks','Buckner_17Networks','Cole_10Networks','SC12_10cluster'
        toPlot=varargin{1}; % 'SC12_10cluster'
        condType=varargin{2}; % 'subset1', 'subset2' etc
        
        % Aesthetics
        CAT.markertype='none';
        CAT.errorwidth=.5;
        CAT.linecolor={'r','k'};
        CAT.errorcolor={'r','k'};
        CAT.linewidth={2, 2};
        CAT.linestyle={'-','-'};
        
        sc1_sc2_HCP('EVAL:PLOT:CURVES',toPlot,condType,'CAT',CAT);
        
        % Labelling
        set(gca,'YLim',[0 0.55],'XLim',[0 35],'FontSize',14,'xtick',[0:5:35],'XTickLabel',{'0','','','','','','','35'}); %
        xlabel('Spatial Distances (mm)');
        ylabel('Activity Correlation (R)');
        %         title(plotName);
        set(gcf,'units','centimeters','position',[5,5,15,15])
        %         axis('auto')
        % do stats
        %         sc1_sc2_ICB('EVAL:STATS:CURVES',toPlot)
    case 'AXES:diff_curves' % make summary graph for diff curves for all maps
        toPlot=varargin{1}; % {'SC12_10cluster','Buckner_17Networks'}
        condType=varargin{2}; % 'subset1', 'subset2' etc
        
        % aesthetics
        CAT.errorwidth=.5;
        CAT.markertype='none';
        CAT.linewidth=3;
        CAT.linestyle={'-','-','-','-','-','-'};
        CAT.linewidth={2, 2, 2, 2, 2, 2};
        errorcolor={'w','r','b','g','r','r'};
        linecolor={'w','r','b','g','r','r'};
        
        %         errorcolor={[0 0 0],[0 50/255 150/255],[44/255 26/255 226/255],[0 150/255 255/255],[185/255 0 54/255],[139/255 0 123/255],[0 158/255 96/255],[0 158/255 96/255]};
        %         linecolor={[0 0 0],[0 50/255 150/255],[44/255 26/255 226/255],[0 150/255 255/255],[185/255 0 54/255],[139/255 0 123/255],[0 158/255 96/255],[0 158/255 96/255]};
        %
        for m=1:length(toPlot),
            CAT.errorcolor=errorcolor{m};
            CAT.linecolor=linecolor{m};
            sc1_sc2_HCP('EVAL:PLOT:DIFF',toPlot{m},condType,'CAT',CAT); % always take crossval + unique
            hold on
        end
        hold off
        
        % Labelling
        set(gca,'YLim',[-0.02 0.08],'FontSize',12,'XLim',[0 35],'xtick',[0:5:35],'XTickLabel',{'0','','','','','','','35'});
        xlabel('Spatial Distances (mm)');
        ylabel('Difference');
        set(gcf,'units','centimeters','position',[5,5,9,12])
        %         legend(plotName,'Location','NorthWest')
        
        % do stats
        %         sc1_sc2_functionalAtlas('EVAL:STATS:DIFF',toPlot,evalNums,'group',1,'unique'); % always take crossval + unique
end


% Local functions
function dircheck(dir)
if ~exist(dir,'dir');
    warning('%s doesn''t exist. Creating one now. You''re welcome! \n',dir);
    mkdir(dir);
end