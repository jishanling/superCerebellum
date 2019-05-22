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
fshem     = {'lh','rh'};
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
            surf_resliceFS2WB(subj_name{i},fsDir,atlasDir,wbDir,'resolution','32k');
        end;
    case 'SURF:map_con'         % STEP 11.5: Map con / ResMS (.nii) onto surface (.gifti)
        sn    = returnSubjs;     % subjNum
        study = [1 2];
        glm  = 4;     % glmNum
        hemis = [1 2];
        resolution = '32k';
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        vararginoptions(varargin,{'sn','study','glm','what','hemis','resolution'});
        
        for s=sn
            for h=hemis
                surfDir = fullfile(wbDir, subj_name{s});
                white=fullfile(surfDir,sprintf('%s.%s.white.%s.surf.gii',subj_name{s},Hem{h},resolution));
                pial=fullfile(surfDir,sprintf('%s.%s.pial.%s.surf.gii',subj_name{s},Hem{h},resolution));
                C1=gifti(white);
                C2=gifti(pial);
                
                for st = study
                    glmDir =fullfile(baseDir,studyDir{st},sprintf('GLM_firstlevel_%d',glm),subj_name{s});
                    T=getrow(D,D.StudyNum==st);
                    filenames={};
                    for i=1:length(T.condNames);
                        filenames{i} = fullfile(glmDir,sprintf('con_%s-rest.nii',T.condNames{i}));
                    end;
                    filenames{i+1} = fullfile(glmDir,'ResMS.nii');
                    T.condNames{i+1}= 'ResMS';
                    outfile = fullfile(surfDir,sprintf('%s.%s.%s.con.%s.func.gii',subj_name{s},Hem{h},studyDir{st},resolution));
                    
                    G=surf_vol2surf(C1.vertices,C2.vertices,filenames,'column_names',T.condNames,'anatomicalStruct',hemname{h});
                    save(G,outfile);
                    
                    fprintf('mapped %s %s %s \n',subj_name{s},Hem{h},studyDir{st});
                end;
            end;
        end
    case 'SURF:combine_con'    %  Does univariate noise normalization (combined) and then combines the contrast files by the common baseline
        sn    = returnSubjs;     % subjNum
        glm  = 4;     % glmNum
        hemis = [1 2];
        resolution = '32k';
        T=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        vararginoptions(varargin,{'sn','glm','what','hemis','resolution'});
        for s=sn
            surfDir = fullfile(wbDir, subj_name{s});
            for h=hemis
                for st=[1 2]
                    G{st}=gifti(fullfile(surfDir,sprintf('%s.%s.%s.con.%s.func.gii',subj_name{s},Hem{h},studyDir{st},resolution)));
                    N=size(G{st}.cdata,1);
                    wcon{st} = [G{st}.cdata(:,2:end-1) zeros(N,1)]; % Add rest
                    Ts = getrow(T,T.StudyNum==st);
                    wcon{st}=bsxfun(@minus,wcon{st},mean(wcon{st}(:,Ts.overlap==1),2));
                end;
                resMS=(G{1}.cdata(:,end)+G{1}.cdata(:,end))/2;
                W = [wcon{1} wcon{2}];
                W = bsxfun(@rdivide,W,resMS); % Common noise normalization
                outfile = fullfile(surfDir,sprintf('%s.%s.wcon.%s.func.gii',subj_name{s},Hem{h},resolution));
                Go=surf_makeFuncGifti(W,'columnNames',T.condNames,'anatomicalStruct',hemname{h});
                save(Go,outfile);
                fprintf('combined %s %s \n',subj_name{s},Hem{h});
            end;
        end;
    case 'SURF:groupFiles'
        hemis=[1 2];
        sn = returnSubjs;
        cd(fullfile(wbDir,'group164k'));
        for h=hemis
            inputFiles = {};
            for s=1:length(sn)
                inputFiles{s}  = fullfile(wbDir, subj_name{sn(s)},sprintf('%s.%s.wcon.func.gii',subj_name{sn(s)},Hem{h}));
                columnName{s} = subj_name{sn(s)};
            end;
            groupfile=sprintf('group.wcon.%s.func.gii',Hem{h});
            outfilenamePattern=sprintf('wcon.%%s.%s.func.gii',Hem{h});
            surf_groupGiftis(inputFiles,'groupsummary',groupfile,'outcolnames',subj_name(sn),'outfilenamePattern',outfilenamePattern);
        end;
    case 'SURF:resample32k'  % Resample functional data from group164 to group32
        hemis=[1 2];
        sn = [2,3,4];
        sourceDir =fullfile(wbDir,'group164k');
        targetDir =fullfile(wbDir,'group32k');
        T=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        for h=hemis
            % wb_command -metric-resample group164K/group.wcon.L.func.gii group164K/fs_LR.164k.L.sphere.surf.gii group32K/fs_LR.32k.L.sphere.surf.gii BARYCENTRIC group32K/group.wcon.L.func.gii
            % wb_command -metric-resample group164K/group.wcon.R.func.gii group164K/fs_LR.164k.R.sphere.surf.gii group32K/fs_LR.32k.R.sphere.surf.gii BARYCENTRIC group32K/group.wcon.R.func.gii
        end;
    case 'PARCEL:annot2labelgii'
        filename = varargin{1};
        for i=1:2
            name = [fshem{i} '.' filename '.annot'];
            [v,label,colorT]=read_annotation(name);
            values=colorT.table(:,5);
            newlabel = zeros(size(label));
            for j=1:length(values)
                newlabel(label==values(j))=j-1;  % New label values starting from 0
            end;
            RGB = colorT.table(:,1:3)/256; % Scaled between 0 and 1
            G=surf_makeLabelGifti(newlabel,'anatomicalStruct',hemname{i},'labelNames',colorT.struct_names,'labelRGBA',[RGB ones(size(RGB,1),1)]);
            save(G,[fshem{i} '.' filename '.label.gii']);
        end;
    case 'PARCEL:fsaverage2FSLR'
        infilename = varargin{1}; % '~/Data/Atlas_templates/CorticalParcellations/Yeo2011/Yeo_JNeurophysiol11_FreeSurfer/fsaverage/label/lh.Yeo2011_17Networks_N1000.label.gii'
        outfilename= varargin{2}; % '~/Data/Atlas_templates/FS_LR_164/Yeo_JNeurophysiol11_17Networks.164k.L.label.gii';
        hem        = varargin{3};
        surf       = varargin{4}; % '164k','32k'
        inSphere   = fullfile(atlasDir,['fs_' Hem{hem}],['fsaverage.' Hem{hem} '.sphere.164k_fs_' Hem{hem} '.surf.gii']);
        outSphere  = fullfile(atlasDir,'resample_fsaverage',['fs_LR-deformed_to-fsaverage.' Hem{hem} '.sphere.' surf '_fs_LR.surf.gii']);
        
        % Convert surface to Gifti
        system(['wb_command -metric-resample ' infilename ' ' inSphere ' ' outSphere ' ADAP_BARY_AREA ' outfilename ' -area-surfs ' inSphere ' ' outSphere]);
    case 'DCBC:computeDistances'
        sn=returnSubjs;
        hem = [1 2];
        resolution = '32k';
        maxradius = 40;
        
        vararginoptions(varargin,{'sn','hem','resolution','maxradius'});
        for s=sn
            for h=hem
                surfDir = fullfile(wbDir, subj_name{s});
                white=fullfile(surfDir,sprintf('%s.%s.white.%s.surf.gii',subj_name{s},Hem{h},resolution));
                pial=fullfile(surfDir,sprintf('%s.%s.pial.%s.surf.gii',subj_name{s},Hem{h},resolution));
                C1=gifti(white);
                C2=gifti(pial);
                vertices = double((C1.vertices' + C2.vertices' )/2); % Midgray surface
                faces    = double(C1.faces');
                numVert=size(vertices,2);
                n2f=surfing_nodeidxs2faceidxs(faces);
                D=inf([numVert numVert],'single');
                for i=1:numVert;
                    [subV, subF, subIndx,vidxs, fidxs]=surfing_subsurface(vertices, faces, i, maxradius, n2f); % construct correct sub surface
                    D(vidxs,i)=single(surfing_dijkstradist(subV,subF,subIndx,maxradius));
                    if (mod(i,100)==0)
                        fprintf('.');
                    end
                end;
                fprintf('\n');
                save(fullfile(surfDir,sprintf('distances.%s.mat',Hem{h})),'D','-v7.3');
            end;
        end;
    case 'DCBC:avrgdistances'
        sn=[2,3,4,6,8,9,10,12,14];
        hem = [1 2];
        resolution = '32k';
        vararginoptions(varargin,{'sn','hem','resolution','maxradius'});
        avrgD=single(zeros(32492,32492));
        N=length(sn);
        for h=hem
            for s=1:length(sn)
                fprintf('.');
                surfDir = fullfile(wbDir, subj_name{sn(s)});
                load(fullfile(surfDir,sprintf('distances.%s.mat',Hem{h})));
                D(isinf(D))=45;
                w=single((s-1)/N);
                avrgD = w.*avrgD+(1-w)*D;
            end;
            fprintf('\n');
            avrgDs=sparse(double(avrgD));
            save(fullfile(wbDir,'group32k',sprintf('distances_sp.%s.mat',Hem{h})),'avgrDs');
        end;
    case 'Eval:DCBC'
        sn=returnSubjs;
        hem = [1 2];   
        resolution = '32k';
        taskSet = [1 2];   
        condType = 'unique'; % Evaluate on all or only unique conditions? 
        bins = [0:5:40];     % Spatial bins in mm 
        parcel = [];         % N*2 matrix for both hemispheres 
        RR=[];
        vararginoptions(varargin,{'sn','hem','bins','parcel','condType','taskSet','resolution'});
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        numBins = numel(bins)-1;
        for h=hem
            load(fullfile(wbDir,'group32k',sprintf('distances_sp.%s.mat',Hem{h})));
            
            % Now find the pairs that we care about to safe memory 
            vertIdx = find(parcel(:,h)>0);
            avrgDs = avrgDs(vertIdx,vertIdx);
            par    = parcel(vertIdx,h); 
            
            [row,col,avrgD]=find(avrgDs);
            inSp = sub2ind(size(avrgDs),row,col);
            sameReg=(bsxfun(@ne,par',par)+1);
            sameReg=sameReg(inSp);
            clear avrgDs par; 
            for ts = taskSet;
                D1=getrow(D,D.StudyNum==ts);
                switch condType,
                    case 'unique'
                        % if funcMap - only evaluate unique tasks in sc1 or sc2
                        condIdx=D1.condNum(D1.overlap==0); % get index for unique tasks
                    case 'all'
                        condIdx=D1.condNum;
                end
                for s=sn
                    A=gifti(fullfile(wbDir,subj_name{s},sprintf('%s.%s.%s.con.%s.func.gii',subj_name{s},Hem{h},studyDir{ts},resolution)));
                    Data = [A.cdata(:,2:end-1) zeros(size(A.cdata,1),1)];
                    Data = bsxfun(@rdivide,Data,sqrt(A.cdata(:,end)));
                    Data = Data(vertIdx,condIdx); % Take the right subset 
                    Data = bsxfun(@minus,Data,mean(Data,2));
                    Data = single(Data');
                    [K,P]=size(Data); 
                    clear A; 
                    
                    SD = sqrt(sum(Data.^2)/K);
                    
                    VAR = (SD'*SD);
                    COV = Data'*Data/K;
                    fprintf('%d',s);
                    for bw=[1 2]
                        for i=1:numBins
                            fprintf('.');
                            in = i+(bw-1)*numBins;
                            inBin = avrgD>bins(i) & avrgD<=bins(i+1) & sameReg==bw;
                            R.SN(in,1)      = s;
                            R.hem(in,1)     = h; 
                            R.studyNum(in,1) = ts; 
                            R.N(in,1)       = sum(inBin(:));
                            R.avrDist(in,1) = mean(avrgD(inBin));
                            R.bwParcel(in,1)= bw-1;
                            R.bin(in,1)     = i;
                            R.distmin(in,1) = bins(i);
                            R.distmax(in,1) = bins(i+1);
                            R.meanVAR(in,1) = full(nanmean(VAR(inSp(inBin))));
                            R.meanCOV(in,1) = full(nanmean(COV(inSp(inBin))));
                        end;
                    end;
                    clear VAR COV;
                    
                    
                    R.corr=R.meanCOV./sqrt(R.meanVAR);
                    fprintf('\n');
                    RR = addstruct(RR,R);
                end;
            end;
        end;
        varargout={RR}; 
    case 'EVAL:doEval' 
%         for h=1:2
%             A=gifti(sprintf('Yeo_JNeurophysiol11_7Networks.32k.%s.label.gii',Hem{h})); 
%             parcel(:,h)=A.cdata; 
%         end; 
%         T=sc1_sc2_neocortical('Eval:DCBC','hem',[1 2],'parcel',parcel,'condType','all');
%         save('Eval_Yeo7_all.mat','-struct','T'); 

%         for h=1:2
%             A=gifti(sprintf('Yeo_JNeurophysiol11_17Networks.32k.%s.label.gii',Hem{h})); 
%             parcel(:,h)=A.cdata; 
%         end; 
%         T=sc1_sc2_neocortical('Eval:DCBC','hem',[1 2],'parcel',parcel,'condType','all');
%         save('Eval_Yeo17_all.mat','-struct','T'); 

%         for h=1:2
%             A=gifti(sprintf('Yeo_JNeurophysiol11_17Networks.32k.%s.label.gii',Hem{h})); 
%             parcel(:,h)=A.cdata; 
%         end; 
%         T=sc1_sc2_neocortical('Eval:DCBC','hem',[1 2],'parcel',parcel,'condType','all');
%         save('Eval_Yeo17_all.mat','-struct','T'); 

%         for h=1:2
%             A=gifti(sprintf('Yeo_JNeurophysiol11_17Networks.32k.%s.label.gii',Hem{h})); 
%             parcel(:,h)=A.cdata; 
%         end; 
%         T=sc1_sc2_neocortical('Eval:DCBC','hem',[1 2],'parcel',parcel,'condType','all');
%         save('Eval_Yeo17_all.mat','-struct','T'); 
        
        for h=1:2 
            D=ft_read_cifti(sprintf('Glasser_2016.32k.%s.dlabel.nii',Hem{h})); 
            parcel(:,h)=D.x1; 
        end; 
        parcel(isnan(parcel))=0; 
        T=sc1_sc2_neocortical('Eval:DCBC','hem',[1 2],'parcel',parcel,'condType','all');
        save('Eval_Glasser_all.mat','-struct','T'); 

    case 'EVAL:plotEval' 
        toPlot={'Yeo17','Yeo7'}; 
        CAT.linecolor={'r','b','k','g'}; 
        CAT.linestyle={'-','-','-','-'}; 
        CAT.linewidth=2; 
        CAT.markertype='none'; 
        CAT.errorcolor={'r','b','k','g'}; 
        condType='all'; 
        vararginoptions(varargin,{'toPlot','condType'}); 
        T=[]; 
        for i=1:length(toPlot)
            D=load(sprintf('Eval_%s_%s.mat',toPlot{i},condType)); 
            TT=tapply(D,{'bin','SN','distmin','distmax'},{'corr','mean','name','corrB','subset',D.bwParcel==0},...
                                                         {'corr','mean','name','corrW','subset',D.bwParcel==1}); 
            TT.parcel=repmat(toPlot(i),length(TT.SN),1); 
            T=addstruct(T,TT); 
        end;
        T.DCBC=T.corrB-T.corrW; 
        T.binC = (T.distmin+T.distmax)/2; 
        lineplot(T.binC,T.DCBC,'CAT',CAT,'split',T.parcel,'leg',toPlot); 
end;