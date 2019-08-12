function varargout=sc1sc2_connectivity(what,varargin)
% Matlab script to do connectivity analyses for sc1 and sc2 experiments
% Does time series extraction
%==========================================================================
type=[];
% % (1) Directories
% rootDir           = '/Users/jdiedrichsen/Data/super_cerebellum';
rootDir         = '/Volumes/MotorControl/data/super_cerebellum_new';
wbDir           = fullfile(rootDir,'sc1','surfaceWB');
sc1Dir          = [rootDir '/sc1'];
sc2Dir          = [rootDir '/sc2'];
behavDir        = ['data'];
imagingDir      = ['imaging_data'];
imagingDirRaw   = ['imaging_data_raw'];
dicomDir        = ['imaging_data_dicom'];
anatomicalDir   = ['anatomicals'];
suitDir         = ['suit'];
caretDir        = ['surfaceCaret'];
regDir          = ['RegionOfInterest'];
connDir         = ['connectivity_cerebellum'];

%==========================================================================

% % (2) Hemisphere and Region Names
Hem       = {'L','R'};
hemname   = {'CortexLeft','CortexRight'};
subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11',...
    's12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24',...
    's25','s26','s27','s28','s29','s30','s31'};
returnSubjs=[2,3,4,6,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31];
expStr    = {'sc1','sc2'};

%==========================================================================
switch(what)
    case 'make_bucknermargin'   % Generates the Buckner Margin ROI
        kernel = 8;
        threshold = 0.05;
        for s=2:22
            Vcort=spm_vol(fullfile(sc1Dir,regDir,'data',subj_name{s},'cortical_mask_grey_corr.nii'));
            Vcereb=spm_vol(fullfile(sc1Dir,suitDir,'anatomicals',subj_name{s},'cereb_prob_corr_grey.nii'));
            % smooth cerebellar image using 6 mm kernel
            Vcort.dat = spm_read_vols(Vcort);
            X = spm_read_vols(Vcereb);
            Vcereb.dat = zeros(size(X));
            spm_smooth(X,Vcereb.dat,kernel);
            
            linind = find(Vcort.dat>0);
            [i,j,k]=ind2sub(Vcort.dim,linind);
            [i1,j1,k1]=spmj_affine_transform(i,j,k,inv(Vcereb.mat)*Vcort.mat);
            Vcereb.dt =[64 0]; Vcereb.pinfo=[1 0 0]';
            Data = spm_sample_vol(Vcereb,i1,j1,k1,1);
            Vcort.dat(linind(Data<threshold))=0; % Set these parts to zero
            Vcort.fname =  (fullfile(sc1Dir,regDir,'data',subj_name{s},'bucknermargin.nii'));
            spm_write_vol(Vcort,Vcort.dat);
        end;
    case 'TS_correlation_map'   % Seed based correlation analysis for 1 location - for visualization
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        location=[-8 -50 -11]; % [x,y,z]
        
        glmDir =[sc1Dir sprintf('/GLM_firstlevel_%d',glm)];
        subjs=length(sn);
        
        for s=1:subjs,
            glmDirSubj=fullfile(glmDir, subj_name{sn(s)});
            load(fullfile(glmDirSubj,'SPM_light.mat'));
            T=load(fullfile(glmDirSubj,'SPM_info.mat'));
            
            % load data
            tic;
            SPM=spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{sn(s)}));
            % Get the raw data files
            V=SPM.xY.VY;
            fprintf('Moved %d\n',toc);
            
            % Determine locations in image
            VresMS = spm_vol(fullfile(glmDirSubj,'ResMS.nii'));
            XresMS =spm_read_vols(VresMS);
            linIndx = find(XresMS>0);
            [i,j,k]=ind2sub(VresMS.dim,linIndx);
            SS1 = zeros(2,length(i));
            SS2 = zeros(2,1);
            SSC = zeros(2,length(i));
            
            [i1,j1,k1]=spmj_affine_transform(location(1),location(2),location(3),inv(VresMS.mat));
            tar=[round(i1) round(j1) round(k1)];
            target = findrow([i j k],tar);
            if (isempty(target))
                error('target not in volume');
            end;
            % Now loop over runs and get the
            tic;
            for r=1:16
                fprintf('%d.',r);
                row  = SPM.Sess(r).row;
                col  = SPM.Sess(r).col;  % Regressors of interest
                col1 = [SPM.Sess(r).col SPM.xX.iB(r)];
                Y = zeros(length(row),length(i));
                for t=1:length(row)
                    Y(t,:) = spm_sample_vol(V(row(t)),i,j,k,0)';  % Data is N x P
                end;
                Yfilt = SPM.xX.W(row,row)*Y;
                B = SPM.xX.pKX(col1,row)*Yfilt;                             %-Parameter estimates
                Yhat = SPM.xX.xKXs.X(row,col1)*B; %- predicted values
                Yres = Yfilt - Yhat; %- residuals
                
                Ya{1}=bsxfun(@minus,Yhat,mean(Yhat));
                Ya{2}=bsxfun(@minus,Yres,mean(Yres));
                
                for ty=1:2
                    SS1(ty,:)=SS1(ty,:)+sum(Ya{ty}.^2);
                    SS2(ty,:)=SS2(ty,:)+sum(Ya{ty}(:,target).^2);
                    SSC(ty,:)=SSC(ty,:)+Ya{ty}(:,target)'*Ya{ty};
                end;
            end;
            DE = sqrt(bsxfun(@times,SS1,SS2));
            CR = SSC./DE;
            for ty=1:2
                Vo= VresMS;
                Vo.fname = fullfile(glmDirSubj,sprintf('correlation_map_%d.nii',ty));
                Vo.dt = [16 0];
                Xo = zeros(Vo.dim);
                Xo (linIndx)=CR(ty,:);
                spm_write_vol(Vo,Xo);
            end;
        end;
    case 'TS_get_meants'        % Get univariately pre-whitened mean times series for each region
        % sc1_connectivity('TS_get_meants',[2 3 6 8 9 10 12 17:22],'sc2',4,'162_tessellation_hem');
        % sc1_connectivity('TS_get_meants',[2:22],'sc1',4,'162_tessellation_hem');
        sn=varargin{1}; % subjNum
        exper=varargin{2}; % Experiment string 'sc1' or 'sc2'
        glm=varargin{3}; % glmNum
        type=varargin{4}; % 'cortical_lobes','whole_brain','yeo','desikan','cerebellum','yeo_cerebellum'
        
        if glm == 6;     % Select the right directory for GLMs different from GLM4
            imagingDir=[imagingDir '_aggr'];
        end
        
        B = [];
        glmDir =fullfile(rootDir,exper,sprintf('/GLM_firstlevel_%d',glm));
        subjs=length(sn);
        
        for s=1:subjs,
            glmDirSubj=fullfile(glmDir, subj_name{sn(s)});
            load(fullfile(glmDirSubj,'SPM_light.mat'));
            %            load(fullfile(glmDirSubj,'SPM.mat'));
            T=load(fullfile(glmDirSubj,'SPM_info.mat'));
            
            % load data
            tic;
            load(fullfile(sc1Dir,regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',type))); % 'regions' are defined in 'ROI_define'
            SPM=spmj_move_rawdata(SPM,fullfile(rootDir,'sc1',imagingDir,subj_name{sn(s)}));
            
            % Get the raw data files
            V=SPM.xY.VY;
            VresMS = spm_vol(fullfile(glmDirSubj,'ResMS.nii'));
            fprintf('Moved %d\n',toc);
            % Get time series data
            tic;
            Y = region_getdata(V,R);  % Data is N x P
            resMS = region_getdata(VresMS,R);
            fprintf('Data %d\n',toc);
            
            % Spatially prewhiten and average
            Data=zeros(numel(V),numel(R));
            for r=1:numel(R), % R is the output 'regions' structure from 'ROI_define'
                Y{r}=bsxfun(@rdivide,Y{r},sqrt(resMS{r}));
                Data(:,r)=nanmean(Y{r},2);
            end;
            
            clear Y;
            
            % Redo regression
            reg_interest=[SPM.xX.iH SPM.xX.iC];
            Yfilt = SPM.xX.W*Data;
            B = SPM.xX.pKX*Yfilt;                             %-Parameter estimates
            Yres = spm_sp('r',SPM.xX.xKXs,Yfilt);             %-Residuals
            
            % Decompose betas into mean and residuals
            Z=indicatorMatrix('identity',T.cond);
            Bm = Z*pinv(Z)*B(reg_interest,:);   % Mean betas
            Br = B(reg_interest,:)-Bm;
            Yhatm = SPM.xX.xKXs.X(:,reg_interest)*Bm; %- predicted values
            Yhatr = SPM.xX.xKXs.X(:,reg_interest)*Br; %- predicted values
            filename=(fullfile(rootDir,exper,regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('ts_%s.mat',type)));
            save(filename,'Data','Yres','Yhatm','Yhatr','B');
            
            fprintf('ts saved for %s (%s) for %s \n',subj_name{sn(s)},sprintf('glm%d',glm),type);
        end
    case 'TS_get_ts'            % Get univariately pre-whitened time series for each voxel of a region
        % sc1_connectivity('TS_get_ts',[2:22],'sc2',4,'Cerebellum_grey');
        % sc1_connectivity('TS_get_ts',[2 3 6 8 9 10 12 17:22],'sc2',4,'Cerebellum_grey');
        sn=varargin{1}; % subjNum
        exper=varargin{2}; % Experiment string 'sc1' or 'sc2
        glm=varargin{3}; % glmNum
        type=varargin{4}; % 'Cerebellum_grey'
        
        B = [];
        glmDir =fullfile(rootDir,exper,sprintf('GLM_firstlevel_%d',glm));
        subjs=length(sn);
        
        for s=1:subjs,
            glmDirSubj=fullfile(glmDir, subj_name{sn(s)});
            load(fullfile(glmDirSubj,'SPM_light.mat'));
            %             load(fullfile(glmDirSubj,'SPM.mat'));
            T=load(fullfile(glmDirSubj,'SPM_info.mat'));
            
            % load data
            tic;
            load(fullfile(sc1Dir,regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',type))); % 'regions' are defined in 'ROI_define'
            SPM=spmj_move_rawdata(SPM,fullfile(rootDir,'sc1',imagingDir,subj_name{sn(s)}));
            
            % Get the raw data files
            V=SPM.xY.VY;
            VresMS = spm_vol(fullfile(glmDirSubj,'ResMS.nii'));
            fprintf('Moved %d\n',toc);
            % Get time series data
            tic;
            Y = region_getdata(V,R{1});  % Data is N x P
            resMS = region_getdata(VresMS,R{1});
            fprintf('Data %d\n',toc);
            
            % Spatially prewhiten
            Y=bsxfun(@rdivide,Y,sqrt(resMS));
            
            %             % Redo regression
            reg_interest=[SPM.xX.iH SPM.xX.iC];
            Yfilt = SPM.xX.W*Y;
            B = SPM.xX.pKX*Yfilt;                             %-Parameter estimates
            Yres = spm_sp('r',SPM.xX.xKXs,Yfilt);             %-Residuals
            
            % Decompose betas into mean and residuals
            Z=indicatorMatrix('identity',T.cond);
            Bm = Z*pinv(Z)*B(reg_interest,:);     % Mean betas
            Br = B(reg_interest,:)-Bm;
            Yhatm = SPM.xX.xKXs.X(:,reg_interest)*Bm; %- predicted values
            Yhatr = SPM.xX.xKXs.X(:,reg_interest)*Br; %- predicted values
            filename=(fullfile(rootDir,exper,regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('ts_%s.mat',type)));
            save(filename,'Yres','Yhatm','Yhatr','B');
            fprintf('ts saved for %s (%s) for %s \n',subj_name{sn(s)},sprintf('glm%d',glm),type);
        end
    case 'TS_allsubj'           % Make a structure of all cortical time series of all subject - also seperate out residual from
        sn=goodsubj;  % Take only good subjects
        glm=4; % glmNum
        type='162_tessellation_hem'; % 'cortical_lobes','whole_brain','yeo','desikan','cerebellum','yeo_cerebellum'
        session_specific = 1;
        expDir = sc1Dir;
        
        vararginoptions(varargin,{'sn','glm4','type','session_specific'});
        
        glmDir =fullfile(expDir,sprintf('/GLM_firstlevel_%d',glm));
        numSubj=length(sn);
        
        for s=1:numSubj
            fprintf('%d\n',s);
            % load condition data
            glmDirSubj=fullfile(glmDir, subj_name{sn(s)});
            load(fullfile(glmDirSubj,'SPM_light.mat'));
            T=load(fullfile(glmDirSubj,'SPM_info.mat'));
            % load data
            filename=(fullfile(expDir,regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('ts_%s.mat',type)));
            S=load(filename);
            
            reg_interest=[SPM.xX.iH SPM.xX.iC];
            % Decompose betas into mean and residuals
            if (session_specific)
                Z=indicatorMatrix('hierarchicalI',[T.sess T.cond]);
            else
                Z=indicatorMatrix('identity',T.cond);
            end;
            B(:,:,s)=S.B;
            Yres(:,:,s)=S.Yres;
            Bm = Z*pinv(Z)*S.B(reg_interest,:);   % Mean betas
            Br = S.B(reg_interest,:)-Bm;
            Yhatm(:,:,s) = SPM.xX.xKXs.X(:,reg_interest)*Bm; %- predicted values
            Yhatr(:,:,s) = SPM.xX.xKXs.X(:,reg_interest)*Br; %- predicted values
        end;
        if (session_specific)
            save(fullfile(expDir,regDir,sprintf('glm%d',glm),sprintf('ts_%s_all_se.mat',type)),'Yres','Yhatm','Yhatr','B');
        else
            save(fullfile(expDir,regDir,sprintf('glm%d',glm),sprintf('ts_%s_all.mat',type)),'Yres','Yhatm','Yhatr','B');
        end;
    case 'get_meanBeta_cerebellum' % Take the time series structure and extract mean beta for cross-session modeling
        % This stores the betas like the occurr in the GLM
        type    = '162_tessellation_hem';  % 'Cerebellum_grey';
        exper   = {'sc1','sc2'};
        sn      = goodsubj;
        glm     = 4;
        vararginoptions(varargin,{'sn','glm','type'});
        
        % Loop over experiments
        for e=1:2
            myRegDir = fullfile(rootDir,exper{e},regDir);
            for s=1:length(sn);
                fprintf('%d\n',sn(s));
                filename = (fullfile(myRegDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('ts_%s.mat',type)));
                D          = load(filename);
                glmDirSubj = fullfile(rootDir,exper{e},sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)});
                T          = load(fullfile(glmDirSubj,'SPM_info.mat'));
                
                if (strcmp(type,'162_tessellation_hem'))
                    Tcort=load(fullfile(rootDir,'sc1',regDir,'data','162_reorder.mat'));
                    Tcort=getrow(Tcort,Tcort.good==1);
                    D.B = D.B(:,Tcort.regIndx);
                    D.Yres = D.Yres(:,Tcort.regIndx);
                end;
                
                % Now generate mean estimates per session
                % Note that in the new SPM_info Instruction is coded as cond =0
                T.cond=T.cond+1;            % This is only for generating the avergaging matrix
                for se=1:2
                    X=indicatorMatrix('identity_p',T.cond.*(T.sess==se));
                    B{s,se}=pinv(X)*D.B(1:size(X,1),:);  % This calculates average beta across runs, skipping intercepts
                    ResMS{s,1} = sum(D.Yres.^2)./size(D.Yres,1); % Residual mean square for this voxel/region
                end;
            end;
            save(fullfile(myRegDir,sprintf('glm%d',glm),sprintf('mbeta_%s_all.mat',type)),'B','ResMS');
        end;
    case 'sample_DesignMatrix'
        exper=varargin{1};
        glm = varargin{2};
        wdir = fullfile(rootDir,expStr{exper},sprintf('GLM_firstlevel_%d',glm));
        load(fullfile(wdir,'s02','SPM_light.mat'));
        Sess = SPM.Sess;    % Session Structure
        X = SPM.xX.X;      % Design matrix
        wX = SPM.xX.xKXs.X; % Filtered design matrix
        save(fullfile(wdir,'SampleDesignMatrix.mat'),'Sess','X','wX');
        
    case 'cortical_covariances'                 % Covariances between cortical areas in predicted time series
        sn = goodsubj;
        glm = 4;
        type = '162_tessellation_hem';
        vararginoptions(varargin,{'sn','glm','type'});
        
        % Figure out the right ones to use
        X=load(fullfile(sc1Dir,regDir,'data','162_reorder.mat'));
        Xx=getrow(X,X.newIndx);
        Xx=getrow(Xx,Xx.good==1);
        
        % Calculate correlations
        for s=1:length(sn)
            % load Individual data
            filename=(fullfile(sc1Dir,regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('ts_%s.mat',type)));
            load(filename);
            
            Y = Yhatm + Yhatr;          % Predicted timeseries(Both mean and variable)
            predT=Y(:,Xx.regIndx2);     % Remove the bad regions
            
            % Covariance and correlation
            COV(:,:,s)=cov(predT);
            COR(:,:,s)=corrcov(COV(:,:,s));
            
            % Variance inflation factor for each
            T.SN(s,1)=s;
            T.VAR(s,:)=diag(COV(:,:,s)); % Variance of the task-related betas: How much activity could we induce?
            T.VIF(s,:)=diag(inv(COR(:,:,s)))'; % Variance inflation factor: How much does the region suffer from co-linearity?
            T.VARB(s,:)=diag(inv(COV(:,:,s))); % Combination of the last two: How well would we estimate connectivity weights from this region?
        end;
        
        % Plot group mean in flatmap
        data{1}=sqrt(mean(T.VAR,1));
        data{2}=sqrt(mean(T.VIF,1));
        data{3}=sqrt(mean(T.VARB,1));
        plottitle{1}='Variance';
        plottitle{2}='Variance Inflation Factor';
        plottitle{3}='Combined';
        
        figure;
        set(gcf,'color','w');
        ppos=1;
        for i=1:3
            subplot(3,2,ppos);
            ii=find(Xx.hem==1)';
            sc1sc2_connectivity('plot_cortex',1,data{i}(ii));
            axis off;
            tmp=get(gca,'position');
            set(gca,'position',[tmp(1) tmp(2) 0.9*tmp(3) tmp(4)]);
            title(plottitle{i})
            ppos=ppos+1;
            subplot(3,2,ppos);
            ii=find(Xx.hem==2)';
            sc1sc2_connectivity('plot_cortex',2,data{i}(ii));
            axis off;
            tmp=get(gca,'position');
            set(gca,'position',[tmp(1)-0.1 tmp(2) 1.2*tmp(3) tmp(4)]);
            ppos=ppos+1;
            colorbar;
        end
        
        set(gcf,'paperposition',[10 10 7 7])
        wysiwyg
                
    case 'get_wcon_cerebellum'           %  Gets the standardized activation maps for each session and experiment
        exper = {'sc1','sc2'};
        glm   = 4;
        yname = 'Cerebellum_grey';
        inclInstr = 0;                           % Include instruction?
        vararginoptions(varargin,{'exper','glm','yname','xname','inclInstr'});
        
        % Load the betas for all the tasks
        for e=1:2
            myRegDir = fullfile(rootDir,exper{e},regDir);
            YD{e}=load(fullfile(myRegDir,sprintf('glm%d',glm),sprintf('mbeta_%s_all.mat',yname)));
        end;
        
        % Make an integrated structure, either including instruction or
        % rest
        T=dload(fullfile(rootDir,'sc1_sc2_taskConds.txt'));
        
        % Make an integrated Structure
        T1=T;
        T2=T;
        T1.sess=ones(length(T.condNum),1)*1;
        T2.sess=ones(length(T.condNum),1)*2;
        T=addstruct(T1,T2);
        for s=1:size(YD{1}.B,1)
            Y{s}=[];
            for se=1:2
                for e=1:2
                    YD{e}.B{s,se}=[YD{e}.B{s,se}(2:end,:);zeros(1,size(YD{e}.B{s,se},2))];
                    YD{e}.B{s,se}=bsxfun(@minus,YD{e}.B{s,se},mean(YD{e}.B{s,se}));
                end;
                Y{s}=[Y{s};YD{1}.B{s,se};YD{2}.B{s,se}];
            end;
        end;
        varargout = {Y,T};
    case 'get_wcon_cortex'  % Gets the avarage wcontrast - from the t=regions_beta 
        exper = [1 2];
        glm   = 4;
        xres = 42;
        inclInstr = 0; 
        center = 1; 
        sn   = returnSubjs; 
        avrgResMS = 1;   % Averager residual MS across experiments? 
        vararginoptions(varargin,{'exper','glm','xres','inclInstr','avrgResMS'});
        
        % Load the betas for all the tasks
        Isoname =  fullfile(wbDir,'group32k',sprintf('Icosahedron-%d.32k.L.label.gii',xres));
        G=gifti(Isoname); 

        % Make an integrated structure, either including instruction or
        T=dload(fullfile(rootDir,'sc1_sc2_taskConds.txt'));
        if (inclInstr)
        
        end;
        
        % Make an integrated Structure
        T1=T;
        T2=T;
        T1.sess=ones(length(T.condNum),1)*1;
        T2.sess=ones(length(T.condNum),1)*2;
        T=addstruct(T1,T2);

        for s=1:length(sn) 
            load(fullfile(rootDir,'sc1','RegionOfInterest','data',subj_name{sn(s)},sprintf('regions_cortex.mat'))); % 'regions' are defined in 'ROI_define'
            ResMS={}; 
            voxDat={}; 
            for e=exper
                glmDirSubj = fullfile(rootDir,expStr{e},sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)});
                spmT          = load(fullfile(glmDirSubj,'SPM_info.mat'));
            
                % Load the betas 
                myRegDir = fullfile(rootDir,expStr{e},regDir,sprintf('glm%d',glm),subj_name{sn(s)});
                load(fullfile(myRegDir,'betas_cortex.mat'));        
                for h=1:2 
                    ResMS{h}(e,:)=B{h}.resMS;  % Save away residual mean square for averaging 
            
                    % Now generate mean estimates per session
                    % Note that in the new SPM_info Instruction is coded as cond =0
                    if (inclInstr) 
                        spmT.cond=spmT.cond+1;            % This is only for generating the avergaging matrix
                    end;
                    for se=1:2
                        X=indicatorMatrix('identity_p',spmT.cond.*(spmT.sess==se));
                        voxDat{h,e,se}=pinv(X)*B{h}.betasNW(1:size(X,1),:);  % This calculates average beta across runs, skipping intercepts
                    end; 
                end; 
            end;
        
        
            for h=1:2  % Look over hemispheres and condense the data
                if (avrgResMS) 
                    ResMS{h}(1,:) = mean(ResMS{h}); 
                    ResMS{h}(2,:) = ResMS{h}(1,:); 
                end; 
                
                % Figure out mapping from Nodes to voxels in region 
                N = length(R{h}.linvoxidxs); 
                MAP = nan(size(R{h}.location2linvoxindxs)); 
                for i=1:N
                    MAP(R{h}.location2linvoxindxs==R{h}.linvoxidxs(i))=i;
                end;
                numReg = max(G.cdata); % Numner of cortical regions. 
                regIndx = G.cdata(R{h}.location); 
                
                for se=1:2
                    for e=1:2 
                        voxDat{h,e,se}=[voxDat{h,e,se};zeros(1,size(voxDat{h,e,se},2))];  % Add rest 
                        voxDat{h,e,se}=bsxfun(@rdivide,voxDat{h,e,se},sqrt(ResMS{h}(h,:))); % Prewhiten 
                        voxDat{h,e,se}=bsxfun(@minus,voxDat{h,e,se},mean(voxDat{h,e,se}));  % Remove the mean from each session 
                        
                        % Summarize over cortical parcels 
                        for r=1:numReg 
                            voxIndx=(MAP(regIndx == r,:));
                            regData{h,e,se}(:,r)=nanmean(voxDat{h,e,se}(:,unique(voxIndx(~isnan(voxIndx)))),2); 
                        end; 
                    end; 
                end; 
                reg(h).region = [1:numReg]'; 
                reg(h).hemis  = ones(numReg,1)*h;
            end; 

            D{s}=[[regData{1,1,1} regData{2,1,1}];...
                [regData{1,2,1} regData{2,2,1}];...
                [regData{1,1,2} regData{2,1,2}];...
                [regData{1,2,2} regData{2,2,2}]]; 
            fprintf('Subject %d\n',sn(s)); 
        end;                 
        Reg = addstruct(reg(1),reg(2)); 
        outname = fullfile(rootDir,'sc1','connectivity_cerebellum','wcon_all',sprintf('wcon_cortex_%d.mat',xres)); 
        save(outname,'D','T','Reg');
        varargout = {D,T,Reg};
    case 'mbeta_reliability' % Get beta reliability only on common task conditions
        % In the cerebellum (Y)
        glm =4;
        vararginoptions(varargin,{'glm'});
        [X,Y,T]=sc1sc2_connectivity('get_mbeta_all','glm',glm);
        indx=find(T.overlap==1);
        T=getrow(T,indx);
        [condn,ind] = unique(T.condNumUni);
        condNames = T.condNames(ind);
        RR=[];
        for s=1:length(Y)               % Loop over subjects
            D=Y{s}(indx,:);             % Get the data
            B=indicatorMatrix('identity',T.StudyNum*2+T.sess-2);
            D=D-B*pinv(B)*D;            % Subtract the mean pattern within each study and session
            in=~isnan(sum(D));
            for c=1:length(condn)
                R.cond=condn(c);
                R.sn=s;
                i11 = find(T.condNumUni==condn(c) & T.StudyNum==1 & T.sess==1);
                i12 = find(T.condNumUni==condn(c) & T.StudyNum==1 & T.sess==2);
                i21 = find(T.condNumUni==condn(c) & T.StudyNum==2 & T.sess==1);
                i22 = find(T.condNumUni==condn(c) & T.StudyNum==2 & T.sess==2);
                R.Rw1  =     [corrN(D(i11,in)',D(i12,in)')];
                R.Rw2  =     [corrN(D(i21,in)',D(i22,in)')];
                R.Rbetween = [corrN(D(i11,in)',D(i21,in)') corrN(D(i11,in)',D(i22,in)') ...
                    corrN(D(i12,in)',D(i21,in)') corrN(D(i12,in)',D(i22,in)')];
                RR=addstruct(RR,R);
            end;
        end;
        RR.Rw = ssqrt(RR.Rw1.*RR.Rw2);
        RR.Rb = mean(RR.Rbetween,2);
        subplot(1,2,1);
        scatterplot(RR.Rw,RR.Rb,'identity');
        xlabel('Within-experiment reliability');
        xlabel('Between-experiment reliability');
        subplot(1,2,2);
        barplot(RR.cond,[RR.Rw1 RR.Rw2 RR.Rb],'leg',{'Exp1','Exp2','crossExp'});
        ylabel('Within-subject Reliability');
        set(gca,'XTickLabel',condNames,'XTickLabelRotation',60);
        varargout={RR,T};
    case 'conn_mbeta'              % Estimate connectivty model on the meanbetas
        % L2=[100 500 1000 2000 3000 4000];
        % sc1sc2_connectivity('conn_mbeta','method','ridgeFixed','lambdaL1',L2*0,'lambdaL2',L2,'name','mb4_162_ridge','exper',1);
        % sc1sc2_connectivity('conn_mbeta','method','ridgeFixed','lambdaL1',L2*0,'lambdaL2',L2,'name','mb4_162_ridge','exper',2);
        % L2=[100 1000 2000];
        % sc1sc2_connectivity('conn_mbeta','method','cplexqpL1L2','lambdaL1',L2*0,'lambdaL2',L2,'name','mb4_162_nn','exper',1);
        
        T=dload(fullfile(rootDir,'sc1_sc2_taskConds.txt'));
        
        sn=returnSubjs;
        RR= [];
        method= 'ridgeFixed';  % linRegress, ridgeFixed, nonNegExp, cplexqp, lasso, winnerTakeAll
        yname = 'Cerebellum_grey';
        xname = 'cortex_42';
        name  = 'mb4_42_nneg';  % mb4_162_ridge or ts4_yeo_nneg or mb5_162_lasso
        trainMode = 'crossed';    % training can be done in a crossed or uncrossed fashion
        lambdaL1 = 2500;
        lambdaL2 = 2500;
        numReg = NaN; 
        overwrite = 0;  % Overwrite old results or add to existing structure?
        exper  = 1;     % Task set to train
        inclInstr = 0;   % For future use: include Instruction
        vararginoptions(varargin,{'sn','xname','type','exper','method','glm',...
            'name','trainMode','lambdaL1','lambdaL2','incInstr','overwrite','numReg'});
        
        % Get data
        [Y,S]=sc1sc2_connectivity('get_wcon_cerebellum','inclInstr',inclInstr);
        fname = fullfile(rootDir,'sc1','connectivity_cerebellum','wcon_all',sprintf('wcon_%s.mat',xname)); 
        load(fname); 
        
        % Save the fit
        outDir = fullfile(rootDir,expStr{exper},connDir,name);
        dircheck(outDir);
        
        % Find the data that we want to use for fitting the connectivity
        % model
        SI1 = find(ismember(S.StudyNum,exper) & S.sess==1);
        SI2 = find(ismember(S.StudyNum,exper) & S.sess==2);
        trainXindx=[SI1;SI2];
        switch (trainMode)
            case 'crossed'
                trainYindx=[SI2;SI1];
            case 'uncrossed'
                trainYindx=[SI1;SI2];
        end;
        
        % Estimate the model and store the information attached to it
        for s=1:length(sn)
            % Find relevant file
            fprintf('%d\n',sn(s));
            outName=fullfile(outDir,sprintf('%s_%s.mat',name,subj_name{sn(s)}));
            if (exist(outName) && overwrite==0);
                RR=load(outName);
            else
                RR=[];
            end;
            
            % Get data
            ss = find(sn(s)==returnSubjs); 
            xx = D{ss}(trainXindx,:);
            yy = Y{ss}(trainYindx,:);
            
            % Run for all regularisation parameters
            for l=1:length(lambdaL1)
                R=[]; 
                [W,fR2m,fRm] = sc_connect_fit(yy,xx,method,...
                    'lambda',[lambdaL1(l) lambdaL2(l)],'numReg',numReg);
                for i=1:size(W,3)
                    R.SN(i,1)     = sn(s);
                    R.W{i,1}      = W(:,:,i);
                    R.lambda(i,:) = [lambdaL1(l) lambdaL2(l)];
                    R.method{i,1} = method;
                    R.incInstr(i,1) = inclInstr;
                    R.trainMode{i,1} = trainMode;
                    R.xname{i,1}  = xname;
                    R.numReg(i,1) = i; 
                    R.fR2m(i,1)   = fRm(i); 
                    R.fRm(i,1)    = fR2m(i); 
                end; 
                RR=addstruct(RR,R);
            end;
            save(outName,'-struct','RR');
        end;  
        varargout={RR}; 
    case 'run_connB'
        L1=[10 25 50 100];
        L2=[250 500 1000 1750 3500]; 
        sc1sc2_connectivity('conn_mbeta','method','cplexqpL1L2','lambdaL1',L1,'lambdaL2',L1*0,'name','mb4_162_nn','exper',1);
        sc1sc2_connectivity('conn_mbeta','method','cplexqpL1L2','lambdaL1',L1,'lambdaL2',L1*0,'name','mb4_162_nn','exper',2);
        sc1sc2_connectivity('conn_mbeta','method','cplexqpL1L2','lambdaL1',L2*0,'lambdaL2',L2,'name','mb4_162_nn','exper',1);
        sc1sc2_connectivity('conn_mbeta','method','cplexqpL1L2','lambdaL1',L2*0,'lambdaL2',L2,'name','mb4_162_nn','exper',2);
        sc1sc2_connectivity('conn_mbeta','method','winnerTakeAll_nonNeg','name','mb4_362_wtan','xname','cortex_362','exper',1,'overwrite',1,'numReg',1);
        sc1sc2_connectivity('conn_mbeta','method','nonNegStepwise','name','mb4_362_nnStep','xname','cortex_362','exper',1,'overwrite',1,'numReg',7);
        sc1sc2_connectivity('conn_mbeta','method','nonNegStepwise','name','mb4_362_nnStep','xname','cortex_362','exper',2,'overwrite',1,'numReg',7);
    
    case 'conn_ts'                 % Run encoding model on time series
        % sc1_connectivity('conn_ts',[2:22],'method','winnerTakeAll_nonNeg','name','glm4_162_nnWTA','lambdaL1',0,'lambdaL2',0);
        % sc1_connectivity('conn_ts',goodsubj,'method','winnerTakeAll','name','glm4_162_WTA','lambdaL1',0,'lambdaL2',0);
        % sc1_connectivity('conn_ts',[2:22],'method','ridgeFixed','name','glm4_162_ridge','lambdaL1',[0 0 0 0 0 0],'lambdaL2',[0 500 1000 2500 5000 7500]);
        % sc1_connectivity('conn_ts',[2:22],'method','cplexqpL1L2','name','glm4_162_L1L2c','lambdaL1',[2500 7500 7500],'lambdaL2',[2500 7500 2500]);
        sn=goodsubj;%varargin{1};
        glm=4;                          % usually 4
        method= 'ridgeFixed';           % cplexqpL1L2, linRegress, nonNegExp, cplexqp, lasso, winnerTakeAll
        yname = 'cerebellum_grey';
        xname = '162_tessellation_hem';
        name  = 'glm4_162_ridgeFixed';  % Name of the output
        type  = {'Y','Yhat'};           % Run with three Y
        exper   = 'sc1';
        lambdaL1 = 0;                   %[2500 7500 7500];    % Initially set to 0
        lambdaL2 = 2500;                %[2500 7500 2500];    % Inititally select one value
        partial = [];
        vararginoptions({varargin{2:end}},{'xname','type','method','glm','partial','name','lambdaL1','lambdaL2'});
        
        outDir = fullfile(rootDir,exper,connDir,name);
        dircheck(outDir);
        myRegDir = fullfile(rootDir,exper,regDir);
        
        % Build matrix to remove block mean
        SPM=load(fullfile(rootDir,exper,'GLM_firstlevel_4','SampleDesignMatrix.mat'));
        for b=1:16
            B(SPM.Sess(b).row,b)=1;
        end;
        
        for s=sn
            TT=[];
            T=[];
            
            % Load X
            file = fullfile(myRegDir,sprintf('glm%d',glm),subj_name{s},sprintf('ts_%s.mat',xname));
            XX=load(file);
            
            if (strcmp(xname,'162_tessellation_hem'))
                Tcort=load(fullfile(myRegDir,'data','162_reorder.mat'));
                Tcort=getrow(Tcort,Tcort.good==1);
                XX.Yhatm = XX.Yhatm(:,Tcort.regIndx,:);
                XX.Yhatr = XX.Yhatr(:,Tcort.regIndx,:);
                XX.Yres =  XX.Yres(:,Tcort.regIndx,:);
            end;
            
            if (strcmp(xname,'yeo'))
                XX.Yhatm = XX.Yhatm(:,2:18,:);
                XX.Yhatr = XX.Yhatr(:,2:18,:);
                XX.Yres =  XX.Yres(:,2:18,:);
            end;
            
            XX.Yhat = XX.Yhatm + XX.Yhatr;
            XX.Y    = XX.Yres + XX.Yhat;
            
            
            if (~isempty(partial))
                PP=load(fullfile(rootDir,exper,regDir,regDir,'glm4',subj_name{sn(s)},sprintf('ts_%s.mat',partial)));
                PP.Yhat = PP.Yhatm + PP.Yhatr;
                PP.Y = PP.Yhat + PP.Yres;
            end;
            
            % Load Y
            Y=load(fullfile(myRegDir,'glm4',subj_name{s},sprintf('ts_%s.mat',yname)));
            Y.Yhat = Y.Yhatm + Y.Yhatr;
            Y.Y    = Y.Yhat  + Y.Yres;
            
            % Loop over the different types
            for t=1:length(type);
                yy = Y.(type{t});
                xx = XX.(type{t});
                
                % remove mean
                yy = yy - B*pinv(B)*yy;
                xx = xx - B*pinv(B)*xx;
                
                % Normalize to unity
                %ySS = sum((Y.Y).^2);
                % yy = bsxfun(@rdivide,yy,sqrt(ySS));
                %xSS = sum((X.Y).^2);
                % xx = bsxfun(@rdivide,X.Yhat,sqrt(xSS));
                % X.Yres = bsxfun(@rdivide,X.Yres,sqrt(xSS));
                
                % Regress out the partial ROI
                if (~isempty(partial))
                    pp = PP.(type{t});
                    pp = pp - B*pinv(B)*pp;
                    yy = yy - pp*pinv(pp)*yy;  % Subtract from cerebellar data only
                end;
                
                for l=1:length(lambdaL1);
                    T.SN = s;
                    T.typeNum = t;
                    T.type    = {type{t}};
                    T.lambda  = [lambdaL1(l) lambdaL2(l)];
                    
                    fprintf('Fitting s%2.2d %s lambdaL1 = %2.2f lambdaL2 = %2.2f\n',s,type{t},lambdaL1(l),lambdaL2(l));
                    [W,T.fR2m,T.fRm,T.fR2v,T.fRv]=sc_connect_fit(yy,xx,method,'lambda',[lambdaL1(l) lambdaL2(l)]);
                    T.W = {W};
                    TT=addstruct(TT,T);
                end;
            end;
            outName=fullfile(outDir,sprintf('%s_s%2.2d.mat',name,s));
            save(outName,'-struct','TT');
            fprintf('encode model: cerebellar voxels predicted for %s \n',subj_name{s});
        end
        
    case 'numberOfConditions'      % Build predictive models on a random assortment of conditions - evaluate them on the other experiment
        T = dload(fullfile(rootDir,'sc1_sc2_taskConds.txt'));
        RR=[];
        for nc=2:29
            for n=1:50
                CondSet = sample_wor([1:29],nc,1);
                M=sc1_connectivity('conn_mbeta','trainMode','uncrossed','subset',T.StudyNum==1 & ismember(T.condNum,CondSet));
                R=sc1_connectivity('evaluate',M,'subset',T.StudyNum==2,'splitby',ismember(T.condNumUni,CondSet));
                vec = ones(length(R.SN),1);
                R.exp      = 2*vec;
                R.numSet  = n*vec;
                R.numCond = nc*vec;
                R.CondSet = repmat(ismember([1:29],CondSet),length(R.SN),1);
                RR=addstruct(RR,R);
                R=sc1_connectivity('conn_evaluate',M,'subset',T.StudyNum==1,'splitby',ismember(T.condNumUni,CondSet));
                vec = ones(length(R.SN),1);
                R.exp      = vec;
                R.numSet  = n*vec;
                R.numCond = nc*vec;
                R.CondSet = repmat(ismember([1:29],CondSet),length(R.SN),1);
                RR=addstruct(RR,R);
            end;
        end;
        save(fullfile(sc1Dir,'connectivity_cerebellum','evaluation','numberOfConditions_glm4_1_2.mat'),'-struct','RR');
        varargout={RR};
    case 'singleTaskPlot'
        R=varargin{1};
        T=dload(fullfile(rootDir,'sc1_sc2_taskConds.txt'));
        T=getrow(T,T.StudyNum==2);
        [D.Rcv,D.cond]=pivottable(R.cond,[],R.Rcv,'mean');
        [D.Ry,D.cond]=pivottable(R.cond,[],R.Ry,'mean');
        [D.Rp,D.cond]=pivottable(R.cond,[],R.Rp,'mean');
        barplot(D.cond,D.Rcv./sqrt(D.Ry.*D.Rp));
        set(gca,'XTickLabel',T.condNames,'XTickLabelRotation',70);
        keyboard;
        varargout={R,D};
    case 'map_cortex'              % Map of where projections come from
        sn=varargin{1};        % [2:22]
        name = varargin{2};
        
        glm=4;       % usually 4
        parcel = 'yeo17';
        type  = {'Yres','Yhat'};
        vararginoptions({varargin{3:end}},{'type','glm','parcel'});
        
        lambda=0;
        encodeDirSp = fullfile(encodeDir,sprintf('glm%d',glm),'encode_4');
        
        for s=1:length(sn)
            T=load(fullfile(encodeDirSp,sprintf('encode_%s_s%2.2d.mat',name,sn(s))));
            for t=1:2
                data{t}(s,:) = mean(T.W{t}');
            end;
        end;
        for t=1:2
            md =mean(data{t},1);
            max(md)
            for h=1:2
                subplot(2,2,h+(t-1)*2);
                switch(parcel)
                    case '162'
                        dat = md([1:142]+(h-1)*142);
                        cscale = [0 0.005];
                    case 'yeo17'
                        dat = md;
                        cscale = [0 0.05];
                end;
                sc1_imana_jd('plot_cortex',h,dat,'parcel',parcel,'cscale',cscale);
            end;
        end;
    case 'plot_cortex'
        h=varargin{1};   %
        data=varargin{2};
        
        cscale = [];
        parcel = '162';
        
        vararginoptions({varargin{3:end}},{'cscale','parcel'});
        switch (parcel)
            case '162'
                D=caret_load(fullfile(sc1Dir,caretDir,'fsaverage_sym',hemName{h},sprintf([hem{h} '.tessel162_new.metric'])));
                Tcort=load(fullfile(sc1Dir,regDir,'data','162_reorder.mat'));
                Tcort=getrow(Tcort,Tcort.good==1 & Tcort.hem==h);
                indx = Tcort.regIndx-158*(h-1);
            case 'yeo17'
                D=caret_load(fullfile(sc1Dir,caretDir,'fsaverage_sym',hemName{h},sprintf([hem{h} '.yeo17.paint'])));
                indx=[2:18]';
        end;
        
        FData=nan(size(D.data));
        for i=1:length(data)
            FData(D.data==indx(i),1)=data(i);
        end;
        
        coord=fullfile(sc1Dir,caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
        topo=fullfile(sc1Dir,caretDir,'fsaverage_sym',hemName{h},[hem{h} '.CUT.topo']);
        if (h==1)
            caret_plotflatmap('coord',coord,'topo',topo,'data',FData,'xlims',[-200 150],'ylims',[-140 140],'cscale',cscale);
        else
            caret_plotflatmap('coord',coord,'topo',topo,'data',FData,'xlims',[-150 200],'ylims',[-140 140],'cscale',cscale);
        end;
    case 'map_cerebellum'          % Maps certain stats from an individual to SUIT and then to flat map
        data    = varargin{1};
        sn      = varargin{2};
        
        % Determine the voxels we want to resample in SUIT space
        V=spm_vol(fullfile(sc1Dir,suitDir,'anatomicals','cerebellarGreySUIT.nii'));
        X= spm_read_vols(V);
        grey_threshold = 0.1; % gray matter threshold
        linIn1=find(X>grey_threshold);
        [i1,j1,k1]= ind2sub(V.dim,linIn1');
        [x1,y1,z1]= spmj_affine_transform(i1,j1,k1,V.mat);
        
        % Determine voxel locations from the original ROI
        load(fullfile(sc1Dir,regDir,'data',subj_name{sn},'regions_cerebellum_grey.mat')); % 'regions' are defined in 'ROI_define'
        Vmask = spm_vol(fullfile(sc1Dir,suitDir,'glm4',subj_name{sn},'maskbrainSUITGrey.nii'));
        [i3,j3,k3]=spmj_affine_transform(R{1}.data(:,1),R{1}.data(:,2),R{1}.data(:,3),inv(Vmask.mat));
        linIn3 = sub2ind(Vmask.dim,round(i3),round(j3),round(k3));
        
        % transform SUIT coords into anatomical space of the individual
        flowfield = fullfile(sc1Dir,suitDir,'anatomicals',subj_name{sn},'u_a_c_anatomical_seg1.nii');
        affine    = fullfile(sc1Dir,suitDir,'anatomicals',subj_name{sn},'Affine_c_anatomical_seg1.mat');
        [Def,Aff]=spmdefs_get_dartel(flowfield,affine);
        [x2,y2,z2]=spmdefs_transform(Def,Aff,x1,y1,z1);
        [i2,j2,k2]=spmj_affine_transform(x2,y2,z2,inv(Vmask.mat));
        
        % resample the weights into SUIT space
        for r=1:size(data,1)
            Vout=Vmask;
            Vout.dat = zeros(Vout.dim);
            Vout.dat(linIn3)=data(r,:);
            Vout.dt = [64 0];
            Vout.pinfo = [1 0 0]';
            DataSUIT(r,:)=spm_sample_vol(Vout,i2,j2,k2,1);
            V.dat = zeros(V.dim);
            Vres(r)=V;
            Vres(r).dat(linIn1)=DataSUIT(r,:);  % Offset by one to account for 1 being medial wall
            Vres(r).fname=sprintf('data_%2.2d.nii',r);
            Vres(r).pinfo=[1 0 0]';
        end;
        
        % Now map to surface-based representation
        D = suit_map2surf(Vres,'stats','mean');
        varargout={D,Vres};
    case 'evaluate' % Evaluates predictive performance of connectivity model
        % M: is the model structure with the connectivity weights in it (M.W)
        % 'subset': what data should the evaluation be based upon?
        % 'splitby': by which variable should the evaluation be split?
        % 'meanSub': Mean pattern subtraction before evaluation?
        M = varargin{1};        % This is a structure with the fitted data
        subset = [];
        splitby= [];
        yname = 'Cerebellum_grey';
        xname = 'cortex_42';
        inclInstr= 0; 
        meanSub = 0; % Mean pattern subtraction before prediction?
        vararginoptions(varargin(2:end),{'subset','splitby','meanSub','xname','yname'});
        
        % Get all the mean betas and prepare the evaulation data
        [Y,S]=sc1sc2_connectivity('get_wcon_cerebellum','inclInstr',inclInstr);
        fname = fullfile(rootDir,'sc1','connectivity_cerebellum','wcon_all',sprintf('wcon_%s.mat',xname)); 
        load(fname);
        X=D; 
        T=dload(fullfile(rootDir,'sc1_sc2_taskConds.txt'));
        if (isempty(subset))
            subset = T.StudyNum==2; 
        end; 
        S.subset= [subset;subset];
        if (isempty(splitby))
            splitby = ones(length(T.StudyNum),1);
        end;
        S.splitby=[splitby;splitby];  % Double for the two conditions
        sS = getrow(S,S.subset);
        splits = unique(sS.splitby);
        
        % Loop over models and evaulate
        RR=[];
        numModels = size(M.SN,1);
        for m=1:numModels
            s=find(returnSubjs==M.SN(m));
            goodindx = sum(abs(Y{s}(:,:)))>0 & ~isnan(sum(Y{s}(:,:))) & ~isnan(sum(M.W{m}));
            for sp =1:length(splits);
                if(meanSub)
                    for sess=1:2
                        indx=find(S.subset & S.sess==sess & S.splitby==splits(sp));
                        X{s}(indx,:)=bsxfun(@minus,X{s}(indx,:),mean(X{s}(indx,:)));
                        Y{s}(indx,:)=bsxfun(@minus,Y{s}(indx,:),mean(Y{s}(indx,:)));
                    end;
                end;
                testAindx=[find(S.subset & S.sess==1 & S.splitby==splits(sp));...
                    find(S.subset & S.sess==2 & S.splitby==splits(sp))];
                testBindx=[find(S.subset & S.sess==2 & S.splitby==splits(sp));...
                    find(S.subset & S.sess==1 & S.splitby==splits(sp))];
                predY   = X{s}(testAindx,:)*M.W{m};                    % Predicted Y using crossvalidation
                predYnc = X{s}(testBindx,:)*M.W{m};                    % Predicted Y not crossvalidated
                SSP   = sum(sum(predY(:,goodindx).^2));             % Sum of square of predictions
                SSY   = sum(sum(Y{s}(testBindx,goodindx).^2));               % Sum of squares of data
                SSCp  = sum(sum(predY(:,goodindx).*predYnc(:,goodindx))); % Covariance of PRedictions
                SSCy  = sum(sum(Y{s}(testAindx,goodindx).*Y{s}(testBindx,goodindx)));  % Covariance of Y's
                SSCn  = sum(sum(predYnc(:,goodindx).*Y{s}(testBindx,goodindx)));   % Covariance of non-cross prediction and data
                SSCc  = sum(sum(predY(:,goodindx).*Y{s}(testBindx,goodindx)));   % Covariance of cross prediction and data
                R.SN    = M.SN(m);
                R.lambda = M.lambda(m,:);
                R.method = M.method(m);
                R.trainMode = M.trainMode(m);
                R.xname = M.xname(m);
                R.Rcv   = SSCc ./ sqrt(SSY.*SSP); % Double-Crossvalidated predictive correlation
                R.Rnc   = SSCn ./ sqrt(SSY.*SSP); % Not double-crossvalidated predictive correlation
                R.Ry    = SSCy ./ SSY;            % Reliability of data: noise ceiling
                R.Rp    = SSCp ./ SSP;            % Reliability of prediction
                R.split = splits(sp);
                % Calucate Sparseness measures
                Ws = sort(abs(M.W{m}));
                Wss= bsxfun(@rdivide,Ws,sum(Ws));% standardized coefficients (to the sum overall
                R.spIdx = nanmean(Wss(end,:)); % Largest over the sum of the others
                N=size(Wss,1);
                w=(N-[1:N]'+0.5)/N;
                ginni = 1-2*sum(bsxfun(@times,Wss,w));
                R.ginni = nanmean(ginni); % Ginni index: mean over voxels
                R.numReg = M.numReg(m); 
                RR = addstruct(RR,R);
            end;
        end;
        varargout={RR,Y{s}(testBindx,:),predY};
    case 'evaluate_all'                     % Evaluates different sets of Connnectivity models
        whatAna=varargin{1};                % Which particular models / experiments / modes do you want to compare
        outname = whatAna;
        xnames = {'cortex_42','cortex_162','cortex_362'};
        sn = returnSubjs;

        D=dload(fullfile(rootDir,'sc1_sc2_taskConds.txt'));
        switch(whatAna)
            case 'mb4_nnStep_all'
                name   = {'mb4_42_nnStep';'mb4_42_nnStep';'mb4_162_nnStep';'mb4_162_nnStep';'mb4_362_nnStep'};
                traindata = [1 2 1 2 1]';
                subset    = [D.StudyNum==2 D.StudyNum==1 D.StudyNum==2 D.StudyNum==1 D.StudyNum==2]; % Evaluate on the other experiment
                xnIn      = [1 1 2 2 3]; % use appropriate tesselation 
                ysplit    = ones(length(D.StudyNum),1); % No splitting
            case 'mb4_wtan_all'
                name   = {'mb4_42_wtan';'mb4_42_wtan';'mb4_162_wtan';'mb4_162_wtan';'mb4_362_wtan';'mb4_362_wtan';};
                traindata = [1 2 1 2 1 2]';
                subset    = [D.StudyNum==2 D.StudyNum==1 D.StudyNum==2 D.StudyNum==1 D.StudyNum==2 D.StudyNum==1]; % Evaluate on the other experiment
                xnIn      = [1 1 2 2 3 3]; % use appropriate tesselation 
                ysplit    = ones(length(D.StudyNum),1); % No splitting
            case 'mb4_nneg_all'
                name   = {'mb4_162_nneg';'mb4_162_nneg'};
                traindata = [1 2 1 2 1 2]';
                subset    = [D.StudyNum==2 D.StudyNum==1]; % Evaluate on the other experiment
                xnIn      = [2 2 ]; % use appropriate tesselation 
                ysplit    = ones(length(D.StudyNum),1); % No splitting
        end;
        RR=[];
        for i=1:length(name)
            for s=sn
                fprintf('%d\n',s);
                
                % Load the connectivity weights
                M = load(fullfile(rootDir,expStr{traindata(i)},'connectivity_cerebellum',name{i},sprintf('%s_%s.mat',name{i},subj_name{s})));
                R=  sc1sc2_connectivity('evaluate',M,'xname',xnames{xnIn(i)},'subset',subset(:,i),'splitby',ysplit);
                n = length(R.SN);
                R.traindata = ones(n,1)*traindata(i);
                RR=addstruct(RR,R);
            end;
        end;
        save(fullfile(rootDir,'sc1','connectivity_cerebellum','evaluation',sprintf('eval_%s.mat',outname)),'-struct','RR');
        varargout={RR};
    case 'make_obsPredicted'
        mapping = 'glm4_162_nn';
        sn = varargin{1};
        T = dload(fullfile(rootDir,'sc1_sc2_taskConds.txt'));
        M = load(fullfile(rootDir,'sc1','connectivity_cerebellum',mapping,sprintf('%s_%s.mat',mapping,subj_name{sn})));
        [R,Y,Yp]=sc1_connectivity('evaluate',getrow(M,4),'subset',T.StudyNum==2,'splitby',[]);
        S=getrow(T,T.StudyNum==2);
        S= addstruct(S,S); % Two sessions
        
        % Average patterns across runs
        X = indicatorMatrix('identity',S.condNum);
        Y=pinv(X)*Y;
        Yp=pinv(X)*Yp;
        
        % Map the patterns to
        Data = [Y;Yp];
        % Resorting of data
        index = reshape([1:64],32,2)';
        index = index(:);
        Data = Data(index,:);
        S=getrow(S,index);
        [DA,V] = sc1_connectivity('map_cerebellum',Data,sn);
        DA(isnan(DA))=0;
        M=caret_struct('metric','data',DA);
        M.column_name = S.condNames';
        filename = fullfile(sc1Dir,caretDir,'suit_flat','glm4',sprintf('obsPred_s%2.2d.metric',sn));
        coordfile = fullfile(sc1Dir,caretDir,'suit_flat','FLAT.coord.gii');
        topofile  = fullfile(sc1Dir,caretDir,'suit_flat','CUT.topo.gii');
        caret_save(filename,M);
        sfilename = caret_smooth(filename,'coord',coordfile,'topo',topofile);
    case 'evaluate_crossval_test'           % Checks whether the cross-session crossvalidation works appropriately
        P1 = 20;   % Areas cortex
        P2 = 1000; % voxels cerebellum
        K = 30;            %conditions
        sigma_s_1 = 1;    % Systematic signal in the cortex
        sigma_s_2 = 0;    % Additive signal in the cerebellum
        sigma_n_1 = 1;    % Systematic signal in the cortex
        sigma_n_2 = 1;  % Additive noise in the cerebellum
        sigma_w_s = 1;    % Connectivity weights for the signal
        sigma_w_n = 0;    % Connectivity weights for the noise
        Y.sess = kron([1;2],ones(K,1));
        Y.run = kron([1;2],ones(K,1));
        Y.cond = kron([1;1],[1:K]');
        T.SN=2;
        RR=[];
        for n=1:100
            WS = normrnd(0,sigma_w_s,P1,P2); WS(WS<0)=0; % non-negative connectivity weights for noise
            WN = normrnd(0,sigma_w_n,P1,P2); WN(WN<0)=0; % non-negative connectivity weights for signal
            
            N1  = normrnd(0,1,2*K,P1);
            N2  = N1*WN + normrnd(0,1,2*K,P2);
            S1  = normrnd(0,sigma_s_1,K,P1);
            S2  = S1*WS + normrnd(0,sigma_s_2,K,P2);
            S1=[S1;S1];
            S2=[S2;S2];
            
            T.WC{1} = (S1'*S1)\S1'*S2;  % Weights determined only on noise processes
            T.WC{1} = (N1'*N1)\N1'*N2;  % Weights determined only on noise processes
            T.lambda = 0;
            T.type = {[1]};
            T.typeNum = 1;
            
            Y.data = N2+S2;
            X.B    = N1+S1;
            R=sc1_imana_jd('encoding_crossval_evaluate',T,'Xdata',X,'Ydata',Y,'xname','other');
            RR=addstruct(RR,R);
        end;
        myboxplot([],[RR.Rcv RR.Rnc RR.R]);
        set(gca,'YLim',[-0.2 1]);
        drawline(0,'dir','horz');
        keyboard;
    case 'evaluate_plot_numberOfConditions'
        T=load(fullfile(sc1Dir,'connectivity_cerebellum','evaluation','numberOfConditions_glm4_1_2.mat'));
        lineplot(T.numCond,T.Rcv./abs(ssqrt(T.Ry.*T.Rp)),'subset',T.exp==2,'split',T.split,'style_thickline',...
            'linecolor',{'b','b'},'markercolor',{'b','b'},'linestyle',{'-',':'},'errorcolor',{'b','b'})
        set(gca,'YLim',[0 0.7]);
        set(gcf,'PaperPosition',[2 2 5 4]);
        
        drawline(0.319,'dir','horz');
        wysiwyg;
        
    case 'evaluate_plot_method'   % Does sc1 and sc2 with BM or not
        name = {'eval_mb4_162_all.mat'};
        yfield = 'Rcv';
        xfield = 'ginni';
        
        T=[];
        for n=1:length(name)   % Loop over files
            TT=load(fullfile(sc1Dir,'connectivity_cerebellum','evaluation',name{n}));
            T=addstruct(T,TT);
        end;
        
        [methStr,~,T.methNum] = unique(T.method);
        [lambda,b,T.lamCat]=unique([T.lambda],'rows');
        % Do an xyplot normalized to an upper noise ceiling, which only is
        % determined by the relability of the data
        xyplot(T.(xfield),T.(yfield)./sqrt(T.Ry),T.lamCat,'split',T.methNum,'style_thickline','leg',methStr,'subset',T.traindata==2);
        
        % Now determine a lower noise ceiling.... This is determined by the
        % reliability of the prediction, which is model dependent. We here
        % use the average across models
        noiseCeil = mean(T.Rp);
        drawline(noiseCeil,'dir','horz');
        set(gca,'YLim',[0.3 1]);
        set(gcf,'PaperPosition',[2 2 6 6]);
        varargout={T};
        % wysiwyg;
    case 'evaluate_plot_tessel'   % Does sc1 and sc2 with BM or not
        name = {'eval_mb4_nnStep_all.mat','eval_mb4_wtan_all.mat','eval_mb4_nneg_all.mat'};
        yfield = 'Rcv';
        xfield = 'spIdx';
        
        T=[];
        for n=1:length(name)   % Loop over files
            TT=load(fullfile(sc1Dir,'connectivity_cerebellum','evaluation',name{n}));
            T=addstruct(T,TT);
        end;
        
        [methStr,~,T.methNum] = unique(T.method);
        [lambda,b,T.lamCat]=unique([T.lambda],'rows');
        for i=1:length(T.SN) 
            [~,T.numTessel(i,1)] = textscan(T.xname{i},'%s%d','delimiter','_'); 
        end; 
        % Do an xyplot normalized to an upper noise ceiling, which only is
        % determined by the relability of the data
        xyplot(T.(xfield),T.(yfield)./sqrt(T.Ry),T.lamCat,'split',T.methNum,'style_thickline','leg',methStr,'subset',T.traindata==2);
        
        % Now determine a lower noise ceiling.... This is determined by the
        % reliability of the prediction, which is model dependent. We here
        % use the average across models
        noiseCeil = mean(T.Rp);
        drawline(noiseCeil,'dir','horz');
        set(gca,'YLim',[0.3 1]);
        set(gcf,'PaperPosition',[2 2 6 6]);
        varargout={T};
        % wysiwyg;
    case 'TS_subspace_overlap'   % Determine relative eigenvalues after projection
        load(fullfile(regDir,'glm4','Covariance_by_session.mat'));
        [numReg,~,numSubj,numSess]=size(C{1});
        type ={'Yhatm','Yres'};
        typeIn = [1 3];
        D=[];
        for i=1:numSubj
            for k=1:2 % Sessions
                [T.ev(1,:),T.ev(5,:)]=sc1_subspace_overlap(C{typeIn(1)}(:,:,i,k),C{typeIn(1)}(:,:,i,3-k));
                [T.ev(2,:),T.ev(6,:)]=sc1_subspace_overlap(C{typeIn(2)}(:,:,i,k),C{typeIn(2)}(:,:,i,3-k));
                T.ev(3,:)=sc1_subspace_overlap(C{typeIn(2)}(:,:,i,k),C{typeIn(1)}(:,:,i,3-k));
                T.ev(4,:)=sc1_subspace_overlap(C{typeIn(1)}(:,:,i,k),C{typeIn(2)}(:,:,i,3-k));
                T.session1=[k;k;k;k;k;k];
                T.session2=[3-k;3-k;3-k;3-k;k;k];
                T.subj=ones(6,1)*i;
                T.data1 = [1;2;2;1;1;2];
                T.data2 = [1;2;1;2;1;2];
                T.type  = [1;2;3;4;5;6];
                D=addstruct(D,T);
            end;
        end;
        CAT.linestyle={':';':';'-';'-';'--';'--'};
        CAT.linecolor={'r';'b';'r';'b';'r';'b'};
        
        D.cev  = cumsum(D.ev,2);
        D.cevn = bsxfun(@rdivide,D.cev,D.cev(:,end));  % Normalized cumulative eigenvalues
        D.evn  = bsxfun(@rdivide,D.ev,D.cev(:,end));   % Normalized eigenvalues
        A      = 1-[zeros(size(D.cev,1),1) D.cevn(:,1:end-1)];
        D.eevn = bsxfun(@rdivide,A,[numReg:-1:1]);
        
        P = 10;
        subplot(1,2,1);
        traceplot([1:P],D.cevn(:,1:P),'split',D.type,'leg','auto','CAT',CAT);
        subplot(1,2,2);
        traceplot([1:P],D.evn(:,1:P),'split',D.type,'subset',D.type<3);
        hold on;
        traceplot([1:P],D.eevn(:,1:P),'split',D.type,'subset',D.type<3,'linestyle',':');
        hold off;
        set(gcf,'PaperPosition',[2 2 8 4]);
        wysiwyg;
        varargout={D};
    case 'crossval_space_overlap_simulate'   % Check out a couple of different scenarious
        P=200;
        N=1000;
        sigm=0.5;
        eigC=zeros(1,P);
        eigX=zeros(1,P);eigX(1:3)=[100 30 10];
        eigY=zeros(1,P);eigY(1:30)=1;
        D=[];
        for n=1:10
            A=randn(P,P);
            [Wc,l]=eig(A*A');
            C=normrnd(0,1,N,P)*diag(sqrt(eigC))*Wc;
            A=randn(P,P);
            [Wx,l]=eig(A*A');
            X=C+normrnd(0,1,N,P)*diag(sqrt(eigX))*Wx;
            X1=X+normrnd(0,sigm,N,P);
            X2=X+normrnd(0,sigm,N,P);
            
            A=randn(P,P);
            [Wy,l]=eig(A*A');
            Y=C+normrnd(0,1,N,P)*diag(sqrt(eigY))*Wy;
            Y1=Y+normrnd(0,sigm,N,P);
            Y2=Y+normrnd(0,sigm,N,P);
            Cx1=cov(X1);
            Cx2=cov(X2);
            Cy1=cov(Y1);
            Cy2=cov(Y2);
            
            T.ev(1,:)=sc1_subspace_overlap(Cx2,Cx1);
            T.ev(2,:)=sc1_subspace_overlap(Cy2,Cy1);
            T.ev(3,:)=sc1_subspace_overlap(Cy1,Cx1);
            T.ev(4,:)=sc1_subspace_overlap(Cx1,Cy1);
            T.ev(5,:)=sc1_subspace_overlap(Cx1,Cx1);
            T.ev(6,:)=sc1_subspace_overlap(Cy1,Cy1);
            T.data1 = [1;2;2;1;1;2];
            T.data2 = [1;2;1;2;1;2];
            T.type  = [1;2;3;4;5;6];
            D=addstruct(D,T);
        end;
        CAT.linestyle={':';':';'-';'-';'--';'--'};
        CAT.linecolor={'r';'b';'r';'b';'r';'b'};
        D.cev  = cumsum(D.ev,2);
        D.cevn = bsxfun(@rdivide,D.cev,D.cev(:,end));  % Normalized cumulative eigenvalues
        D.evn  = bsxfun(@rdivide,D.ev,D.cev(:,end));   % Normalized eigenvalues
        A      = 1-[zeros(size(D.cev,1),1) D.cevn(:,1:end-1)];
        D.eevn = bsxfun(@rdivide,A,[P:-1:1]);
        
        subplot(1,2,1);
        traceplot([1:P],D.cevn(:,1:P),'split',D.type,'leg','auto','CAT',CAT);
        subplot(1,2,2);
        traceplot([1:P],D.evn(:,1:P),'split',D.type,'subset',D.type<3);
        hold on;
        traceplot([1:P],D.eevn(:,1:P),'split',D.type,'subset',D.type<3,'linestyle',':');
        hold off;
        set(gcf,'PaperPosition',[2 2 8 4 ]);
        wysiwyg;
        varargout={D};
    case 'make_reorder_mat' % consolidates the by-lobe reordering and removal of medial wall
        DD=[];
        for h=1:2
            P=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.tessel162.paint']));
            PL=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.lobes.paint']));
            PD=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.desikan.paint']));
            PY=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.Yeo17.paint']));
            D.lobe=pivottable(P.data,[],PL.data,'mode');
            D.desikan=pivottable(P.data,[],PD.data,'median');
            D.yeo17=pivottable(P.data,[],PY.data,'mode');
            D.numVert=pivottable(P.data,[],PD.data,'length');
            D.good = (D.numVert<800 | D.numVert>3000)+1;
            D.hem=ones(158,1)*h;
            DD=addstruct(DD,D);
        end;
        [Y,DD.newIndx]=sortrows([DD.hem DD.good DD.lobe DD.desikan]);
        DD.regIndx = [1:316]';
        DD.regIndx2 = [1:148 NaN 149:305 NaN 306:314]'; % This is the index into the reduced structure!
        save(fullfile(regDir,'data','162_reorder.mat'),'-struct','DD');
    case 'Correlation_by_task'
        T = load(fullfile(regDir,'glm4','ts_162_tessellation_hem_all.mat'));
        numsubj = size(T.Yres,3);
        N = size(T.Yres,1);
        X=load(fullfile(regDir,'data','162_reorder.mat'));
        ii = X.regIndx(X.good==1);
        
        type = {'Yhatm','Yhatr','Yres'};
        blockSize = N/16;
        Z=kron(eye(16),ones(blockSize,1));
        
        for t=1:3
            Data = T.(type{t})(:,ii,:);
            for s=1:21
                Data(:,:,s) = Data(:,:,s) - Z*pinv(Z)*Data(:,:,s);  % subtract run means
                
                D=dload(fullfile(behavDir,subj_name{s+1},sprintf('sc1_%s.dat',subj_name{s+1})));
                D=getrow(D,D.runNum>=51);
                [taskName,~,D.task]=unique(D.taskName);
                for i=1:length(D.task)
                    D.run(i,1) = D.runNum(i)-50;
                    D.start(i,1) = D.TRreal(i)+(D.run(i)-1)*N/16+5;
                    D.timeIndx(i,:)=[D.start(i,1):D.start(i,1)+28];
                end;
                D.sess = (D.run>8)+1;
                for se=1:2
                    for task=1:17
                        idx=find(D.sess==se & D.task==task);
                        tidx=vec(D.timeIndx(idx,:));
                        C{t}(:,:,task,s,se)=cov(Data(tidx,:,s));
                    end;
                end;
            end;
        end;
        save(fullfile(regDir,'data','correlation_by_task.mat'),'C','taskName');
    case 'MDS_plot'
        load((fullfile(regDir,'data','correlation_by_task.mat')));
        [numROI,numTask,numSubj,numSess]=size(M);
        for s=1:numSubj
            G(:,:,s)=(Mhat(:,:,s,1)'*Mhat(:,:,s,2)+Mhat(:,:,s,2)'*Mhat(:,:,s,1))/2;
        end;
        Gm=mean(G,3);
        C=indicatorMatrix('allpairs',[1:17]');
        Y=pcm_classicalMDS(Gm,'contrast',C,'rotation','varimax');
        scatterplot3(Y(:,1),Y(:,2),Y(:,3),'label',taskName);
    case 'Between_task_differences'
        load(fullfile(regDir,'data','correlation_by_task.mat'));
        [numROI,~,numTask,numSubj,numSess]=size(C{1});
        T=[];
        titstr={'Yhatm','Yhatr','Yres'};
        for t=1:3
            for s=1:numSubj
                for i=1:17
                    CC{t}(:,i) = rsa_vectorizeRDM(C{t}(:,:,i,s,1))';
                    CC{t}(:,i+17)= rsa_vectorizeRDM(C{t}(:,:,i,s,2))';
                end;
                R{t}(:,:,s)=corr(CC{t});
            end;
        end;
        
        for t=1:3
            subplot(3,2,t);
            imagesc(mean(R{t},3));
            drawline(17.5,'dir','horz');
            drawline(17.5);
            title(titstr{t});
            
            for s=1:numSubj
                D.type=[t;t];
                D.SN  =[s;s];
                A = R{t}(18:end,1:17,s);
                K=size(A,1);
                D.avrgR(1,1) = mean(diag(A));
                D.avrgR(2,1) = (sum(sum(A))-trace(A))/(K*(K-1));
                D.wb   = [1;2];
                T=addstruct(T,D);
            end;
        end;
        
        subplot(3,2,[5 6]);
        barplot([T.type T.wb],T.avrgR);
        set(gcf,'PaperPosition',[2 2 4 7]);
        wysiwyg;
        varargout={T};
    case 'tsVarianceByTask'
        load(fullfile(regDir,'data','correlation_by_task.mat'));
        [numReg,~,numTask,numSubj,numSess]= size(Chat);
        
        % Before combining tasks, see what the variance introduced by each
        % task in the different regions
        T=[];
        for s=1:numSubj
            for t=1:numTask
                D.SN = s;
                D.task = t;
                mChat = mean(Chat(:,:,t,s,:),5);
                D.varHat = mean(diag(mChat)); % variance introduced by task
                mCres = mean(Cres(:,:,t,s,:),5);
                D.varRes = mean(diag(mCres)); % variance introduced by task
                T=addstruct(T,D);
            end;
        end;
        lineplot(T.task,[T.varHat T.varRes]);
    case 'tsWhiteStim'
        type='mean';
        load(fullfile(regDir,'data','correlation_by_task.mat'));
        [numReg,~,numTask,numSubj,numSess]= size(Chat);
        
        for s=1:numSubj
            for t=1:numTask
                switch (type)
                    case 'mean'
                        B(t,:,s)=mean(Mhat(:,t,s,:),4)';
                end;
            end;
        end;
        
        % Cost: variance Inflation factor
        for s=1:numSubj
            fprintf('subject %d\n',s);
            PR.objective = @(x) estimationVariance(x,B(:,:,s));
            PR.solver = 'fminsearch';
            PR.x0    = zeros(numTask,1);
            PR.options.MaxFunEvals = 2000000;
            PR.options.MaxIter = 2000000;
            [x,f]= fminsearch(PR);
            W(s,:)=exp(x)/sum(exp(x))
        end;
        varargout={W};
    otherwise
        disp('there is no such case.')
end;

function C=intersubj_corr(Y)
numSubj=size(Y,3);
for i=1:numSubj
    for j=1:numSubj
        C(i,j)=nansum(nansum(Y(:,:,i).*Y(:,:,j)))/...
            sqrt(nansum(nansum(Y(:,:,i).*Y(:,:,i)))*...
            nansum(nansum(Y(:,:,j).*Y(:,:,j))));
    end;
end;

function [DiffSubjSameSess,SameSubjDiffSess,DiffSubjDiffSess]=bwse_corr(C);
N=size(C,1);
n=N/2;
sess = [ones(1,n) ones(1,n)*2];
subj = [1:n 1:n];
SameSess = bsxfun(@eq,sess',sess);
SameSubj = bsxfun(@eq,subj',subj);
DiffSubjSameSess=mean(C(~SameSubj & SameSess));
SameSubjDiffSess=mean(C(SameSubj & ~SameSess));
DiffSubjDiffSess=mean(C(~SameSubj & ~SameSess));
ROI
%  Cost function: total estimation variance of stimulation
% after removing the mean activity
function tvar=estimationVariance(x,B)
w=exp(x);
P=size(B,2);
wB=bsxfun(@times,B,w);
m=sum(wB,1)/sum(w);
X=bsxfun(@minus,B,m);
Gr=X'*diag(w)*X/sum(w);
L=eye(P)*0.01;
tvar=trace(inv(Gr+L));

function dircheck(dir)
if ~exist(dir,'dir');
    fprintf('%s doesn''t exist. Creating one now. \n',dir);
    mkdir(dir);
end

