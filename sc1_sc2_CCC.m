function varargout=sc1_sc2_CCC(what,varargin)
% Matlab script to do connectivity analyses for sc1 and sc2 experiments
% Does time series extraction
%==========================================================================
% % (1) Directories
baseDir          = '/Users/maedbhking/Documents/Cerebellum_Cognition';
% baseDir           = '/Volumes/MotorControl/data/super_cerebellum_new/';

atlasDir='/Users/maedbhking/Documents/Atlas_templates/';

studyDir{1}     =fullfile(baseDir,'sc1');
studyDir{2}     =fullfile(baseDir,'sc2');
studyStr        ={'SC1','SC2','SC12'};
behavDir        =['data'];
imagingDir      =['imaging_data'];
suitDir         =['suit'];
caretDir        =['surfaceCaret'];
regDir          =['RegionOfInterest/'];
connDir         =['connectivity'];
encodeDir       =['encoding'];

%==========================================================================

% % (2) Hemisphere and Region Names
funcRunNum = [51,66];  % first and last behavioural run numbers (16 runs per subject)

run = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'};

subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11',...
    's12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24',...
    's25','s26','s27','s28','s29','s30','s31'};

returnSubjs=[2,3,4,6,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31]; % also known as 'good_subjs'

hem       = {'lh','rh'};
hemName   = {'LeftHem','RightHem'};
exper     = {'sc1','sc2'};
glm       = 4;

%==========================================================================

switch(what)
        
    case 'SIGNAL:meanTS'
        % sc1_connectivity('TS_get_meants',[2 3 6 8 9 10 12 17:22],'sc2',4,'162_tessellation_hem');
        % sc1_connectivity('TS_get_meants',[2:22],'sc1',4,'162_tessellation_hem');
        sn=varargin{1}; % subjNum
        exper=varargin{2}; % Experiment string 'sc1' or 'sc2'
        
        B = [];
        glm=4;
        type='Cerebellum_grey';
        glmDir =fullfile(baseDir,exper,sprintf('/GLM_firstlevel_%d',glm));
        subjs=length(sn);
        
        for s=1:subjs,
            glmDirSubj=fullfile(glmDir, subj_name{sn(s)});
            load(fullfile(glmDirSubj,'SPM_light.mat'));
            %            load(fullfile(glmDirSubj,'SPM.mat'));
            T=load(fullfile(glmDirSubj,'SPM_info.mat'));
            
            % load data
            tic;
            load(fullfile(sc1Dir,regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',type))); % 'regions' are defined in 'ROI_define'
            SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'sc1',imagingDir,subj_name{sn(s)}));
            
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
            filename=(fullfile(baseDir,exper,regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('ts_%s.mat',type)));
            save(filename,'Data','Yres','Yhatm','Yhatr','B');
            
            fprintf('ts saved for %s (%s) for %s \n',subj_name{sn(s)},sprintf('glm%d',glm),type);
        end
    case 'SIGNAL:voxelTS'
        % sc1_connectivity('TS_get_ts',[2:22],'sc2',4,'Cerebellum_grey');
        % sc1_connectivity('TS_get_ts',[2 3 6 8 9 10 12 17:22],'sc2',4,'Cerebellum_grey');
        sn=varargin{1}; % subjNum
        exper=varargin{2}; % Experiment string 'sc1' or 'sc2
        
        B = [];
        glm=4;
        type='Cerebellum_grey';
        glmDir =fullfile(baseDir,exper,sprintf('GLM_firstlevel_%d',glm));
        subjs=length(sn);
        
        for s=1:subjs,
            glmDirSubj=fullfile(glmDir, subj_name{sn(s)});
            load(fullfile(glmDirSubj,'SPM_light.mat'));
            %             load(fullfile(glmDirSubj,'SPM.mat'));
            T=load(fullfile(glmDirSubj,'SPM_info.mat'));
            
            % load data
            tic;
            load(fullfile(sc1Dir,regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',type))); % 'regions' are defined in 'ROI_define'
            SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'sc1',imagingDir,subj_name{sn(s)}));
            
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
            filename=(fullfile(baseDir,exper,regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('ts_%s.mat',type)));
            save(filename,'Yres','Yhatm','Yhatr','B');
            fprintf('ts saved for %s (%s) for %s \n',subj_name{sn(s)},sprintf('glm%d',glm),type);
        end
    case 'SIGNAL:meanBeta'
        % In the new version, this always stores the betas like the occurr
        % in the GLM
        type    = '162_tessellation_hem';  % 'Cerebellum_grey';
        exper   = {'sc1','sc2'};
        sn      = goodsubj;
        glm     = 4;
        vararginoptions(varargin,{'sn','glm','type'});
        
        % Loop over experiments
        for e=1:2
            myRegDir = fullfile(baseDir,exper{e},regDir);
            for s=1:length(sn);
                fprintf('%d\n',sn(s));
                filename = (fullfile(myRegDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('ts_%s.mat',type)));
                D          = load(filename);
                glmDirSubj = fullfile(baseDir,exper{e},sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)});
                T          = load(fullfile(glmDirSubj,'SPM_info.mat'));
                
                if (strcmp(type,'162_tessellation_hem'))
                    Tcort=load(fullfile(baseDir,'sc1',regDir,'data','162_reorder.mat'));
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
        
    case 'STRUCTURE:TS'      % Make a structure of all cortical time series of all subject
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
    case 'STRUCTURE:meanBeta' %  Gets the mean betas for each session and experiment
        model=varargin{1}; % '162_tessellation_hem'
        data=varargin{2}; % 'Cerebellum_grey';
        incInstr=varargin{3}; % 0 or 1
        
        %         vararginoptions(varargin,{'exper','glm','yname','xname','incInstr'});
        
        % Load the betas for all the tasks
        for e=1:2,
            myRegDir = fullfile(baseDir,exper{e},regDir);
            YD{e}=load(fullfile(myRegDir,sprintf('glm%d',glm),sprintf('mbeta_%s_all.mat',data)));
            XD{e}=load(fullfile(myRegDir,sprintf('glm%d',glm),sprintf('mbeta_%s_all.mat',model)));
        end;
        
        % Make an integrated structure, either including instruction or
        % rest
        if (incInstr)
            T=dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        else
            T=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        end;
        
        % Make an integrated Structure
        T1=T;
        T2=T;
        T1.sess=ones(length(T.condNum),1)*1;
        T2.sess=ones(length(T.condNum),1)*2;
        T=addstruct(T1,T2);
        for s=1:size(YD{1}.B,1)
            X{s}=[];
            Y{s}=[];
            for se=1:2
                if (incInstr==0)  % If no instruction, add rest and center within session
                    for e=1:2
                        YD{e}.B{s,se}=[YD{e}.B{s,se}(2:end,:);zeros(1,size(YD{e}.B{s,se},2))];
                        XD{e}.B{s,se}=[XD{e}.B{s,se}(2:end,:);zeros(1,size(XD{e}.B{s,se},2))];
                        YD{e}.B{s,se}=bsxfun(@minus,YD{e}.B{s,se},mean(YD{e}.B{s,se}));
                        XD{e}.B{s,se}=bsxfun(@minus,XD{e}.B{s,se},mean(XD{e}.B{s,se}));
                    end;
                end;
                Y{s}=[Y{s};YD{1}.B{s,se};YD{2}.B{s,se}];
                X{s}=[X{s};XD{1}.B{s,se};XD{2}.B{s,se}];
            end;
        end;
        varargout = {X,Y,T};
        
    case 'MODEL:construct'
        type=varargin{1}; % 'cortical_lobes','yeo','desikan','cerebellum','tasks','feature','tasks_saccades_hands','desikan_hem','yeo_hem','yeo_reduced'
        study=varargin{2}; % which study ? 1 or 2
        
        sn=returnSubjs;
        subjs=length(returnSubjs);
        glm=4;
        
        for s=1:subjs,
            X=[];
            connSubjDir = fullfile(studyDir{study},connDir,sprintf('glm%d',glm),subj_name{sn(s)}); dircheck(connSubjDir);
            glmSubjDir =[studyDir{study} sprintf('/GLM_firstlevel_%d/',glm) subj_name{sn(s)}];
            D=load(fullfile(glmSubjDir,'SPM_info.mat'));
            
            % load betas for different cortical surfaces
            load(fullfile(studyDir{study},regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('betas_%s.mat',type)));
            
            switch type,
                case 'cortical_lobes'
                    for r=1:length(B),
                        X.Xx(r,:)=mean(B{r}.betasUW,2);
                        X.idx(r,:)=r;
                    end
                    X.Xx(:,length(D.SN)+1:length(D.SN)+numel(run))=0; % zero out last 'intercept' run
                case 'yeo'
                    networks=[2:18];
                    for r=1:17, % we don't want 'FreeSurfer_Defined_Medial_Wall'
                        X.Xx(r,:)=mean(B{networks(r)}.betasUW,2);
                        X.idx(r,1)=networks(r);
                    end
                    X.Xx(:,length(D.SN)+1:length(D.SN)+numel(run))=0; % zero out last 'intercept' run
                case 'yeo_hem'
                    networks=[2:18,20:36]; % not interested in freesurfer defined medial wall
                    paintFileNums=repmat([2:18],1,2);
                    for r=1:34,
                        X.Xx(r,:)=mean(B{networks(r)}.betasUW,2);
                        X.idx(r,1)=paintFileNums(r);
                    end
                    X.Xx(:,length(D.SN)+1:length(D.SN)+numel(run))=0; % zero out last 'intercept' run
                case 'desikan'
                    networks=[2:36]; % remove 'medial wall'
                    for r=1:35,
                        X.Xx(r,:)=mean(B{networks(r)}.betasUW,2);
                        X.idx(r,1)=networks(r);
                    end
                    X.Xx(:,length(D.SN)+1:length(D.SN)+numel(run))=0; % zero out last 'intercept' run
                case 'desikan_hem'
                    networks=[2:36,38:72];
                    paintFileNums=repmat([2:36],1,2);
                    for r=1:70,
                        X.Xx(r,:)=mean(B{networks(r)}.betasUW,2);
                        X.idx(r,1)=paintFileNums(r);
                    end
                    X.Xx(:,length(D.SN)+1:length(D.SN)+numel(run))=0; % zero out last 'intercept' run
                case '162_tessellation'
                    tessels=[1:148,150:158];
                    for r=1:157, % we don't want 'FreeSurfer_Defined_Medial_Wall' - tessel number 149
                        X.Xx(r,:)=mean(B{tessels(r)}.betasUW,2);
                        X.idx(r,1)=tessels(r);
                    end
                    X.Xx(:,length(D.SN)+1:length(D.SN)+numel(run))=0; % zero out last 'intercept' run
                case '162_tessellation_hem'
                    I=load(fullfile(regDir,'data','162_reorder.mat'));
                    for r=1:length(B),
                        % determine good and bad tessels
                        if I.good(r)==1,
                            X.Xx(r,:)=mean(B{r}.betasUW,2);
                        end
                        X.idx(r,1)=r;
                        goodIdx(r,1)=I.good(r);
                    end
                    X.Xx(:,size(X.Xx,2)-numel(run)+1:size(X.Xx,2))=0; % zero out last 'intercept' run
                    X=getrow(X,goodIdx==1);
                case 'tasks'
                    load(fullfile(encodeDir,subj_name{sn(s)},sprintf('Y_info_glm%d_grey.mat',glm)));
                    X.idx=[1:29]';
                    X.Xx=indicatorMatrix('identity_p',D.cond)';
                    X.Xx(:,length(D.SN)+1:length(D.SN)+numel(run))=0; % zero out last 'intercept' run
            end
            outFile=sprintf('%s_glm%d_model.mat',type,glm);
            save(fullfile(connSubjDir,outFile),'X');
            fprintf('%s-weighted betas (glm%d) (X) have been computed for %s \n',type,glm,subj_name{sn(s)});
        end
    case 'MODEL:covariance'
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
        
    case 'HYPERPARAMS:run'
        step=varargin{1}; % 'weights' or 'eval' ?
        method=varargin{2}; % 'l1' 'l2' 'l1l2'
        
        lambdas=sc1_sc2_CCC('getLambda',method);
        
        switch step,
            case 'weights'
                %                 for i=1:2,
                for l=1:length(lambdas),
                    sc1_sc2_CCC('HYPERPARAMS:weights',returnSubjs,'ridgeFixed',2,'l2',lambdas(l))
                end
                %                 end
            case 'eval'
                for l=1:length(lambdas),
                    sc1_sc2_CCC('HYPERPARAMS:eval',returnSubjs,'ridgeFixed',2,'l2',lambdas(l))
                    sc1_sc2_CCC('HYPERPARAMS:eval',returnSubjs,'ridgeFixed',1,'l2',lambdas(l))
                end
        end
    case 'HYPERPARAMS:weights' % Determine optimal lambdas for subsequent modelling
        sn=varargin{1}; % subjNum
        method=varargin{2}; % 'ridgeFixed', other options: nonNegative etc
        study=varargin{3}; % building on which study ? 1 or 2 [1]
        
        % vararginoptions
        meanSub=1;      % Mean Pattern Subtract?
        incInstr=0;     % Include instruction ?
        l1=0;
        l2=0;
        signal='fitted';
        model='162_tessellation_hem';
        data='Cerebellum_grey';
        trainMode='crossed';
        
        vararginoptions(varargin(4:end),{'l1','l2','signal','model','data','trainMode'});
        
        % how many samples for bootstrap
        numSamples=1;
        
        T=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % which study are we taking ?
        subset = T.StudyNum==study; % Indices of the experiments / conditions that we want to use
        subset = [subset;subset]; % Make the subset out onto both session
        
        % get X and Y data
        [X,Y,S]=sc1_sc2_CCC('STRUCTURE:meanBeta',model,data,incInstr);
        
        % 'weights' fileName
        fileName=sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d.mat',model,signal,method,study,round(l1),round(l2));
        
        % loop over # of subjects
        for s=1:length(sn),
            RR=[];
            
            hyperParamDir=fullfile(studyDir{2},connDir,'glm4','weights_hyperParams',subj_name{sn(s)});dircheck(hyperParamDir);
            
            if(meanSub)
                for sess=1:2
                    indx=find(subset & S.sess==sess);
                    X{s}(indx,:)=bsxfun(@minus,X{s}(indx,:),mean(X{s}(indx,:)));
                    Y{s}(indx,:)=bsxfun(@minus,Y{s}(indx,:),mean(Y{s}(indx,:)));
                end;
            end;
            
            % crossed or uncrossed across sessions ?
            trainXindx=[find(subset & S.sess==1);find(subset & S.sess==2)];
            switch (trainMode)
                case 'crossed'
                    trainYindx=[find(subset & S.sess==2);find(subset & S.sess==1)];
                case 'uncrossed'
                    trainYindx=[find(subset & S.sess==1);find(subset & S.sess==2)];
            end;
            
            xx = X{s}(trainXindx,:);
            yy = Y{s}(trainYindx,:);
            
            % loop over # of bootstraps
            for b=1:numSamples,
                
                % resample the y data
                %                 yy_boot=datasample(yy,size(yy,1));
                
                R.SN=sn(s);
                W=sc_connect_fit(yy,xx,method,'lambda',[l1 l2]);
                R.W={W};
                R.lambda=[l1 l2];
                R.model={model};
                R.signal={signal};
                R.trainMode={trainMode};
                R.study=study;
                R.method={method};
                %                 R.boot=b;
                RR=addstruct(RR,R);
                clear W yy_boot
                fprintf('subj%d-bootstrap%d-lambda%2.2d %2.2d \n',sn(s),b,[l1 l2]);
            end
            % save connectivity weights
            save(fullfile(hyperParamDir,fileName),'-struct','RR')
            clear RR R
        end
    case 'HYPERPARAMS:eval'    % Evaluates predictive performance of connectivity model
        %  M: is the model structure with the connectivity weights in it (M.W)
        % 'subset': what data should the evaluation be based upon?
        % 'splitby': by which variable should the evaluation be split?
        % 'meanSub': Mean pattern subtraction before evaluation?
        sn=varargin{1}; % subjNum
        method=varargin{2}; % 'ridgeFixed' or 'nonNegative'
        study=varargin{3}; % evaluating on which study, 1 or 2 ? [1,2]
        % example:
        % sc1_sc2_CCC('CONNECTIVITY:evaluate','162_tessellation_hem','Cerebellum_grey','fitted','ridgeFixed')
        
        % if the connectivity weights are built on sc1, then evaluate on
        % sc2
        if study==1,
            studyModel=2; % connectivity weights: which study ?
        elseif  study==2,
            studyModel=1;
        end
        
        l1=0;
        l2=0;
        signal='fitted';
        model='162_tessellation_hem';
        data='Cerebellum_grey';
        trainMode='crossed';
        condType='unique';
        incInstr=0;
        meanSub = 0; % Mean pattern subtraction before prediction?
        
        T = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        vararginoptions(varargin(4:end),{'l1','l2','subset','meanSub','incInstr','data','trainMode','model','signal','condType'}); % 'condType' is either 'unique or 'all'
        
        % Get all the mean betas and prepare the evaluation data (Y)
        [X,Y,S]=sc1_sc2_CCC('STRUCTURE:meanBeta',model,data,incInstr);
        
        % what are we evaluating ?
        subset = T.StudyNum==study;
        switch condType,
            case 'unique'
                splitby = (T.condNum & T.overlap==0).*T.condNum;
            case 'all'
                splitby = T.condNum;
            case 'shared'
                splitby = (T.condNum & T.overlap==1).*T.condNum;
        end
        S.subset= [subset;subset];
        if (isempty(splitby))
            splitby = ones(length(T.StudyNum),1);
        end;
        S.splitby=[splitby;splitby];  % Double for the two conditions
        sS = getrow(S,S.subset);
        splits = unique(sS.splitby); % related to one task set
        splits(splits==0)=[]; % remove zero
        
        % loop over subjects
        for s=1:length(sn) % subjects
            
            RR=[];
            
            hyperParamWeightDir=fullfile(studyDir{2},connDir,'glm4','weights_hyperParams',subj_name{sn(s)});dircheck(hyperParamWeightDir);
            hyperParamEvalDir=fullfile(studyDir{2},connDir,'glm4','eval_hyperParams',subj_name{sn(s)});dircheck(hyperParamEvalDir);
            
            % load in the connectivity weights (X): option to load in multiple
            % models (each can have different methods, trainings etc)
            % M=sc1_sc2_CCC('CONNECTIVITY:runModel',returnSubjs,model,data,signal,method,'trainMode',trainMode,studyModel,'lambdaL1',lambdaL1,'lambdaL2',lambdaL2);
            M=load(fullfile(hyperParamWeightDir,sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d.mat',model,signal,method,studyModel,round(l1),round(l2))));
            
            % calculate predictions for each bootstrap
            numSamples=length(unique(M.boot));
            for b=1:1:numSamples,
                
                Mm=getrow(M,M.boot==b & strcmp(M.trainMode,trainMode));
                
                indx = sum(abs(Y{s}(:,:)))>0 & ~isnan(sum(Y{s}(:,:))) & ~isnan(sum(Mm.W{1})); % why is the absolute value being taken ?
                %                     indx = ~isnan(sum(Y{s}(:,:))) & ~isnan(sum(Mm.W{s}));
                for sp=1:length(splits); % task conditions
                    if(meanSub)
                        for sess=1:2 % sessions
                            indx=find(S.subset & S.sess==sess & S.splitby==splits(sp));
                            X{s}(indx,:)=bsxfun(@minus,X{s}(indx,:),mean(X{s}(indx,:)));
                            Y{s}(indx,:)=bsxfun(@minus,Y{s}(indx,:),mean(Y{s}(indx,:)));
                        end;
                    end;
                    testAindx=[find(S.subset & S.sess==1 & S.splitby==splits(sp));...
                        find(S.subset & S.sess==2 & S.splitby==splits(sp))];
                    testBindx=[find(S.subset & S.sess==2 & S.splitby==splits(sp));...
                        find(S.subset & S.sess==1 & S.splitby==splits(sp))];
                    predY   = X{s}(testAindx,:)*Mm.W{1};                    % Predicted Y using crossvalidation
                    predYnc = X{s}(testBindx,:)*Mm.W{1};                    % Predicted Y not crossvalidated
                    SSP   = sum(sum(predY(:,indx).^2));             % Sum of square of predictions
                    SSY   = sum(sum(Y{s}(testBindx,indx).^2));               % Sum of squares of data
                    SSCp  = sum(sum(predY(:,indx).*predYnc(:,indx))); % Covariance of Predictions
                    SSCy  = sum(sum(Y{s}(testAindx,indx).*Y{s}(testBindx,indx)));  % Covariance of Y's
                    SSCn  = sum(sum(predYnc(:,indx).*Y{s}(testBindx,indx)));   % Covariance of non-cross prediction and data
                    SSCc  = sum(sum(predY(:,indx).*Y{s}(testBindx,indx)));   % Covariance of cross prediction and data
                    R.SN    = Mm.SN;
                    R.lambda = Mm.lambda(:,:);
                    R.method = Mm.method;
                    R.trainMode = Mm.trainMode;
                    R.model = Mm.model;
                    R.Rcv   = SSCc ./ sqrt(SSY.*SSP); % Double-Crossvalidated predictive correlation
                    R.Rnc   = SSCn ./ sqrt(SSY.*SSP); % Not double-crossvalidated predictive correlation
                    R.Ry    = SSCy ./ SSY;            % Reliability of data: noise ceiling
                    R.Rp    = SSCp ./ SSP;            % Reliability of prediction
                    R.splitby = splits(sp);
                    R.split = {condType};
                    R.boot  = b;
                    RR = addstruct(RR,R);
                end;
                fprintf('pred done for %s-boot%d \n',subj_name{sn(s)},b)
            end;
            save(fullfile(hyperParamEvalDir,sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d_%s.mat',model,signal,method,study,round(l1),round(l2),condType)),'-struct','RR')
            clear RR Mm
        end
    case 'HYPERPARAMS:plot'    % Plots reg-perf curve for each subject
        sn=varargin{1};
        model=varargin{2};
        signal=varargin{3};
        method=varargin{4};
        study=varargin{5};
        l1=varargin{6};
        l2=varargin{7};
        condType=varargin{8};
        
        vararginoptions({varargin{9:end}},{'CAT'}); % option if doing individual map analysis
        
        %
        l2=round(l2);
        l1=round(l1);
        
        TT=[];
        for s=1:length(sn),
            for i=1:length(l1),
                for ii=1:length(l2),
                    T=load(fullfile(studyDir{2},connDir,'glm4','eval_hyperParams',subj_name{sn(s)},...
                        sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d_%s.mat',model,signal,method,study,l1(i),l2(ii),condType)));
                    TT=addstruct(TT,T);
                end
            end
        end
        
        if exist('CAT'),
            xyplot(TT.lambda(:,2),TT.Rcv,TT.lambda(:,2),'split',TT.lambda(:,2),'CAT',CAT);
        else
            xyplot(TT.lambda(:,2),TT.Rcv,TT.lambda(:,2),'split',TT.lambda(:,2),'style_thickline');
        end
        
    case 'CONNECT:run'
        step=varargin{1}; % 'weights' or 'eval' ?
        method=varargin{2}; % 'l1' 'ridgeFixed' 'l1l2'
        
        lambdas=sc1_sc2_CCC('getLambda',method);
        
        switch step,
            case 'weights'
                for i=1:2,
                    for l=1:length(lambdas),
                        sc1_sc2_CCC('CONNECT:weights',returnSubjs,method,i,'l2',lambdas(l))
                    end
                end
            case 'eval'
                for i=1:2,
                    for l=1:length(lambdas),
                        sc1_sc2_CCC('CONNECT:eval',returnSubjs,method,i,'l2',lambdas(l))
                    end
                end
        end
    case 'getLambda'
        method=varargin{1}; % 'L1','L2','Lasso'
        
        switch method,
            case 'L1'
            case 'L2'
                n=20;
                lambdas=exp(linspace(log(1000),log(100000),n)); % 20 log spaced values from 10 to 1000 (or 1000 to 10000)
            case 'Lasso'
        end
        
        varargout={lambdas};
    case 'CONNECT:weights'
        sn=varargin{1}; % subjNum
        method=varargin{2}; % 'L2', other options: nonNegative etc
        study=varargin{3}; % building on which study ? 1 or 2 [1]
        
        T=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % vararginoptions
        meanSub=1;      % Mean Pattern Subtract?
        incInstr=0;     % Include instruction ?
        l1=0;
        l2=0;
        signal='fitted';
        model='162_tessellation_hem';
        data='Cerebellum_grey';
        trainMode={'crossed'};
        
        vararginoptions(varargin(4:end),{'l1','l2','signal','model','data','trainMode'});
        
        % 'weights' fileName
        modelName=fullfile(studyDir{2},connDir,'glm4','weights',method,...
            sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d.mat',model,signal,method,study,round(l1),round(l2)));
        
        subset = T.StudyNum==study; % Indices of the experiments / conditions that we want to use
        subset = [subset;subset]; % Make the subset out onto both session
        
        switch signal,
            case 'fitted'
                [X,Y,S]=sc1_sc2_CCC('STRUCTURE:meanBeta',model,data,incInstr);
            case 'residual'
                [X,Y,S]=sc1_sc2_CCC('');
            case 'totalTS'
                [X,Y,S]=sc1_sc2_CCC('');
        end
        
        trainXindx=[find(subset & S.sess==1);find(subset & S.sess==2)];
        
        RR=[];
        for s=1:length(sn)
            if(meanSub)
                for sess=1:2
                    indx=find(subset & S.sess==sess);
                    X{s}(indx,:)=bsxfun(@minus,X{s}(indx,:),mean(X{s}(indx,:)));
                    Y{s}(indx,:)=bsxfun(@minus,Y{s}(indx,:),mean(Y{s}(indx,:)));
                end;
            end;
            
            for t=1:length(trainMode),
                
                switch (trainMode{t})
                    case 'crossed'
                        trainYindx=[find(subset & S.sess==2);find(subset & S.sess==1)];
                    case 'uncrossed'
                        trainYindx=[find(subset & S.sess==1);find(subset & S.sess==2)];
                end;
                
                xx = X{s}(trainXindx,:);
                yy = Y{s}(trainYindx,:);
                
                R.SN = sn(s);
                [W,R.fR2m,R.fRm] = sc_connect_fit(yy,xx,method,'lambda',[l1 l2]);
                R.W={W};
                R.lambda = [l1 l2];
                R.model={model};
                R.signal={signal};
                R.trainMode={trainMode{t}};
                R.trainIdx=t;
                R.study=study;
                R.method  = {method};
                RR=addstruct(RR,R);
                clear W
            end
            fprintf('%d\n',sn(s));
        end
        % save connectivity weights
        save(modelName,'-struct','RR','-v7.3')
    case 'CONNECT:eval' % Evaluates predictive performance of connectivity model
        %  M: is the model structure with the connectivity weights in it (M.W)
        % 'subset': what data should the evaluation be based upon?
        % 'splitby': by which variable should the evaluation be split?
        % 'meanSub': Mean pattern subtraction before evaluation?
        sn=varargin{1}; % subjNum
        method=varargin{2}; % 'L2' or 'nonNegative'
        study=varargin{3}; % evaluating on which study, 1 or 2 ? [1,2]
        % example:
        % sc1_sc2_CCC('CONNECTIVITY:evaluate','162_tessellation_hem','Cerebellum_grey','fitted','ridgeFixed')
        
        % if the connectivity weights are built on sc1, then evaluate on
        % sc2
        if study==1,
            studyModel=2; % connectivity weights: which study ?
        elseif  study==2,
            studyModel=1;
        end
        
        l1=0;
        l2=0;
        signal='fitted';
        model='162_tessellation_hem';
        data='Cerebellum_grey';
        trainMode={'crossed'};
        condType='unique';
        incInstr=0;
        meanSub = 0; % Mean pattern subtraction before prediction?
        voxels=0;
        
        T = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        vararginoptions(varargin(4:end),{'l1','l2','subset','meanSub','incInstr','data','trainMode','model','signal','condType','voxels'}); % 'condType' is either 'unique or 'all'
        
        % Get all the mean betas and prepare the evaluation data (Y)
        [X,Y,S]=sc1_sc2_CCC('STRUCTURE:meanBeta',model,data,incInstr);
        
        % what are we evaluating ?
        subset = T.StudyNum==study;
        switch condType,
            case 'unique'
                splitby = (T.condNum & T.overlap==0).*T.condNum;
            case 'all'
                splitby = T.condNum;
            case 'shared'
                splitby = (T.condNum & T.overlap==1).*T.condNum;
        end
        S.subset= [subset;subset];
        if (isempty(splitby))
            splitby = ones(length(T.StudyNum),1);
        end;
        S.splitby=[splitby;splitby];  % Double for the two conditions
        sS = getrow(S,S.subset);
        splits = unique(sS.splitby); % related to one task set
        splits(splits==0)=[]; % remove zero
        
        weightDir=fullfile(studyDir{2},connDir,'glm4','weights',method);dircheck(weightDir);
        evalDir=fullfile(studyDir{2},connDir,'glm4','eval');dircheck(evalDir);
        
        % load in the connectivity weights (X): option to load in multiple
        % models (each can have different methods, trainings etc)
        % M=sc1_sc2_CCC('CONNECTIVITY:runModel',returnSubjs,model,data,signal,method,'trainMode',trainMode,studyModel,'lambdaL1',lambdaL1,'lambdaL2',lambdaL2);
        M=load(fullfile(weightDir,sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d.mat',model,signal,method,studyModel,round(l1),round(l2))));
        
        % eval file name - output and check to see if connectivity eval exists [load if it does exist]
        if exist('voxels'),
            evalFile=sprintf('%s_%s_%s_sc%d_%s_voxels.mat',model,signal,method,study,condType);
            RR=[];
        else
            evalFile=sprintf('%s_%s_%s_sc%d_%s.mat',model,signal,method,study,condType);
            RR=load(fullfile(evalDir,evalFile));
        end

        for s=1:length(sn) % subjects
            
            Mm=getrow(M,M.SN==sn(s) & strcmp(M.trainMode,trainMode));
            numVox=size(Mm.W{1},2);
            
            %             indx = sum(abs(Y{s}(:,:)))>0 & ~isnan(sum(Y{s}(:,:))) &
            %             ~isnan(sum(Mm.W{1})); % introduces diff number of voxels --
            %             introduces complications later on when reslicing into suit
            %             space
            indx=ones(1,numVox);
            
            for sp=1:length(splits); % task conditions
                if(meanSub)
                    for sess=1:2, % sessions
                        indx=find(S.subset & S.sess==sess & S.splitby==splits(sp));
                        X{s}(indx,:)=bsxfun(@minus,X{s}(indx,:),mean(X{s}(indx,:)));
                        Y{s}(indx,:)=bsxfun(@minus,Y{s}(indx,:),mean(Y{s}(indx,:)));
                    end;
                end;
                testAindx=[find(S.subset & S.sess==1 & S.splitby==splits(sp));...
                    find(S.subset & S.sess==2 & S.splitby==splits(sp))];
                testBindx=[find(S.subset & S.sess==2 & S.splitby==splits(sp));...
                    find(S.subset & S.sess==1 & S.splitby==splits(sp))];
                predY   = X{s}(testAindx,:)*Mm.W{1};                    % Predicted Y using crossvalidation
                predYnc = X{s}(testBindx,:)*Mm.W{1};                    % Predicted Y not crossvalidated
                
                switch voxels
                    case 1
                        SSP   = sum(predY(:,indx).^2);             % Sum of square of predictions
                        SSY   = sum(Y{s}(testBindx,indx).^2);               % Sum of squares of data
                        SSCp  = sum(predY(:,indx).*predYnc(:,indx)); % Covariance of Predictions
                        SSCy  = sum(Y{s}(testAindx,indx).*Y{s}(testBindx,indx));  % Covariance of Y's
                        SSCn  = sum(predYnc(:,indx).*Y{s}(testBindx,indx));   % Covariance of non-cross prediction and data
                        SSCc  = sum(predY(:,indx).*Y{s}(testBindx,indx));   % Covariance of cross prediction and data
                        R.Rcv   = {SSCc ./ sqrt(SSY.*SSP)}; % Double-Crossvalidated predictive correlation
                        R.Rnc   = {SSCn ./ sqrt(SSY.*SSP)}; % Not double-crossvalidated predictive correlation
                        R.Ry    = {SSCy ./ SSY};            % Reliability of data: noise ceiling
                        R.Rp    = {SSCp ./ SSP};            % Reliability of prediction
                    case 0
                        SSP   = sum(sum(predY(:,indx).^2));             % Sum of square of predictions
                        SSY   = sum(sum(Y{s}(testBindx,indx).^2));               % Sum of squares of data
                        SSCp  = sum(sum(predY(:,indx).*predYnc(:,indx))); % Covariance of Predictions
                        SSCy  = sum(sum(Y{s}(testAindx,indx).*Y{s}(testBindx,indx)));  % Covariance of Y's
                        SSCn  = sum(sum(predYnc(:,indx).*Y{s}(testBindx,indx)));   % Covariance of non-cross prediction and data
                        SSCc  = sum(sum(predY(:,indx).*Y{s}(testBindx,indx)));   % Covariance of cross prediction and data
                        R.Rcv   = SSCc ./ sqrt(SSY.*SSP); % Double-Crossvalidated predictive correlation
                        R.Rnc   = SSCn ./ sqrt(SSY.*SSP); % Not double-crossvalidated predictive correlation
                        R.Ry    = SSCy ./ SSY;            % Reliability of data: noise ceiling
                        R.Rp    = SSCp ./ SSP;            % Reliability of prediction
                end
                R.SN    = Mm.SN;
                R.lambda = Mm.lambda(:,:);
                R.method = Mm.method;
                R.trainMode = Mm.trainMode;
                R.model = Mm.model;
                R.splitby = splits(sp);
                R.split = {condType};
                RR = addstruct(RR,R);
            end;
            fprintf('pred done for %s \n',subj_name{sn(s)})
            clear R
        end
        save(fullfile(evalDir,evalFile),'-struct','RR')
    case 'CONNECT:plot' % loads in multiple methods and compares them
        method=varargin{1}; % {'L2'}, {'nonNegative','L2} etc
        study=varargin{2}; % evaluating on which study, 1 or 2 ? [1,2]
        
        signal='fitted';
        model='162_tessellation_hem';
        trainMode={'crossed'};
        condType={'unique'};
        whatToPlot='bestMethod';
        
        RR=[];
        S=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % Aesthetics
        CAT.markertype='none';
        CAT.errorwidth=1;
        CAT.linecolor={'r','k','b','g'};
        CAT.errorcolor={'r','k','b','g'};
        CAT.linewidth={3, 3, 3, 3};
        CAT.linestyle={'-','-','-','-'};
        
        for m=1:length(method),
            % figure out L1,L2, or lasso
            if strcmp(method{m},'L2'),
                M1=2;
            elseif strcmp(method{m},'L1'),
                M1=1;
            elseif strcmp(method{m},'winnerTakeAll'),
                M1=2;
            end
            for c=1:length(condType),
                T=load(fullfile(studyDir{2},connDir,'glm4','eval',sprintf('%s_%s_%s_sc%d_%s.mat',model,signal,method{m},study,condType{c})));
                T.methodNum=repmat(m,size(T.SN,1),1);
                T.condNum=repmat(c,size(T.SN,1),1);
                RR=addstruct(RR,T);
            end
        end
        switch whatToPlot,
            case 'condType'
                % 'unique' versus 'all'
                figure(1)
                lineplot(RR.lambda(:,M1),RR.Rcv,'split',RR.split,'leg','auto','CAT',CAT,...
                    'subset',strcmp(RR.method,method{1}) & strcmp(RR.trainMode,'crossed'))
                ylabel('R')
                xlabel('lambda')
                title(sprintf('%s-%s',signal,method{1}))
                set(gca,'YLim',[.1 0.4],'FontSize',14,'ytick',[.1,.15,.2,.25,.3, .35,.4]);
            case 'crossUncross'
                % 'crossed' versus 'uncrossed'
                figure(1)
                lineplot(RR.lambda(:,M1),RR.Rcv,'split',RR.trainMode,'leg','auto','CAT',CAT,...
                    'subset',strcmp(RR.method,method{1}) & strcmp(RR.split,'unique'))
                ylabel('R')
                xlabel('lambda')
                title(sprintf('%s-%s',signal,method{1}))
                set(gca,'YLim',[.1 0.3],'FontSize',14,'ytick',[.1,.15,.2,.25,.3]);
            case 'taskConds'
                
                % split by taskConds (unique)
                condNames=S.condNames(S.StudyNum==study & S.overlap==0);
                figure(1)
                barplot(RR.splitby,RR.Rcv,'subset',...
                    strcmp(RR.method,method{1}) & strcmp(RR.split,'unique') & strcmp(RR.trainMode,'crossed'))
                ylabel('R')
                title(sprintf('%s-%s-unique',signal,method{1}))
                set(gca,'YLim',[0 0.4],'FontSize',8,'XTickLabel',condNames);
                
                % split by taskConds (all)
                condNames=S.condNames(S.StudyNum==study);
                figure(2)
                barplot(RR.splitby,RR.Rcv,'subset',strcmp(RR.method,method{1}) & strcmp(RR.split,'all') & strcmp(RR.trainMode,'crossed') & strcmp(RR.trainMode,'crossed'))
                ylabel('R')
                title(sprintf('%s-%s-all',signal,method{1}))
                set(gca,'YLim',[0 0.55],'FontSize',8,'XTickLabel',condNames);
                
                % split by taskConds (shared)
                condNames=S.condNames(S.StudyNum==study);
                figure(3)
                barplot(RR.splitby,RR.Rcv,'subset',strcmp(RR.method,method{1}) & strcmp(RR.split,'shared') & strcmp(RR.trainMode,'crossed') & strcmp(RR.trainMode,'crossed'))
                ylabel('R')
                title(sprintf('%s-%s-shared',signal,method{1}))
                set(gca,'YLim',[0 0.55],'FontSize',8,'XTickLabel',condNames);
            case 'RcvRnc'
                % split Rcv versus Rnc
                figure(1)
                lineplot(RR.lambda(:,M1),RR.Rcv,'CAT',CAT,'subset',strcmp(RR.method,method{1}) & strcmp(RR.split,'unique') & strcmp(RR.trainMode,'crossed'))
                hold on
                CAT.linecolor={'k'};
                CAT.linestyle={'-'};
                CAT.errorcolor={'k'};
                lineplot(RR.lambda(:,M1),RR.Rnc,'CAT',CAT,'subset',strcmp(RR.method,method{1}) & strcmp(RR.split,'unique') & strcmp(RR.trainMode,'crossed'))
                ylabel('R')
                xlabel('lambda')
                title(sprintf('%s-%s',signal,method{1}))
                set(gca,'YLim',[.1 0.3],'FontSize',14,'ytick',[.1,.15,.2,.25,.3]);
            case 'bestMethod'
                % bestMethod
                xyplot(RR.lambda(:,M1),RR.Rcv,[RR.lambda(:,M1)],'split',[RR.methodNum],...
                    'subset',strcmp(RR.split,'unique') & strcmp(RR.trainMode,'crossed'),...
                    'style_thickline','markersize',10,'leg','auto','markersize',10);
                ylabel('R')
                xlabel('lambdas')
                
                A=getrow(RR,strcmp(RR.split,'unique') & strcmp(RR.trainMode,'crossed'));
                idx=find(A.Rcv==max(A.Rcv));
                bestLambda=A.lambda(idx,:);
                fprintf('bestLambda is %2.2f \n',bestLambda(2))
%                 set(gca,'FontSize',14,'xtick',[]);
            case 'bestSubj'
                for s=1:length(returnSubjs),
                    A=getrow(RR,RR.SN==returnSubjs(s) & strcmp(RR.split,'unique') & strcmp(RR.trainMode,'crossed'));
                    idx=find(A.Rcv==max(A.Rcv));
                    bestLambda=A.lambda(idx,:);
                    fprintf('bestLambda for s%d is %2.2f \n',returnSubjs(s),bestLambda(2))
                end
        end
        keyboard;
    case 'CONNECT:bestLambda'
        method=varargin{1}; % 'L2' etc
        study=varargin{2}; % evaluating on which study, 1 or 2 ? [1,2]
        
        signal='fitted';
        model='162_tessellation_hem';
        trainMode='crossed';
        condType='unique';
        
        RR=[];
        
        T=load(fullfile(studyDir{2},connDir,'glm4','eval',sprintf('%s_%s_%s_sc%d_%s.mat',model,signal,method,study,condType)));

        A=getrow(T,strcmp(T.trainMode,trainMode));
        idx=find(A.Rcv==max(A.Rcv));
        bestLambda=A.lambda(idx,:);
        fprintf('bestLambda is %2.2f \n',bestLambda(2))
        
        varargout={bestLambda}; 
        
    case 'MAP:cortex'              % Map of where projections come from
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
    case 'MAP:plotCortex'
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
    case 'MAP:cerebellum'          % Maps certain stats from an individual to SUIT and then to flat map
        data    = varargin{1};
        sn      = varargin{2};
        
        % Determine the voxels we want to resample in SUIT space
        V=spm_vol(fullfile(studyDir{1},suitDir,'anatomicals','cerebellarGreySUIT.nii'));
        X= spm_read_vols(V);
        grey_threshold = 0.1; % gray matter threshold
        linIn1=find(X>grey_threshold);
        [i1,j1,k1]= ind2sub(V.dim,linIn1');
        [x1,y1,z1]= spmj_affine_transform(i1,j1,k1,V.mat);
        
        % Determine voxel locations from the original ROI
        load(fullfile(studyDir{1},regDir,'data',subj_name{sn},'regions_cerebellum_grey.mat')); % 'regions' are defined in 'ROI_define'
        Vmask = spm_vol(fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn},'maskbrainSUITGrey.nii'));
        [i3,j3,k3]=spmj_affine_transform(R{1}.data(:,1),R{1}.data(:,2),R{1}.data(:,3),inv(Vmask.mat));
        linIn3 = sub2ind(Vmask.dim,round(i3),round(j3),round(k3));
        
        % transform SUIT coords into anatomical space of the individual
        flowfield = fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn},'u_a_c_anatomical_seg1.nii');
        affine    = fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn},'Affine_c_anatomical_seg1.mat');
        [Def,Aff]=spmdefs_get_dartel(flowfield,affine);
        [x2,y2,z2]=spmdefs_transform(Def,Aff,x1,y1,z1);
        [i2,j2,k2]=spmj_affine_transform(x2,y2,z2,inv(Vmask.mat));
        
        % resample the weights into SUIT space
        for r=1:size(data,1)
            Vout=Vmask;
            Vout.dat = zeros(Vout.dim);
            Vout.dat(linIn3)=data(r,:);
            Vout.dat=data(r,:);
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
    case 'MAP:suit_reslice'
        data=varargin{1};
        sn=varargin{2}; 
        
        load(fullfile(studyDir{1},regDir,'data',subj_name{sn},'regions_Cerebellum_grey.mat')); 
        
        % Determine the voxels we want to resample in SUIT space
        V=spm_vol(fullfile(studyDir{1},suitDir,'anatomicals','cerebellarGreySUIT.nii'));
        X= spm_read_vols(V);
        grey_threshold = 0.1; % gray matter threshold
        linIn1=find(X>grey_threshold);
        [i1,j1,k1]= ind2sub(V.dim,linIn1');
        [x1,y1,z1]= spmj_affine_transform(i1,j1,k1,V.mat);

        % load cerebellar mask in individual func space
        Vi=spm_vol(fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn},'maskbrainSUITGrey.nii'));
        X=spm_read_vols(Vi);
        indx=find(X>0);
        
        % make volume
        for b=1:size(data,1),
            Yy=zeros(1,Vi.dim(1)*Vi.dim(2)*Vi.dim(3));
            Yy(1,indx)=data(b,:);
            Yy=reshape(Yy,[Vi.dim(1),Vi.dim(2),Vi.dim(3)]);
            Yy(Yy==0)=NaN;
            Vi.fname=fullfile(studyDir{2},'connectivity','glm4','eval',sprintf('temp_cereb_%2.4d.nii',b));
            spm_write_vol(Vi,Yy);
            clear Yy
            filenames{b}=Vi.fname;
        end
        % reslice data into suit space
        job.subj.affineTr = {fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn},'Affine_c_anatomical_seg1.mat')};
        job.subj.flowfield= {fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn},'u_a_c_anatomical_seg1.nii')};
        job.subj.resample = filenames';
        job.subj.mask     = {fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn},'cereb_prob_corr_grey.nii')};
        job.vox           = [2 2 2];
        job.outFile       = 'mat';
        D=suit_reslice_dartel(job);
        % delete temporary files
        deleteFiles=dir(fullfile(studyDir{2},'connectivity','eval',subj_name{sn},'*temp*'));
        for b=1:length(deleteFiles),
            delete(char(fullfile(studyDir{2},'connectivity','eval',subj_name{sn},deleteFiles(b).name)));
        end
        
        % Now map to surface-based representation
        D = suit_map2surf(Vres,'stats','mean');
        varargout={D,Vres};
    case 'MAP:plotCerebellum'      % Plots weights or evaluations onto flatmap for subj or group
        method=varargin{1}; % 'L2' or 'nonNegative'
        study=varargin{2}; % evaluating on which study, 1 or 2 ? [1,2]
        
        % set up parameters
        l1=0;
        l2=0;
        signal='fitted';
        model='162_tessellation_hem';
        data='Cerebellum_grey';
        trainMode='crossed';
        condType='unique'; 
        toPlot='eval'; % other options - 'eval'
        splitby='no'; % get pred maps for each task ?
        sn=[];
        
        vararginoptions(varargin(3:end),{'l1','l2','signal','model','data','trainMode','condType','splitby','sn'}); % 'condType' is either 'unique or 'all'

        % get best lambda
        bestLambda=sc1_sc2_CCC('CONNECT:bestLambda',method,study);
        
        sc1_sc2_CCC('CONNECT:eval',returnSubjs,'L2',1,'l2',bestLambda(2),'voxels',1)
        
        switch toPlot,
            case 'weights'
                weightDir=fullfile(studyDir{2},connDir,'glm4','weights',method);dircheck(weightDir);
                M=load(fullfile(weightDir,sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d.mat',model,signal,method,study,round(bestLambda(1)),round(bestLambda(2)))));
            case 'eval'
                 evalDir=fullfile(studyDir{2},connDir,'glm4','eval');dircheck(evalDir);
                 T=load(fullfile(evalDir,sprintf('%s_%s_%s_sc%d_%s_voxels.mat',model,signal,method,study,condType)));     
        end
        
        subjs=unique(T.SN);
        taskSplits=unique(T.splitby);
        % loop over subjects
        for s=1:length(subjs),
            
            A=getrow(T,T.SN==subjs(s) & strcmp(T.trainMode,trainMode) & strcmp(T.split,condType));
            
            for sp=1:length(taskSplits),
                % get flatmap coordinates
                [D,Vres]=sc1_sc2_CCC('MAP:suit_reslice',A.Rcv{taskSplits(sp)},subjs(s));
            end
        end
        keyboard;
 
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

