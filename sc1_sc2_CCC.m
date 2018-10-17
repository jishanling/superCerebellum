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
studyStr        = {'SC1','SC2','SC12'};
behavDir        =['data'];
imagingDir      =['imaging_data'];
suitDir         =['suit'];
caretDir        =['surfaceCaret'];
regDir          =['RegionOfInterest/'];
connDir         =['connectivity'];
encodeDir       ='/encoding';

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
exper = {'sc1','sc2'};
glm   = 4;

%==========================================================================

switch(what)
    
    case 'ACTIVITY:get_data'
        sn=varargin{1}; % Subj numbers to include
        study=varargin{2}; % 1 or 2 or [1,2]
        type=varargin{3}; % 'build' or 'eval'. For build - we get group data. For eval - we get indiv data
        hemi=varargin{4}; % 'LeftHem', or 'RightHem'
        
        vararginoptions({varargin{5:end}},{'sess'}); % fracture further into sessions [either 1 or 2]
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % input could be 'group' or <subjNums>
        if strcmp(sn,'group'),
            sn=returnSubjs;
        end
        
        % left or right hemi ?
        if strcmp(hemi,'LeftHem'),
            hem='lh';
        else
            hem='rh';
        end
        
        % load data
        UFullAvrgAll=[];
        for f=1:length(study),
            load(fullfile(studyDir{study(f)},encodeDir,'glm4',sprintf('cortex_%s_avrgDataStruct_vert.mat',hem)));
            
            D1=getrow(D,D.StudyNum==f);
            idx=D1.condNum(D1.overlap==1); % get index for unique tasks
            % format data
            for s=1:length(sn),
                if exist('sess'),
                    indx=T.SN==sn(s) & T.sess==sess;
                    UFullAvrg=T.data(indx,:);
                else
                    indx=T.SN==sn(s);
                    UFull=T.data(indx,:);
                    [numConds,numVox]=size(UFull);
                    % get average across sessions
                    numConds=numConds/2;
                    for c=1:numConds,
                        UFullAvrg(c,:,s)=nanmean(UFull([c,c+numConds],:),1);
                    end
                end
            end
            
            switch type,
                case 'build'
                    % if group - get mean
                    UFull=nanmean(UFullAvrg,3);
                    % remove mean of shared tasks
                    UFullAvrg_C=bsxfun(@minus,UFull,mean(UFull(idx,:)));
                case 'eval'
                    % remove mean of shared tasks
                    UFullAvrg_C=bsxfun(@minus,UFullAvrg,mean(UFullAvrg(idx,:,:)));
            end
            
            % if func1+func2 - concatenate
            if length(study)>1,
                UFullAvrgAll=[UFullAvrgAll;UFullAvrg_C];
            else
                UFullAvrgAll=UFullAvrg_C;
            end
        end
        
        % center the data (remove overall mean)
        X_C=bsxfun(@minus,UFullAvrgAll,mean(UFullAvrgAll));
        varargout={X_C,sn};
    case 'ACTIVITY:make_model' % make X matrix (feature models)
        study=varargin{1}; % 1, 2, or [1,2]
        type=varargin{2};  % 'yes' or 'no' to modelling each run separately ?
        
        F=dload(fullfile(baseDir,'motorFeats.txt')); % load in motor features
        
        % sort out which study we're taking (or both) ?
        if length(study)>1,
            Fs=F;
        else
            Fs=getrow(F,F.studyNum==study);
        end
        numConds=length(Fs.studyNum);
        
        % make feature model
        x=[eye(numConds) Fs.lHand./Fs.duration Fs.rHand./Fs.duration Fs.saccades./Fs.duration];
        featNames=Fs.condNames;
        featNames{numConds+1}='lHand';
        featNames{numConds+2}='rHand';
        featNames{numConds+3}='saccades';
        
        switch type,
            case 'yes'
                % make cond x run x features
                for f=1:size(x,2),
                    X.x(:,f)=repmat(x(:,f),numel(run),1);
                    X.idx=repmat(Fs.condNum,numel(run),1);
                end
            case 'no'
                % do nothing
                X.x=x;
                X.idx=[1:size(x,1)]';
        end
        
        % normalise features
        X.x   = bsxfun(@minus,X.x,mean(X.x));
        X.x   = bsxfun(@rdivide,X.x,sqrt(sum(X.x.^2)));  % Normalize to unit length vectors
        
        % add rest
        if strcmp(type,'yes')
            % add rest to the end (better way of doing this !)
            rest=X.x(X.idx==numConds,:);
            X.x(X.idx==numConds,:)=[];
            X=[X.x; rest];
        else
            X=X.x;
        end
        
        varargout={X,featNames,numConds,F};
    case 'ACTIVITY:patterns'   % estimate each of the task conditions (motor features are subtracted)
        sn=varargin{1}; % 'group' or <subjNum>
        study=varargin{2}; % 1 or 2 ?
        taskType=varargin{3}; % 'allConds', 'averageConds','averageTasks'
        hemi=varargin{4}; % 'LeftHem' or 'RightHem'
        
        lambda=.01;
        
        % group or individual ?
        if ~strcmp(sn,'group'),
            % load in map
            if length(study)>1,
                outDir=fullfile(studyDir{2},caretDir,sprintf('x%s/%s',subj_name{sn},hemi)); dircheck(outDir)
            else
                outDir=fullfile(studyDir{study},caretDir,sprintf('x%s/%s',subj_name{sn},hemi)); dircheck(outDir);
            end
            subjs=sn;
        else
            if length(study)<2,
                outDir=fullfile(studyDir{study},caretDir,'fsaverage_sym',hemi,'glm4');
            else
                outDir=fullfile(studyDir{2},caretDir,'fsaverage_sym',hemi,'glm4');
            end
            subjs=returnSubjs;
        end
        
        % load in task structure file
        F=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % load in activity patterns
        [data]=sc1_sc2_CCC('ACTIVITY:get_data',returnSubjs,study,'eval',hemi);
        
        % get feature model
        [X,featNames,numConds]=sc1_sc2_CCC('ACTIVITY:make_model',study,'no'); % load in model
        
        % rest
        if length(study)>1,
            rest=[29,61];
        else
            rest=numConds;
        end
        numFeat=size(X,2)-numConds;
        
        % we're not regressing out motor features against rest
        X(rest,numConds+1:numConds+3)=X(rest,numConds+1:numConds+3)*-1;
        
        % regress out motor features
        for s=1:length(subjs),
            B(:,s,:)=(X'*X+eye(numConds+numFeat)*lambda)\(X'*data(:,:,s)); % was subjs(s) for some reason ?
            fprintf('ridge regress done for subj%d done \n',returnSubjs(s))
        end;
        clear data
        
        % subtract baseline
        baseline=nanmean(B,1);
        B=bsxfun(@minus,B,baseline);
        
        % z score the activity patterns
        B=zscore(B);
        
        % group or individual ?
        if strcmp(sn,'group'),
            indices=permute(B,[3 1 2]);
            indices=nanmean(indices,3);
        else
            indices=permute(B,[3 1 2]);
        end
        
        % make metric file
        S=caret_struct('metric','data',indices,'column_name',featNames);
        
        switch taskType,
            case 'allConds'
                condName='unCorr_allTaskConds';
            case 'averageConds'
                condNumUni=[F.condNumUni;62;63;64];
                X1=indicatorMatrix('identity_p',condNumUni);
                uniqueTasks=S.data*X1; % try pinv here ?
                % get new condNames (unique only)
                condNames=[F.condNames(F.StudyNum==1);F.condNames(F.StudyNum==2 & F.overlap==0)];
                condNames{length(condNames)+1}='lHand';
                condNames{length(condNames)+1}='rHand';
                condNames{length(condNames)+1}='saccades';
                S.data=uniqueTasks;
                S.column_name=condNames';
                S.num_cols=size(S.column_name,2);
                S.column_color_mapping=S.column_color_mapping(1:S.num_cols,:);
                condName='unCorr_avrgTaskConds'; % average of certain tasks
        end
        
        % save out metric
        if strcmp(sn,'group'),
            outName=condName;
        else
            outName=sprintf('%s_%s',subj_name{sn},condName);
        end
        caret_save(fullfile(outDir,sprintf('%s.metric',outName)),S);
        %
        %         varargout={B,featNames};
    case 'ACTIVITY:reliability_overall'
        glm=varargin{1};
        type=varargin{2}; % 'cerebellum' or 'cortex' 'basalGanglia'
        study=varargin{3}; % studyNum = 1 or 2
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        load(fullfile(studyDir{2},encodeDir,sprintf('glm%d',glm),sprintf('allVox_sc1_sc2_sess_%s.mat',type)));
        Y=Yy;clear Yy;
        
        numSubj=length(Y);
        
        S=[];
        for subj=1:numSubj,
            R=[];
            idx=1;
            D1=getrow(D,D.StudyNum==study);
            %                 sharedConds=D1.condNumUni.*D1.overlap;
            overallConds=D1.condNumUni.*D1.StudyNum;
            if study==1,
                condNames=D1.condNames(find(overallConds));
            end
            % sharedConds=sharedConds(randperm(numel(sharedConds{2}))); % Shuffle
            cN=condNum{study}-1;  % Important: In the allVox file, instruction is still included!
            pN=partNum{study};    % Partition Numner
            sN=(pN>8)+1;          % Sessions Number
            for se=1:2,
                X1=indicatorMatrix('identity_p',cN.*(sN==se));  % This one is the matrix that related trials-> condition numbers
                X2=indicatorMatrix('identity_p',overallConds); % THis goes from condNum to shared condNumUni
                Yf(:,:,idx,subj)=pinv(X1*X2)*Y{subj}{study};
                Yf(:,:,idx,subj)=bsxfun(@minus,Yf(:,:,idx,subj),nanmean(Yf(:,:,idx,subj)));
                idx=idx+1;
            end;
            CORRMatrix=corr(Yf(:,:,1,subj),Yf(:,:,2,subj));
            CORR(:,subj)=diag(CORRMatrix);
            %             for c=1:size(Yf,1),
            %                 CORR(c,:,:,subj)=interSubj_corr_voxel(Yf(c,:,:,subj));
            %                 T.SN      = returnSubjs(subj);
            %                 T.within1 = CORR(c,1,2,subj);
            %                 T.within2 = CORR(c,3,4,subj);
            %                 T.across  = nanmean(nanmean(CORR(c,1:2,3:4,subj)));
            %                 T.condNum = c;
            %                 T.condNames={condNames{c}};
            %                 R=addstruct(R,T);
            %                 clear T
            %             end;
            fprintf('subj%d done',returnSubjs(subj));
        end;
        save(fullfile(studyDir{study},regDir,'glm4','patternReliability_voxel.mat'),'CORR')
    case 'PLOT:reliabilityA'
        % load relability
        load(fullfile(studyDir{2},regDir,'glm4','patternReliability_cerebellum.mat'));
        
        %         figure();lineplot(S.condNum,[S.within1,S.within2,S.across],'leg',{'within1','within2','across'})
        A=tapply(S,{'SN'},{'across'},{'within1'},{'within2'});
        
        % within & between-dataset reliability
        myboxplot([],[A.within1 A.within2 A.across],'style_twoblock','plotall',1);
        drawline(0,'dir','horz');
        ttest(sqrt(A.within1.*A.within2),A.across,2,'paired');
        
        x1=nanmean(A.within1);x2=nanmean(A.within2);x3=nanmean(A.across);
        SEM1=std(A.within1)/sqrt(length(returnSubjs));SEM2=std(A.within2)/sqrt(length(returnSubjs));SEM3=std(A.across)/sqrt(length(returnSubjs));
        fprintf('average corr for set A is %2.3f; CI:%2.3f-%2.3f \n average corr for set B is %2.3f; CI:%2.3f-%2.3f and average corr across sets A and B is %2.3f; CI:%2.3f-%2.3f \n',...
            x1,x1-(1.96*SEM1),x1+(1.96*SEM1),x2,...
            x2-(1.96*SEM2),x2+(1.96*SEM2),...
            x3,x3-(1.96*SEM3),x3+(1.96*SEM3));
    case 'PLOT:reliability_voxel'
        
        for sess=1:2,
            load(fullfile(studyDir{sess},regDir,'glm4','patternReliability_voxel.mat'))
            data(:,sess)=nanmean(CORR,2);
            clear CORR
        end
        
        % get average across task sets
        data_average=nanmean(data,2);
        data_average=data_average';
        
        [~,~,~,idx]=sc1_sc2_functionalAtlas('PREDICTIONS:datasets','R','run');
        
        V=spm_vol(fullfile(studyDir{1},suitDir,'anatomicals','cerebellarGreySUIT.nii'));
        C{1}.dim=V.dim;
        C{1}.mat=V.mat;
        
        Yy=zeros(size(data_average,1),V.dim(1)*V.dim(2)*V.dim(3));
        
        % make vol
        Yy(:,idx)=data_average;
        
        % get avrg across subjs
        indices=nanmean(Yy,1);
        
        % map vol2surf
        data=reshape(indices,[V.dim(1),V.dim(2),V.dim(3)]);
        C{1}.dat=data;
        
        M=caret_suit_map2surf(C,'space','SUIT','stats','nanmean');  % MK created caret_suit_map2surf to allow for output to be used as input to caret_save
        
        % save out metric
        caret_save(fullfile(studyDir{2},caretDir,'suit_flat','glm4','voxel_reliability.metric'),M);
        
    case 'RELIABILITY:get_spatialFreq'
        study=varargin{1};
        frequencyBands  = [0 0.5 1 1.5 2 inf];
        load(fullfile(studyDir{study},'encoding','glm4','cereb_avrgDataStruct.mat'));
        RR=[];
        sn=unique(T.SN);
        for s = 1:length(sn)
            fprintf('subject %d\n',sn(s));
            for se=1:2
                S=getrow(T,T.SN == sn(s) & T.sess==se);
                for c=1:max(T.cond)
                    X=zeros(V.dim);
                    X(volIndx)=S.data(c,:);
                    X(isnan(X))=0;
                    % Y=mva_frequency3D(X,frequencyBands,'Voxelsize',[2 2 2],'plotSlice',15);
                    Y=mva_frequency3D(X,frequencyBands,'Voxelsize',[2 2 2]);
                    R=getrow(S,c);
                    for f=1:size(Y,4);
                        YY=Y(:,:,:,f);
                        R.data=YY(volIndx);
                        R.freq = f;
                        R.freqLow = frequencyBands(f);
                        RR=addstruct(RR,R);
                    end;
                end;
            end;
        end;
        save(fullfile(studyDir{study},'encoding','glm4','cereb_avrgDataStruct_freq.mat'),'-struct','RR');
    case 'RELIABILITY:spatialFreqCorr'
        study = [1 2]; % %experiment
        glm   = 'glm4';
        
        vararginoptions(varargin,{'study','glm'});
        C=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        A=[];
        for e=study
            T=load(fullfile(studyDir{e},'encoding',glm,'cereb_avrgDataStruct_freq.mat'));
            D=load(fullfile(studyDir{e},'encoding',glm,'cereb_avrgDataStruct.mat'));
            D.T.freq = ones(length(D.T.SN),1)*0;
            D.T.freqLow = -1*ones(length(D.T.SN),1);
            T=addstruct(T,D.T);
            T.study = ones(length(T.SN),1)*e;
            % Recenter Data and combine
            commonCond = C.condNum(C.StudyNum==e & C.overlap==1);
            for sn = unique(T.SN)'
                for se = unique(T.sess)'
                    for f = unique(T.freq)'
                        i = find(T.SN==sn & T.sess==se & T.freq==f);
                        j = find(T.SN==sn & T.sess==se & T.freq==f & ismember(T.cond,commonCond));
                        T.data(i,:)=bsxfun(@minus,T.data(i,:),nanmean(T.data(j,:)));
                    end;
                end;
            end;
            A=addstruct(A,T);
        end;
        
        % before - code below was computing corr on study 2 only (structure
        % was T instead of A)
        D=[];
        sn=unique(A.SN);
        numSubj = length(sn);
        SS=[];
        for st=1:2, % loop over studies
            RR=[];
            for f=unique(A.freq)', % loop over frequencies
                for s = 1:numSubj  % loop over subjects
                    for se=1:2     % loop over sessions
                        temp = A.data(A.study==st & A.SN==sn(s) & A.sess==se & A.freq==f,:);
                        %                         temp = bsxfun(@minus,temp,mean(temp));
                        D(:,:,s+(se-1)*length(sn))=temp;
                    end;
                end;
                C=intersubj_corr(D);
                R.sess = [ones(1,numSubj) ones(1,numSubj)*2]';
                R.subj = [1:numSubj 1:numSubj]';
                R.subj = sn(R.subj);
                R.freq = f*ones(numSubj*2,1);
                R.study= st*ones(numSubj*2,1);
                SameSess = bsxfun(@eq,R.sess',R.sess);
                SameSubj = bsxfun(@eq,R.subj',R.subj);
                for i=1:numSubj*2;
                    R.withinSubj(i,:)=C(i,SameSubj(i,:) & ~SameSess(i,:));
                    R.betweenSubj(i,:)=mean(C(i,~SameSubj(i,:)));
                    R.totSS(i,1) = nansum(nansum(D(:,:,i).^2));
                end;
                RR=addstruct(RR,R);
            end;
            clear temp D R
            SS=addstruct(SS,RR);
        end
        save(fullfile(studyDir{2},'encoding','glm4','cereb_spatialCorr_freq.mat'),'-struct','SS');
        varargout={SS};
    case 'PLOT:spatialFreqCorr'
        CAT=varargin{1};
        
        % load in spatialCorrFreq struct
        T=load(fullfile(studyDir{2},'encoding','glm4','cereb_spatialCorr_freq.mat'));
        
        xlabels={'overall','0-0.5','0.5-1','1-1.5','1.5-2','>2'};
        
        T=tapply(T,{'subj','freq'},{'withinSubj'},{'betweenSubj'},{'totSS'},'subset',ismember(T.subj,returnSubjs));
        T.freqK = T.freq>0;
        [ss,sn]=pivottable(T.subj,[],T.totSS,'mean','subset',T.freq==0);
        a(sn,1)=ss;
        T.relSS=T.totSS./a(T.subj);
        lineplot([T.freqK T.freq],[T.relSS],'CAT',CAT);
        set(gca,'XTickLabel',xlabels,'YLim',[0 0.35]);
        ylabel('Relative Power');
        xlabel('Cycles/cm');
        title('Relative amount of Variance per Frequency band');
    case 'PLOT:interSubjCorr'
        CAT=varargin{1};
        
        % load in spatialCorrFreq struct
        T=load(fullfile(studyDir{2},'encoding','glm4','cereb_spatialCorr_freq.mat'));
        
        T=tapply(T,{'subj','freq'},{'withinSubj'},{'betweenSubj'},{'totSS'},'subset',ismember(T.subj,returnSubjs));
        T.freqK = T.freq>0;
        lineplot([T.freqK T.freq],[T.withinSubj T.betweenSubj],'CAT',CAT,'leg','auto');
        
    case 'REPRESENTATION:get_distances'
        type=varargin{1}; % 'cerebellum'
        removeMotor=varargin{2}; % 'hands','saccades','all','none'
        taskType=varargin{3}; % 'unique' or 'all' task conditions ?
        
        load(fullfile(studyDir{2},regDir,'glm4',sprintf('G_hat_sc1_sc2_%s.mat',type)))
        subjs=size(G_hat,3);
        
        % load in condName info
        T=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % load feature matrix
        F=dload(fullfile(baseDir,'motorFeats.txt')); % load in motor features
        
        % get unique tasks
        switch taskType,
            case 'unique'
                X1=indicatorMatrix('identity_p',T.condNumUni);
                for s=1:subjs,
                    G(:,:,s)=X1'*(G_hat(:,:,s)*X1);
                end
                condNames=[T.condNames(T.StudyNum==1);T.condNames(T.StudyNum==2 & T.overlap==0)];
            case 'all'
                G=G_hat;
                condNames=T.condNames;
        end
        numDist=size(G,1);
        
        switch removeMotor,
            case 'all'
                X   = [F.lHand./F.duration F.rHand./F.duration F.saccades./F.duration];
            case 'none'
                X   = [];
        end
        
        % get unique taskConds
        if strcmp(taskType,'unique'),
            X   = pivottablerow(T.condNumUni,X,'mean(x,1)');
        end
        
        X   = [X eye(numDist)];
        X   = bsxfun(@minus,X,mean(X));
        X   = bsxfun(@rdivide,X,sqrt(sum(X.^2)));  % Normalize to unit length vectors
        
        % Get RDM
        for s=1:subjs,
            H=eye(numDist)-ones(numDist)/numDist; % centering matrix
            G(:,:,s)=H*G(:,:,s)*H'; % subtract out mean pattern
            IPM=rsa_vectorizeIPM(G(:,:,s));
            con = indicatorMatrix('allpairs',[1:numDist]);
            N = rsa_squareIPM(IPM);
            D = rsa.rdm.squareRDM(diag(con*N*con'));
            fullRDM(:,:,s) = D;
        end
        
        varargout={fullRDM,condNames,X,taskType};
    case 'REPRESENTATION:reliability'
        glm=varargin{1};
        type=varargin{2}; % 'cerebellum'
        % example 'sc1_sc2_imana('CHECK:DIST',4,'cerebellum')
        
        load(fullfile(regDir,sprintf('glm%d',glm),sprintf('G_hat_sc1_sc2_sess_%s.mat',type)));
        D=dload('sc1_sc2_taskConds.txt');
        D1=getrow(D,D.StudyNum==1);
        D2=getrow(D,D.StudyNum==2);
        
        % Look at the shared conditions only
        i1 = find(D1.overlap==1);
        i2 = find(D2.overlap==1);
        [~,b] = sort(D2.condNumUni(i2));     % Bring the indices for sc2 into the right order.
        i2=i2(b);
        numCond = length(i1);
        numSubj = size(G_hat_sc1,4);
        numSess = 2;
        
        C=indicatorMatrix('allpairs',[1:numCond]);
        for i=1:numSubj
            for j=1:numSess,
                dist(:,j  ,i)  = ssqrt(diag(C*G_hat_sc1(i1,i1,j,i)*C'));
                dist(:,j+2,i)  = ssqrt(diag(C*G_hat_sc2(i2,i2,j,i)*C'));
            end;
            CORR(:,:,i)    = corr(dist(:,:,i));
            T.SN(i,1)      = i;
            T.within1(i,1) = CORR(1,2,i);
            T.within2(i,1) = CORR(3,4,i);
            T.across(i,1)  = mean(mean(CORR(1:2,3:4,i)));
        end;
        
        save(fullfile(studyDir{2},regDir,'glm4',sprintf('distanceReliability_%s.mat'),type),'T','dist')
    case 'REPRESENTATION:RDM'
        taskType=varargin{1}; % 'unique' or 'all' tasks
        
        threshold=.001;
        
        condNames={'1.No-Go','2.Go','3.Theory of Mind','4.Action Observation','5.Video Knots','6.IAPS Unpleasant',...
            '7.IAPS Pleasant','8.Math','9.Digit Judgement','10.Objects','11.IAPS Sad','12.IAPS Happy','13.Interval Timing',...
            '14.Motor Imagery','15.Finger Simple','16.Finger Sequence','17.Verbal 0Back','18.Verbal 2Back','19.Object 0Back',...
            '20.Object 2Back','21.Spatial Navigation','22.Stroop Incongruent','23.Stroop Congruent','24.Verb Generation',...
            '25.Word Reading','26.Visual Search Small','27.Visual Search Medium','28.Visual Search Hard','29.Rest','30.CPRO','31.Prediction',...
            '32.Prediction Violated','33.Prediction Scrambled','34.Spatial Map Easy','35.Spatial Map Medium','36.Spatial Map Hard',...
            '37.Nature Movie','38.Animated Movie','39.Landscape Movie','40.Mental Rotation Easy','41.Mental Rotation Medium',...
            '42.Mental Rotation Hard','43.Biological Motion','44.Scrambled Motion','45.Response Alt Easy','46.Response Alt Medium','47.Response Alt Hard'};
        
        % load in fullRDM
        [fullRDM,~]=sc1_sc2_functionalAtlas('REPRESENTATION:get_distances','cerebellum','all',taskType); % remove all motorFeats
        
        numDist=size(fullRDM,1);
        
        % Plot RDM
        switch taskType,
            case 'unique',
                reOrder=[1:5,8:10,14:18,21,24,25,29,30,34:39,6,7,11,12,19,20,13,31:33,26:28,22,23,40:47];
            case 'all',
                reOrder=[1,2,6,7,8,9,10,11,12,13,14,17,18,22,23,3,4,5,15,16,19,20,21,24,25,26,...
                    27,28,29,58,59,60,43,44,49,48,36,34,35,55,56,57,61,30,31,32,33,37,38,39,40,41,42,45,46,47,50,51,52,53,54]'; % reorder
        end
        
        % reorder RDM
        fullRDM=fullRDM(reOrder,reOrder,:);
        
        % threshold RDM
        fullRDM_thresh=reshape(fullRDM,[size(fullRDM,1)*size(fullRDM,2)],[]);
        for dd=1:size(fullRDM_thresh,1),
            [t(dd),p(dd)]=ttest(fullRDM_thresh(dd,:),[],1,'onesample');
        end
        ut=t; % unthresholded distances
        
        % uncorrected p-vals
        t(p>threshold)=0; % thresholded distances
        
        % zero the nan values
        t(isnan(t))=0;
        ut(isnan(ut))=0;
        
        fprintf('%2.2f%% of the pairwise distances are significantly different from zero \n',(length(t(t>0))/length(t))*100);
        
        % visualise thresholded RDM (t-values)
        squareT=tril(reshape(t,[size(fullRDM,1),size(fullRDM,2)]));
        squareUT=reshape(ut,[size(fullRDM,1),size(fullRDM,2)]);
        idxUT=find(triu(squareUT));
        squareT(idxUT)=squareUT(idxUT);
        figure();imagesc_rectangle(abs(squareT),'YDir','reverse')
        caxis([0 1]);
        g=set(gca,'Ytick',[1:numDist]','YTickLabel',condNames(reOrder),'FontSize',7);
        g.Color='white';
        colorbar
    case 'REPRESENTATION:MDS'
        taskType=varargin{1}; % 'unique' or 'all' tasks
        clustering=varargin{2}; % 'distance' or 'region'
        
        % colour
        colour={[1 0 0],[0 1 0],[0 0 1],[0.3 0.3 0.3],[1 0 1],[1 1 0],[0 1 1],...
            [0.5 0 0.5],[0.8 0.8 0.8],[.07 .48 .84],[.99 .76 .21],[.11 .7 .68],...
            [.39 .74 .52],[.21 .21 .62],[0.2 0.2 0.2],[.6 .6 .6],[.3 0 .8],[.8 0 .4],...
            [0 .9 .2],[.1 .3 0],[.2 .4 0],[.63 0 .25],[0 .43 .21],[.4 0 .8]};
        
        vararginoptions({varargin{3:end}},{'CAT','colour'}); % option if doing individual map analysis
        
        % load in fullRDM
        [fullRDM,condNames,X,taskType]=sc1_sc2_CCC('REPRESENTATION:get_distances','162_tessellation_hem','all',taskType); % remove all motorFeats
        condIndx=1:length(condNames);
        
        % load in condName info
        T=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        avrgFullRDM=ssqrt(nanmean(fullRDM,3));
        
        vecRDM = rsa.rdm.vectorizeRDM(avrgFullRDM);
        [Y,~] = rsa_classicalMDS(vecRDM,'mode','RDM');
        B = (X'*X+eye(size(X,2))*0.0001)\(X'*Y); % ridge regression
        Yr    = Y  - X(:,1:3)*B(1:3,:); % remove motor features
        
        % define cluster colour
        [V,L]   = eig(Yr*Yr');
        [l,i]   = sort(diag(L),1,'descend'); % Sort the eigenvalues
        V       = V(:,i);
        X       = bsxfun(@times,V,sqrt(l'));
        X = real(X);
        
        switch clustering,
            case 'distance'
                clustTree = linkage(Yr,'average');
                indx = cluster(clustTree,'cutoff',1);
            case 'region'
                load(fullfile(studyDir{2},'encoding','glm4','groupEval_SC12_10cluster','SNN.mat'));
                cmap=load(fullfile(studyDir{2},'encoding','glm4','groupEval_SC12_10cluster','colourMap.txt'));
                
                % assign each task to a cluster
                if strcmp(taskType,'unique'),
                    bestF=pivottablerow(T.condNumUni,bestF,'mean(x,1)');
                end
                
                [x,indx]=max(bestF,[],2);
                
                % set threshold for tasks close to zero
                indx(x<.155)=11;
                
                colour=num2cell(cmap(:,2:4)/255,2)';
                colour{11}=[.8275 .8275 .8275]; % set grey
        end
        
        CAT.markercolor= {colour{indx}};
        CAT.markerfill = {colour{indx}};
        CAT.labelcolor = {colour{indx}};
        
        X1=X(condIndx,condIndx);
        figure()
        scatterplot3(X1(:,1),X1(:,2),X1(:,3),'split',condIndx','CAT',CAT,'label',condNames);
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'Box','on');
        hold on;
        plot3(0,0,0,'+');
        % Draw connecting lines
        %         for i=1:15,
        %             ind=clustTree(i,1:2);
        %             X(end+1,:)=(X(ind(1),:)+X(ind(2),:))/2;
        %             line(X(ind,1),X(ind,2),X(ind,3));
        %         end;
        hold off;
        clear X1  indxShort
        
    case 'MAP:optimal'     % figure out optimal map for multiple clusters
        % example:sc1_sc2_functionalAtlas('MAP:optimal',<subjNums>,1,6,'group')
        sn=varargin{1};     % 'group' or <subjNum>
        study=varargin{2};  % 1 or 2 or [1,2]
        K=varargin{3};      % K=numClusters (i.e. 5);
        hemi=varargin{4};   % 'LeftHem' or 'RightHem'
        
        numCount=5;         % How often the "same" solution needs to be found
        
        vararginoptions({varargin{5:end}},{'sess'}); % option if doing individual sessions
        
        tol_rand = 0.90;    % Tolerance on rand coefficient to call it the same solution
        maxIter=100; % if it's not finding a similar solution - force stop at 100 iters
        
        % Set the String correctly
        studyStr = sprintf('SC%d',study);
        if length(study)>1
            studyStr='SC12'; % both studies combined
        end
        
        % left or right hemi ?
        if strcmp(hemi,'LeftHem'),
            hem='lh';
        else
            hem='rh';
        end
        
        % Set output filename: group or indiv ?
        if strcmp(sn,'group'), % group
            sn=returnSubjs;
            if exist('sess'),
                outDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s_sess%d_%dcluster_%s',studyStr,sess,K,hem));
            else
                outDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s_%dcluster_%s',studyStr,K,hem));
            end
            dircheck(outDir);
            outName=fullfile(outDir,'SNN.mat');
        else % indiv
            outDir=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn});
            outName=fullfile(outDir,sprintf('SNN_%s_%dcluster.mat',studyStr,K));
        end
        
        % get data
        [X_C] = sc1_sc2_CCC('ACTIVITY:get_data',sn,study,'build',hemi);
        
        [F,G,Info,winner]=semiNonNegMatFac(X_C,K,'threshold',0.01); % get current error
        
        % Intialize iterations[G
        %         bestErr = inf;
        %         bestSol = ones(size(X_C,1),1);
        %         iter=1; % How many iterations
        %         count=0;
        %         while iter<maxIter,
        %             [F,G,Info,winner]=semiNonNegMatFac(X_C,K,'threshold',0.01); % get current error
        %             errors(iter)=Info.error;    % record error
        %             randInd(iter)=RandIndex(bestSol,winner); %
        %
        %             % Check if we have a similar solution
        %             if randInd(iter)>tol_rand % Similar solution
        %                 count=count+1;       % count up one
        %                 if (Info.error<bestErr)  % If we got slightly better - update
        %                     bestErr = Info.error;
        %                     bestSol = winner;
        %                     bestG   = G;
        %                     bestF   = F;
        %                     bestInfo = Info;
        %                 end;
        %             else                     % Different (enough) solution
        %                 if (Info.error<bestErr) % Is this a better solution
        %                     bestErr = Info.error;
        %                     bestSol = winner;
        %                     bestG   = G;
        %                     bestF   = F;
        %                     bestInfo = Info;
        %                     count = 0;         % first time we found this solution: reset counter
        %                 end;
        %             end;
        %             fprintf('Error: %2.2f Rand:%2.2f, Best:%2.2f currently found %d times\n',errors(iter),randInd(iter),bestErr,count);
        %             if count>=numCount || iter>=maxIter,
        %                 fprintf('Existing loop....\n');
        %                 break;
        %             end;
        %             iter=iter+1;
        %         end;
        save(outName,'F','G','winner','Info');
    case 'MAP:toSurf'
        % this function takes any labelled volume (already in SUIT space)
        % and plots to the surface
        K=varargin{1}; % number of clusters
        hem=varargin{2}; % 'lh' or 'rh'
        
        if strcmp(hem,'lh'),
            hemi='LeftHem';
        else
            hemi='RightHem';
        end
        
        inDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC12_%dcluster_%s',K,hem));
        inName=fullfile(inDir,'SNN.mat');
        
        load(inName);
        
        % load in medial wall
        M=caret_load(fullfile(studyDir{1},caretDir,'fsaverage_sym',hemi,sprintf('%s.medialWall.paint',hem)));
        
        % black out medial wall
        winner(M.data==1)=0;
        
        % make metric file
        S=caret_struct('paint','data',winner);
        
        inName=sprintf('%s.MDTB-%dcluster',hem,K);
        inDir=fullfile(studyDir{2},caretDir,'fsaverage_sym',hemi,'glm4');
        caret_save(fullfile(inDir,sprintf('%s.paint',inName)),S);
        
    case 'ENCODE:get_features'
        K=varargin{1}; % number of clusters
        hem=varargin{2}; % 'lh' or 'rh'
        
        D=dload(fullfile(baseDir,'featureTable_functionalAtlas.txt'));
        
        %         D=dload(fullfile(baseDir,'featureTable_jd_updated.txt')); % Read feature table - updated with new features "naturalistic bio motion" and "naturalistic scenes"
        S=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt')); % List of task conditions
        
        load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC12_%dcluster_%s',K,hem),'SNN.mat'));
        bestF=F;
        W=pivottablerow(S.condNumUni,bestF,'mean(x,1)'); % we only want unique taskConds
        
        % get new condNames (unique only)
        condNames=[S.condNames(S.StudyNum==1);S.condNames(S.StudyNum==2 & S.overlap==0)];
        
        % Make the feature matrix
        D.LeftHand    = D.leftHandPresses ./D.duration;
        D.RightHand   = D.rightHandPresses ./D.duration;
        D.Saccade    = D.saccades./D.duration;
        
        % remove superfluous
        D=rmfield(D,{'leftHandPresses','rightHandPresses','saccades','Imagination','LongtermMemory','SceneRecog'});
        %         D=rmfield(D,{'leftHandPresses','rightHandPresses','saccades'});
        
        f=fieldnames(D);
        FeatureNames = f(5:end);
        F=[];
        for d=1:length(FeatureNames)
            F = [F D.(FeatureNames{d})];
        end;
        F= bsxfun(@rdivide,F,sum(F.^2,1));
        numCond = length(D.conditionName);
        numFeat = length(FeatureNames);
        numClusters = size(W,2);
        
        lambda = [0.01 0.001];
        X=bsxfun(@minus,F,mean(F,1));
        
        Y=bsxfun(@minus,W,mean(W,1));
        X=bsxfun(@rdivide,X,sqrt(mean(X.^2)));
        Y=bsxfun(@rdivide,Y,sqrt(mean(Y.^2)));
        XX=X'*X;
        XY=X'*Y;
        A = -eye(numFeat);
        b = zeros(numFeat,1);
        
        for p=1:numClusters,
            %             u(:,p) = cplexqp(XX+lambda(2)*eye(numFeat),ones(numFeat,1)*lambda(1)-XY(:,p),A,b);
            u(:,p) = lsqnonneg(X,Y(:,p));
        end;
        
        % Get corr between feature weights
        C=corr(F,W);
        
        % Present the list of the largest three weights for each
        % cluster
        for i=1:numClusters,
            [a,b]=sort(u(:,i),'descend');
            B.clusters(i,1)=i;
            % get 3 highest corrs
            for f=1:3,
                B.featNames{i,f}=FeatureNames{b(f)};
                B.featIdx(i,f)=b(f);
                B.featCorrs(i,f)=a(f);
            end
            % what % do top 3 make up of overall features ?
            B.relSum(i,1)=(a(1)+a(2)+a(3))/sum(a)*100;
            B.relSuma(i,1)=(a(1))/sum(a)*100;
        end;
        
        %         fprintf('on average, %2.2f%% of all feature weights are accounted by the top 3 features \n with the top feature accounting for %2.2f %% \n',mean(B.relSum),mean(B.relSuma));
        fprintf('on average, %2.2f%% of all feature weights are accounted by the top 3 features \n with the top feature accounting for %2.2f %% \n',mean(B.relSum),mean(B.relSuma));
        
        varargout={B,F,W,u,condNames,FeatureNames,X,Y};
    case 'ENCODE:project_featSpace'
        K=varargin{1}; % number of clusters
        hem=varargin{2}; % 'lh' or 'rh'
        toPlot=varargin{3}; % 'winner' or 'all' or 'featMatrix'
        
        sizeWeight=40;
        
        % get features
        [B,F,~,C,condNames,FeatureNames]=sc1_sc2_CCC('ENCODE:get_features',K,hem);
        
        % get cluster colours
        cmap=load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC12_%dcluster_%s',K,hem),'colourMap.txt'));
        cmap=cmap/255;
        
        switch toPlot,
            case 'featList_all'
                % Make word lists for word map
                DD = C;WORD=FeatureNames;
                DD(DD<0)=0;
                DD=DD./max(DD(:));
                numClusters=size(C,2);
                numFeat=size(C,1);
                for j=1:numClusters,
                    subplot(2,5,j);
                    title(sprintf('Region %d',j),'Color',cmap(j,2:4),'FontSize',18)
                    set(gca,'Xticklabel',[],'Yticklabel',[])
                    for i=1:numFeat
                        if (DD(i,j)>.25)
                            siz=ceil(DD(i,j)*sizeWeight);
                            text(unifrnd(0,1,1),unifrnd(0,1,1),WORD{i},'FontSize',siz,'Color',cmap(j,2:4));
                        end;
                    end;
                end;
                %                 imagesc_rectangle(C);
                %                 caxis([0 1]);
                %                 t=set(gca,'Ytick',[1:length(FeatureNames)]','YTickLabel',FeatureNames');
                %                 t.Color='white';
                %                 colorbar
            case 'featList_winner'
                % figure out where to position features
                for i=1:size(B.featNames,1),
                    subplot(2,5,i);
                    set(gca,'XLim',[0 2.5],'YLim',[-0.2 1.2]);
                    text(0,1.2,sprintf('Region %d',i),'FontSize',20,'Color',cmap(i,2:4));
                    for j=1:size(B.featNames,2),
                        siz=ceil(B.featCorrs(i,j)*20);
                        text(unifrnd(0,1,1),unifrnd(0,1,1),B.featNames{i,j},'FontSize',siz,'Color',cmap(i,2:4));
                    end
                end
            case 'featMatrix'
                imagesc_rectangle(F');
                caxis([0 1]);
                t=set(gca,'Ytick',[1:length(FeatureNames)]','YTickLabel',FeatureNames',...
                    'Xtick',[1:length(condNames)]');
                t.Color='white';
        end
    case 'ENCODE:project_taskSpace'
        K=varargin{1}; % number of clusters
        hem=varargin{2}; % 'lh' or 'rh'
        toPlot=varargin{3}; % 'winner' or 'all' or 'featMatrix'
        
        sizeWeight=35;
        
        % get features
        [B,~,W,~,condNames]=sc1_sc2_CCC('ENCODE:get_features',K,hem);
        
        % project back to task-space
        load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC12_%dcluster_%s',K,hem),'SNN.mat'));
        %         W=bestF;
        L=W'*W;
        I=diag(diag(sqrt(L))); % diag not sum
        X=W/I;
        
        % get cluster colours
        cmap=load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC12_%dcluster_%s',K,hem),'colourMap.txt'));
        cmap=cmap/255;
        
        switch toPlot,
            case 'taskList_all'
                % Make word lists for word map
                DD = W;WORD=condNames;
                DD(DD<0)=0;
                DD=DD./max(DD(:));
                numClusters=size(DD,2);
                numFeat=size(DD,1);
                for j=1:numClusters,
                    subplot(2,5,j);
                    set(gca,'Xticklabel',[],'Yticklabel',[])
                    title(sprintf('Network %d',j),'Color',cmap(j,2:4))
                    set(gca,'FontSize',18);
                    for i=1:numFeat
                        if (DD(i,j)>0.25)
                            siz=ceil(DD(i,j)*sizeWeight);
                            text(unifrnd(0,1,1),unifrnd(0,1,1),WORD{i},'FontSize',siz,'Color',cmap(j,2:4));
                        end;
                    end;
                end;
            case 'taskMatrix'
                imagesc_rectangle(X);
                caxis([0 1]);
                t=set(gca,'Ytick',[1:length(condNames)]','YTickLabel',condNames');
                t.Color='white';
                colorbar
        end
    case 'ENCODE:scatterplot'
        mapType=varargin{1};
        type=varargin{2};
        toPlot=varargin{3}; % 1 is [4,10] - left & right hand tasks etc
        
        CAT.markersize=12;
        CAT.labelsize=18;
        sizeWeight=50;
        %         vararginoptions({varargin{3:end}},{'CAT'});
        
        % get features
        [B,F,W,C,condNames,FeatureNames]=sc1_sc2_functionalAtlas('ENCODE:get_features',mapType);
        
        % project back to task-space
        load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'SNN.mat'));
        %         W=bestF;
        L=W'*W;
        I=diag(diag(sqrt(L))); % diag not sum
        X=W/I;
        
        % load in colourMap
        cmap=load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'colourMap.txt'));
        cmap=cmap(:,2:4)/255;
        
        switch type,
            case 'features'
                net1=C(:,toPlot(1));
                net2=C(:,toPlot(2));
                
                for i=1:size(C,1),
                    if net1(i)>0 & net2(i)<0,
                        colourIdx{i,:}=cmap(toPlot(1),:);
                        siz{i,1}=ceil(net1(i)*sizeWeight);
                    elseif net2(i)>0 & net1(i)<0,
                        colourIdx{i,:}=cmap(toPlot(2),:);
                        siz{i,1}=ceil(net2(i)*sizeWeight);
                    elseif net1(i)>0 & net2(i)>0,
                        colourIdx{i,:}=[0 0 0];
                        siz{i,1}=ceil((net1(i)+net2(i))/2*sizeWeight);
                    elseif net1(i)<0 & net2(i)<0,
                        colourIdx{i,:}=[.7 .7 .7];
                        siz{i,1}=ceil((abs(net1(i)+net2(i)/2))*sizeWeight);
                    else
                        colourIdx{i,:}=[.7 .7 .7];
                        siz{i,1}=ceil((abs(net1(i)+net2(i)/2))*sizeWeight);
                    end
                end
                XY=C;
                names=FeatureNames;
            case 'tasks'
                net1=B.featIdx(toPlot(1),1);
                net2=B.featIdx(toPlot(2),1);
                % assign tasks to features (to colour-code)
                for i=1:size(X,1),
                    if F(i,net1)>0 & F(i,net2)==0, % assign to network1
                        colourIdx{i,:}=cmap(toPlot(1),:);
                    elseif F(i,net2)>0 & F(i,net1)==0, % assign to network2
                        colourIdx{i,:}=cmap(toPlot(2),:);
                    elseif F(i,net1)>0 & F(i,net2)>0,
                        colourIdx{i,:}=[0 0 0]; % tasks that load onto both features
                    elseif F(i,net1)==0 & F(i,net2)==0,
                        colourIdx{i,:}=[.7 .7 .7]; % tasks that don't load onto features - grey out
                    else
                        colourIdx{i,:}=[.7 .7 .7];
                    end
                end
                XY=X;
                names=condNames;
        end
        CAT.markercolor=colourIdx;
        CAT.markerfill=colourIdx;
        CAT.labelcolor=colourIdx;
        CAT.labelsize=siz;
        
        scatterplot(XY(:,toPlot(1)),XY(:,toPlot(2)),'label',names,'intercept',0,'draworig','CAT',CAT);
        xlabel(sprintf('Network%d',toPlot(1)));ylabel(sprintf('Network%d',toPlot(2)));a
        
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
        for e=1:2
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
        
    case 'CONNECTIVITY:runModel'
        sn=varargin{1}; % subjNum
        model=varargin{2}; % {'162_tessellation_hem'}
        data=varargin{3}; % 'Cerebellum_grey'
        signal=varargin{4}; % {'fitted'}: other options are 'residual', 'totalTS'
        method=varargin{5}; % {'ridgeFixed'}, other options: nonNegative etc
        trainMode=varargin{6}; % {'crossed'} or {'uncrossed}
        study=varargin{7}; % building on which study ? 1 or 2 [1]
        
        T=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        RR= [];
        SS=[];
        lambdaL1 = {[0]}; % cell array and for each model there are a set of lambdas. e.g. {[10 100 1000],[10 100]}
        lambdaL2 = {[0]};
        meanSub=1;      % Mean Pattern Subtract?
        incInstr=0;     % Include instruction ?
        vararginoptions(varargin(8:end),{'lambdaL1','lambdaL2','meanSub','incInstr'});
        
        for m=1:length(model),
            
            % how many lambdas are we dealing with ?
            if  isempty(lambdaL1{m}) || ~isempty(lambdaL2{m})
                numLambda=length(lambdaL2{m});
                lambdaL1=repmat(0,1,numLambda); 
            elseif isempty(lambdaL1{m}) || ~isempty(lambdaL1{m})
                numLambda=length(lambdaL1{m});
                lambdaL2=repmat(0,1,numLambda); 
            elseif lambdaL1{m} == lambdaL2{m}
                numLambda=length(lambdaL1{m});
            else
                numLambda=1;
            end
            
            subset = T.StudyNum==study; % Indices of the experiments / conditions that we want to use    
            subset = [subset;subset]; % Make the subset out onto both session
            
            switch signal{m},
                case 'fitted'
                    [X,Y,S]=sc1_sc2_CCC('STRUCTURE:meanBeta',model{m},data,'incInstr',incInstr);
                case 'residual'
                    [X,Y,S]=sc1_sc2_CCC('');
                case 'totalTS'
                    [X,Y,S]=sc1_sc2_CCC('');
            end
            
            trainXindx=[find(subset & S.sess==1);find(subset & S.sess==2)];
            switch (trainMode{m})
                case 'crossed'
                    trainYindx=[find(subset & S.sess==2);find(subset & S.sess==1)];
                case 'uncrossed'
                    trainYindx=[find(subset & S.sess==1);find(subset & S.sess==2)];
            end;
            
            for s=1:length(sn)
                if(meanSub)
                    for sess=1:2
                        indx=find(subset & S.sess==sess);
                        X{s}(indx,:)=bsxfun(@minus,X{s}(indx,:),mean(X{s}(indx,:)));
                        Y{s}(indx,:)=bsxfun(@minus,Y{s}(indx,:),mean(Y{s}(indx,:)));
                    end;
                end;
                fprintf('%d\n',sn(s));
                
                xx = X{s}(trainXindx,:);
                yy = Y{s}(trainYindx,:);
                
                % loop over lambda (if there are multiple per model)
                for l=1:numLambda,
                    
12                    l1=lambdaL1{m}(l); 
                    l2=lambdaL2{m}(l);
                    
                    R.SN = sn(s);
                    [W,R.fR2m,R.fRm] = sc1_encode_fit(yy,xx,method{m},'lambda',[l1 l2]);
                    R.W={W};
                    R.lambda = [l1 l2];
                    R.model=model{m};
                    R.signal=signal{m};
                    R.trainMode=trainMode{m};
                    R.study=study(m);
                    R.type   = {method{m}};
                    RR=addstruct(RR,R);
                    fprintf('lambda %d done ..\n',l)
                end;
               SS=addstruct(SS,RR);  
            end
        end
        varargout={RR};
    case 'CONNECTIVITY:evaluate' % Evaluates predictive performance of connectivity model
        % M: is the model structure with the connectivity weights in it (M.W)
        % 'subset': what data should the evaluation be based upon?
        % 'splitby': by which variable should the evaluation be split?
        % 'meanSub': Mean pattern subtraction before evaluation?
        model=varargin{1}; % {'162_tessellation_hem','yeo_17'}
        data=varargin{2};  % 'Cerebellum_grey'
        signal=varargin{3}; % {'fitted','residual'}
        method=varargin{4}; % {'ridgeFixed','nonNegative'}
        trainMode=varargin{5}; % {'crossed','uncrossed'}
        study=varargin{6}; % evaluating on which study, 1 or 2 ? [1,2]
        % example:
        % sc1_sc2_CCC('CONNECTIVITY:evaluate',{'162_tessellation_hem'},{'Cerebellum_grey'},{'fitted'},{'ridgeFixed'},{'crossed'})

        % if the connectivity weights are built on sc1, then evaluate on
        % sc2
        if study==1,
            studyModel=2; % connectivity weights: which study ?
        elseif  study==2,
            studyModel=1;
        end
          
        T = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        incInstr=0; 
        meanSub = 0; % Mean pattern subtraction before prediction?
        lambdaL1={[0]}; 
        lambdaL2={[0]}; 
        vararginoptions(varargin(7:end),{'lambdaL1','lambdaL2','subset','splitby','meanSub','incInstr'});
        
        % load in the connectivity weights (X): option to load in multiple
        % models (each can have different methods, trainings etc)
        M=sc1_sc2_CCC('CONNECTIVITY:runModel',returnSubjs,model,data,signal,method,trainMode,studyModel,'lambdaL1',lambdaL1,'lambdaL2',lambdaL2);
        
        % Get all the mean betas and prepare the evaulation data (Y)
        [X,Y,S]=sc1_sc2_CCC('STRUCTURE:meanBeta',model,data,incInstr);
        
        % what are we evaluating ? 
        subset = T.StudyNum==study;
        splitby = T.condNum;
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
            s=find(goodsubj==M.SN(m));
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
                RR = addstruct(RR,R);
            end;
        end;
        varargout={RR,Y{s}(testBindx,:),predY};
        
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

