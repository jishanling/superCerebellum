function varargout=sc1_sc2_CCC(what,varargin)
% Matlab script to do connectivity analyses for sc1 and sc2 experiments
% Does time series extraction
%==========================================================================
% % (1) Directories
baseDir          = '/Users/maedbhking/Documents/projects/Cerebellum_Cognition';
% baseDir           = '/Volumes/MotorControl/data/super_cerebellum_new/';
baseDir_other          = '/Volumes/Seagate Backup Plus Drive';

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
        
        % get medial wall
        C=caret_load(fullfile(studyDir{1},caretDir,'fsaverage_sym',hemi,[hem '.cerebral_cortex.paint'])); % freesurfer
        medialIdx=C.index(C.data==6); % 6 is the medial wall
        X_C(:,medialIdx,:)=nan;
        
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
        study=varargin{2}; % 1 or 2 or [1,2]
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
    case 'ACTIVITY:reliability]'
        glm=varargin{1};
        type=varargin{2}; % 'cerebellum' or '162_tessellation_hem' 'basalGanglia'
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
        save(fullfile(studyDir{study},regDir,'glm4',sprintf('patternReliability_%s.mat',type)),'CORR')
    case 'PLOT:reliability'
        type=varargin{1}; % '162_tessellation_hem'
        
        % load relability
        load(fullfile(studyDir{2},regDir,'glm4',sprintf('patternReliability_%s.mat',type)));
        
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
        % Labelling
        set(gca,'FontSize',12)
        ylabel('R')
        set(gca,'XTickLabel',{},'xtick',[])
        
    case 'PREDICTIONS:taskModel'
        sn=varargin{1}; % returnSubjs
        study=varargin{2}; % 1 or 2
        partition=varargin{3}; % session or run
        hemi=varargin{4}; % 'lh' or 'rh'
        
        lambdas=[.01:.1:.5];
        subjs=length(sn);
        l=1;
        
        % load X
        [Xx,~,numConds]=sc1_sc2_functionalAtlas('ACTIVITY:make_model',study,'yes');
        
        % loop over subjects
        Ys=[];
        for s=1:subjs,
            
            encodeSubjDir = fullfile(studyDir{study},encodeDir,'glm4',subj_name{sn(s)}); % set directory
            
            % load Y (per subj)
            load(fullfile(encodeSubjDir,sprintf('Y_info_glm4_cortex_%s.mat',hemi)));
            volIndx=Y.nonZeroInd;
            Y=rmfield(Y,'nonZeroInd');
            Yp=getrow(Y,Y.cond~=0);
            
            switch partition,
                case 'run'
                    part=Yp.run;
                    block=run;
                case 'session'
                    part=Yp.sess;
                    block=[1,2]; % session
            end
            
            % normalise (either by run or by session)
            N = (numConds)*numel(run);
            B = indicatorMatrix('identity',part);
            R  = eye(N)-B*pinv(B);
            X = R*Xx;            % Subtract block mean (from X)
            X=bsxfun(@rdivide,X,sqrt(sum(X.*X)/(size(X,1)-numel(block))));
            Yact = R*Yp.data;    % Subtract block mean (from Y)
            Yact=bsxfun(@rdivide,Yact,sqrt(sum(Yact.*Yact)/(size(Yact,1)-numel(block))));
            
            % run encoding model (with ridge regress)
            [M.R2_vox,M.R_vox]=encode_crossval(Yact,X,part,'ridgeFixed','lambda',lambdas(l));
            
            fprintf('subj%d:model done ...\n',sn(s))
            
            M.SN=sn(s);
            M.idx=volIndx;
            Ys=addstruct(Ys,M);
            clear Yp
        end
        
        % save out results
        save(fullfile(studyDir{study},encodeDir,'glm4',sprintf('encode_taskModel_cortex_%s_%s.mat',hemi,partition)),'Ys','-v7.3');
    case 'PREDICTIONS:datasets' % get predictions across datasets
        stat=varargin{1}; % 'R' or 'R2' ?
        partition=varargin{2}; % cv across 'run' or 'session' ?
        hemi=varargin{3}; % 'lh' or 'rh'
        
        F=dload(fullfile(baseDir,'motorFeats.txt')); % load in motor features
        
        YN=[];
        for i=1:2,
            numCond=length(F.condNum(F.studyNum==i));
            load(fullfile(studyDir{i},encodeDir,'glm4',sprintf('encode_taskModel_cortex_%s_%s.mat',hemi,partition)))
            Ys.study=repmat(i,size(Ys.SN,1),1);
            YN=addstruct(YN,Ys);
        end
        X=indicatorMatrix('identity_p',YN.SN);
        
        % which stat: R or R2 ?
        switch stat,
            case 'R'
                data=YN.R_vox;
            case 'R2'
                data=YN.R2_vox;
        end
        
        % get average of models across datasets
        RVox=pinv(X)*data;
        RVox_sc1=nanmean(nanmean(data(YN.study==1,:)));
        RVox_sc2=nanmean(nanmean(data(YN.study==2,:)));
        
        varargout={RVox,RVox_sc1,RVox_sc2,YN.idx(1,:)};
    case 'PREDICTIONS:mapToCortex' % visualise predictions for motor feats
        stat=varargin{1}; % 'R' or 'R2' ?
        partition=varargin{2}; % cv across 'run' or 'session' ?
        hemi=varargin{3}; % 'lh' or 'rh'
        
        if strcmp(hemi,'lh'),
            hemName='LeftHem';
        else
            hemName='RightHem';
        end
        
        [RVox,RVox_sc1,RVox_sc2,idx]=sc1_sc2_CCC('PREDICTIONS:datasets',stat,partition,hemi);
        
        % what is within task set reliability ?
        fprintf('set A reliability is %2.3f \n',RVox_sc1);
        fprintf('set B reliability is %2.3f \n',RVox_sc2);
        
        % get medial wall
        C=caret_load(fullfile(studyDir{1},caretDir,'fsaverage_sym',hemName,[hemi '.cerebral_cortex.paint'])); % freesurfer
        medialIdx=C.index(C.data==6); % 6 is the medial wall
        RVox(:,medialIdx,:)=nan;
        
        % make metric file
        S=caret_struct('metric','data',nanmean(RVox,1)','column_name',{'average_taskSets'});
        
        % save out metric
        caret_save(fullfile(studyDir{2},caretDir,'fsaverage_sym',hemName,'glm4','taskModel.metric'),S);
        
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
        
    case 'run_all'
        model=varargin{1}; % 'desikan','yeo','162_tessellation','frontal_regions' etc
        method=varargin{2}; % 'winnerTakeAll','L2'
        getLambda=varargin{3}; % 'no' (winnerTakeAll) or 'yes' (L2)
        
        sc1_sc2_CCC('SIGNAL:meanTS',returnSubjs,'sc1',model)
        sc1_sc2_CCC('SIGNAL:meanTS',returnSubjs,'sc2',model)
        sc1_sc2_CCC('SIGNAL:meanBeta',model)
        sc1_sc2_CCC('STRUCTURE:meanBeta',model,'Cerebellum_grey',0)
        %         sc1_sc2_CCC('CONNECT:run','weights',method,model,'no')
        %         sc1_sc2_CCC('CONNECT:run','eval',method,model,getLambda)
    case 'SIGNAL:meanTS'
        % sc1_connectivity('TS_get_meants',[2 3 6 8 9 10 12 17:22],'sc2',4,'162_tessellation_hem');
        % sc1_connectivity('TS_get_meants',[2:22],'sc1',4,'162_tessellation_hem');
        sn=varargin{1}; % subjNum
        exper=varargin{2}; % Experiment string 'sc1' or 'sc2'
        type=varargin{3}; % 'Cerebellum_grey', '162_tessellation_hem','yeo' etc
        % ex.sc1_sc2_CCC('SIGNAL:meanTS',returnSubjs,'sc1')
        
        B = [];
        glm=4;
        %         type='Cerebellum_grey';
        glmDir =fullfile(baseDir,exper,sprintf('/GLM_firstlevel_%d',glm));
        subjs=length(sn);
        
        for s=1:subjs,
            glmDirSubj=fullfile(glmDir, subj_name{sn(s)});
            load(fullfile(glmDirSubj,'SPM_light.mat'));
            %            load(fullfile(glmDirSubj,'SPM.mat'));
            T=load(fullfile(glmDirSubj,'SPM_info.mat'));
            
            % load data
            tic;
            
            load(fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',type))); % 'regions' are defined in 'ROI_define'
            
            SPM=spmj_move_rawdata(SPM,fullfile(baseDir_other,'sc1','imaging_data',subj_name{sn(s)}));
            
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
            load(fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',type))); % 'regions' are defined in 'ROI_define'
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
        type    = varargin{1};  % 'Cerebellum_grey' or 'yeo' or 'desikan'
        exper   = {'sc1','sc2'};
        sn      = returnSubjs;
        glm     = 4;
        
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
        sn=returnSubjs;  % Take only good subjects
        glm=4; % glmNum
        type='yeo'; % 'cortical_lobes','whole_brain','yeo','desikan','cerebellum','yeo_cerebellum'
        session_specific = 1;
        expDir = studyDir{1};
        
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
    case 'STRUCTURE:newROI'  % make new models out of existing models
        newRegion=varargin{1}; % 'preCentral','superFrontal','parsTriang','visualCortex'
        
        
        Tcort=load(fullfile(baseDir,'sc1',regDir,'data','162_reorder.mat'));
        Tcort_good=getrow(Tcort,Tcort.good==1);
        Tcort_good=rmfield(Tcort_good,'newIndx');
        Tcort_good.newIndx=[1:length(Tcort_good.good)]';
        
        switch newRegion,
            case 'preCentral'
                tmp=getrow(Tcort_good,Tcort_good.desikan==25);
                regions=tmp.newIndx;
            case 'superFrontal'
                tmp=getrow(Tcort_good,Tcort_good.desikan==29);
                regions=tmp.newIndx;
            case 'parsTriang'
                tmp=getrow(Tcort_good,Tcort_good.desikan==21);
                regions=tmp.newIndx;
            case 'lingual'
                tmp=getrow(Tcort_good,Tcort_good.desikan==14);
                regions=tmp.newIndx;
            case 'inferiorTemporal'
                tmp=getrow(Tcort_good,Tcort_good.desikan==10);
                regions=tmp.newIndx;
            case 'preCuneus'
                tmp=getrow(Tcort_good,Tcort_good.desikan==26);
                regions=tmp.newIndx;
            case 'postCentral'
                tmp=getrow(Tcort_good,Tcort_good.desikan==23);
                regions=tmp.newIndx;
            case 'superiorTemporal'
                tmp=getrow(Tcort_good,Tcort_good.desikan==31);
                regions=tmp.newIndx;
            case 'superiorParietal'
                tmp=getrow(Tcort_good,Tcort_good.desikan==30);
                regions=tmp.newIndx;
            case 'yeo1'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==1);
                regions=tmp.newIndx;
            case 'yeo2'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==2);
                regions=tmp.newIndx;
            case 'yeo3'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==3);
                regions=tmp.newIndx;
            case 'yeo4'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==4);
                regions=tmp.newIndx;
            case 'yeo5'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==5);
                regions=tmp.newIndx;
            case 'yeo6'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==6);
                regions=tmp.newIndx;
            case 'yeo7'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==7);
                regions=tmp.newIndx;
            case 'yeo8'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==8);
                regions=tmp.newIndx;
            case 'yeo9'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==9);
                regions=tmp.newIndx;
            case 'yeo10'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==10);
                regions=tmp.newIndx;
            case 'yeo11'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==11);
                regions=tmp.newIndx;
            case 'yeo12'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==12);
                regions=tmp.newIndx;
            case 'yeo13'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==13);
                regions=tmp.newIndx;
            case 'yeo14'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==14);
                regions=tmp.newIndx;
            case 'yeo15'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==15);
                regions=tmp.newIndx;
            case 'yeo16'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==16);
                regions=tmp.newIndx;
            case 'yeo17'
                tmp=getrow(Tcort_good,Tcort_good.yeo17==17);
                regions=tmp.newIndx;
                
        end
        
        model='162_tessellation_hem';
        data='Cerebellum_grey';
        incInstr=0;
        
        [X,Y,S]=sc1_sc2_CCC('STRUCTURE:meanBeta',model,data,incInstr);
        
        [tasks features]=size(X{1});
        
        for s=1:length(X),
            N{s}=zeros(tasks,features);
            N{s}(:,regions')=X{s}(:,regions');
        end
        
        varargout={N,Y,S};
        
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
        method=varargin{2}; % ''L1','L2','winnerTakeAll'
        getLambda=varargin{3}; % 'yes' automatically generate lambdas or 'no' do nothing (winner take all)
        
        condType='unique';
        
        vararginoptions(varargin(4:end),{'newRegion','model','condType'}); % 'model': 'yeo','desikan','162_tessellation_hem','frontalRegions';
        % 'newRegion':'motorCortex','frontalCortex'
        
        if strcmp(getLambda,'yes')
            lambdas=sc1_sc2_CCC('getLambda',method);
        else
            lambdas=[0];
        end
        
        switch step,
            case 'weights'
                for i=1:2,
                    for l=1:length(lambdas),
                        if exist('newRegion'),
                            sc1_sc2_CCC('CONNECT:weights',returnSubjs,method,i,'l2',lambdas(l),'newRegion',newRegion)
                        else
                            sc1_sc2_CCC('CONNECT:weights',returnSubjs,method,i,'l2',lambdas(l),'model',model)
                        end
                    end
                end
            case 'eval'
                for i=1:2,
                    for l=1:length(lambdas),
                        if exist('newRegion')
                            sc1_sc2_CCC('CONNECT:eval',returnSubjs,method,i,'l2',lambdas(l),'lambdaIdx',l,'newRegion',newRegion,'condType',condType)
                        else
                            sc1_sc2_CCC('CONNECT:eval',returnSubjs,method,i,'l2',lambdas(l),'lambdaIdx',l,'model',model,'condType',condType)
                        end
                    end
                end
        end
    case 'getLambda'
        method=varargin{1}; % 'L1','L2','Lasso'
        
        switch method,
            case 'L1'
            case 'L2'
                n=40;
                lambdas=exp(linspace(log(10),log(100000),n)); % 20 log spaced values from 10 to 1000 (or 1000 to 10000)
                %                 lambdas=lambdas(15:end);
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
        data='Cerebellum_grey';
        trainMode={'crossed'};
        
        vararginoptions(varargin(4:end),{'l1','l2','signal','model','data','trainMode','newRegion','model'});
        
        % 'weights' fileName
        weightDir=fullfile(studyDir{2},connDir,'glm4','weights',method); dircheck(weightDir);
        
        subset = T.StudyNum==study; % Indices of the experiments / conditions that we want to use
        subset = [subset;subset]; % Make the subset out onto both session
        
        switch signal,
            case 'fitted'
                if exist('newRegion'),
                    [X,Y,S]=sc1_sc2_CCC('STRUCTURE:newROI',newRegion);
                    modelName=fullfile(weightDir,sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d.mat',newRegion,signal,method,study,round(l1),round(l2)));
                    model=newRegion;
                else
                    [X,Y,S]=sc1_sc2_CCC('STRUCTURE:meanBeta',model,data,incInstr);
                    modelName=fullfile(weightDir,sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d.mat',model,signal,method,study,round(l1),round(l2)));
                end
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
                if strcmp(method,'winnerTakeAll'),
                    trainYindx=[find(subset & S.sess==1);find(subset & S.sess==2)];
                    yy = Y{s}(trainYindx,:);
                    Yy=sum(yy.*yy,1);
                    Xx=sum(xx.*xx,1);
                    R.C={(xx'*yy)./sqrt(bsxfun(@times,Yy,Xx'))};
                end
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
        data='Cerebellum_grey';
        trainMode={'crossed'};
        condType='unique';
        incInstr=0;
        meanSub = 0; % Mean pattern subtraction before prediction?
        voxels=0;
        
        T = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        vararginoptions(varargin(4:end),{'l1','l2','subset','meanSub','incInstr','data','trainMode','signal','condType','voxels','lambdaRange','lambdaIdx','newRegion','model'}); % 'condType' is either 'unique or 'all'
        
        % Get all the mean betas and prepare the evaluation data (Y)
        if exist('newRegion'),
            model=newRegion;
            [X,Y,S]=sc1_sc2_CCC('STRUCTURE:newROI',model);
        else
            [X,Y,S]=sc1_sc2_CCC('STRUCTURE:meanBeta',model,data,incInstr);
        end
        
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
        weightFile=fullfile(weightDir,sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d.mat',model,signal,method,studyModel,round(l1),round(l2)));
        M=load(weightFile);
        
        % eval file name - output and check to see if connectivity eval exists [load if it does exist]
        if voxels==1,
            evalFile=sprintf('%s_%s_%s_sc%d_%s_voxels.mat',model,signal,method,study,condType);
            RR=[];
        else
            evalFile=sprintf('%s_%s_%s_sc%d_%s.mat',model,signal,method,study,condType);
            if exist(fullfile(evalDir,evalFile)),
                RR=load(fullfile(evalDir,evalFile));
            else
                RR=[];
            end
        end
        
        for s=1:length(sn) % subjects
            for t=1:length(trainMode),
                
                Mm=getrow(M,M.SN==sn(s) & strcmp(M.trainMode,trainMode{t}));
                numVox=size(Mm.W{1},2);
                numReg=size(Mm.W{1},1);
                
                %                 indx = sum(abs(Y{s}(:,:)))>0 & ~isnan(sum(Y{s}(:,:))) & ~isnan(sum(Mm.W{1})); % introduces diff number of voxels --
                %                 introduces complications later on when reslicing into suit
                %                 space
                indx=1:numVox;
                
                for sp=1:length(splits); % task conditions
                    if(meanSub)
                        for sess=1:2, % sessions
                            indx=find(S.subset & S.sess==sess & S.splitby==splits(sp));
                            X{s}(indx,:)=bsxfun(@minus,X{s}(indx,:),nanmean(X{s}(indx,:)));
                            Y{s}(indx,:)=bsxfun(@minus,Y{s}(indx,:),nanmean(Y{s}(indx,:)));
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
                            SSP   = nansum(predY(:,indx).^2);             % Sum of square of predictions
                            SSY   = nansum(Y{s}(testBindx,indx).^2);               % Sum of squares of data
                            SSCp  = nansum(predY(:,indx).*predYnc(:,indx)); % Covariance of Predictions
                            SSCy  = nansum(Y{s}(testAindx,indx).*Y{s}(testBindx,indx));  % Covariance of Y's
                            SSCn  = nansum(predYnc(:,indx).*Y{s}(testBindx,indx));   % Covariance of non-cross prediction and data
                            SSCc  = nansum(predY(:,indx).*Y{s}(testBindx,indx));   % Covariance of cross prediction and data
                            R.Rcv   = {SSCc ./ sqrt(SSY.*SSP)}; % Double-Crossvalidated predictive correlation
                            R.Rnc   = {SSCn ./ sqrt(SSY.*SSP)}; % Not double-crossvalidated predictive correlation
                            R.Ry    = {SSCy ./ SSY};            % Reliability of data: noise ceiling
                            R.Rp    = {SSCp ./ SSP};            % Reliability of prediction
                            threshold=nanmean(max(abs(Mm.W{1})));
                            R.maxRegv={max(abs(Mm.W{1}))};
                            R.sumRegv={nansum(abs(Mm.W{1}))};
                            R.relMaxRegv={R.maxRegv{1}./R.sumRegv{1}};
                            R.numRegv={nansum(abs(Mm.W{1})>threshold)};
                        case 0
                            SSP   = nansum(nansum(predY(:,indx).^2));             % Sum of square of predictions
                            SSY   = nansum(nansum(Y{s}(testBindx,indx).^2));               % Sum of squares of data
                            SSCp  = nansum(nansum(predY(:,indx).*predYnc(:,indx))); % Covariance of Predictions
                            SSCy  = nansum(nansum(Y{s}(testAindx,indx).*Y{s}(testBindx,indx)));  % Covariance of Y's
                            SSCn  = nansum(nansum(predYnc(:,indx).*Y{s}(testBindx,indx)));   % Covariance of non-cross prediction and data
                            SSCc  = nansum(nansum(predY(:,indx).*Y{s}(testBindx,indx)));   % Covariance of cross prediction and data
                            R.Rcv   = SSCc ./ sqrt(SSY.*SSP); % Double-Crossvalidated predictive correlation
                            R.Rnc   = SSCn ./ sqrt(SSY.*SSP); % Not double-crossvalidated predictive correlation
                            R.Ry    = SSCy ./ SSY;            % Reliability of data: noise ceiling
                            R.Rp    = SSCp ./ SSP;            % Reliability of prediction
                            R.maxRegm=nanmean(max(abs(Mm.W{1})));
                            R.sumRegm=nanmean(nansum(abs(Mm.W{1})));
                            R.relMaxRegm=R.maxRegm./R.sumRegm; % mean(maxReg)
                    end
                    R.SN    = Mm.SN;
                    R.lambda = Mm.lambda(:,:);
                    if voxels==0,
                        R.modelIdx=lambdaIdx;
                    end
                    R.method = Mm.method;
                    R.trainMode = Mm.trainMode;
                    R.model = Mm.model;
                    R.splitby = splits(sp);
                    R.split = {condType};
                    RR = addstruct(RR,R);
                end;
                fprintf('pred done for %s for l1-%2.2f, l2-%2.2f \n',subj_name{sn(s)},l1,l2)
                clear R
            end
        end
        save(fullfile(evalDir,evalFile),'-struct','RR')
    case 'CONNECT:checks' % does some checks: crossed vs uncrossed, unique vs all vs shared etc for one method
        method=varargin{1}; % {'L2'}, {'nonNegative','L2} etc
        study=varargin{2}; % evaluating on which study, 1 or 2 ? [1,2]
        signal=varargin{3}; % 'fitted'
        model=varargin{4}; % '162_tessellation_hem'
        trainMode=varargin{5}; % 'crossed'
        condType=varargin{6}; % {'unique','all','shared'} or {'unique'} etc
        condIdx=varargin{7}; % 1, 2, or 3 [1-unique,2-all,3-shared]
        whatToPlot=varargin{8}; % 'condType','crossUcross','taskConds','RcvRnc' etc
        
        vararginoptions({varargin{9:end}},{'CAT'}); % option if plotting individual map analysis.
        
        RR=[];
        S=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % figure out L1,L2, or lasso
        if strcmp(method,'L2'),
            M1=2;
        elseif strcmp(method,'L1'),
            M1=1;
        elseif strcmp(method,'winnerTakeAll'),
            M1=2;
        end
        for c=1:length(condType),
            T=load(fullfile(studyDir{2},connDir,'glm4','eval',sprintf('%s_%s_%s_sc%d_%s.mat',model,signal,method,study,condType{c})));
            T.condNum=repmat(c,size(T.SN,1),1);
            RR=addstruct(RR,T);
        end
        switch whatToPlot,
            case 'condType'
                % 'unique' versus 'all'
                lineplot(RR.lambda(:,M1),RR.Rcv,'split',RR.split,'leg','auto',...
                    'subset',strcmp(RR.method,method) & strcmp(RR.trainMode,trainMode))
                ylabel('R')
                xlabel('lambda')
                title(sprintf('%s-%s',signal,method))
                set(gca,'YLim',[.1 0.4],'FontSize',14,'ytick',[.1,.15,.2,.25,.3, .35,.4]);
            case 'crossUncross'
                % 'crossed' versus 'uncrossed'
                lineplot(RR.lambda(:,M1),RR.Rcv,'split',RR.trainMode,'leg','auto',...
                    'subset',strcmp(RR.method,method) & strcmp(RR.split,condType{condIdx}))
                ylabel('R')
                xlabel('lambda')
                title(sprintf('%s-%s',signal,method))
                set(gca,'YLim',[.1 0.3],'FontSize',14,'ytick',[.1,.15,.2,.25,.3]);
            case 'taskConds'
                
                % split by taskConds (unique)
                condNames=S.condNames(S.StudyNum==study & S.overlap==0);
                figure(1)
                barplot(RR.splitby,RR.Rcv,'subset',...
                    strcmp(RR.method,method) & strcmp(RR.split,'unique') & strcmp(RR.trainMode,trainMode),'leg','auto')
                ylabel('R')
                title(sprintf('%s-%s-unique',signal,method))
                set(gca,'YLim',[0 0.4],'FontSize',8,'XTickLabel',condNames);
                
                % split by taskConds (all)
                condNames=S.condNames(S.StudyNum==study);
                figure(2)
                barplot(RR.splitby,RR.Rcv,'subset',...
                    strcmp(RR.method,method) & strcmp(RR.split,'all') & strcmp(RR.trainMode,trainMode),'leg','auto')
                ylabel('R')
                title(sprintf('%s-%s-all',signal,method))
                set(gca,'YLim',[0 0.55],'FontSize',8,'XTickLabel',condNames);
                
                % split by taskConds (shared)
                condNames=S.condNames(S.StudyNum==study);
                figure(3)
                barplot(RR.splitby,RR.Rcv,'subset',...
                    strcmp(RR.method,method) & strcmp(RR.split,'shared') & strcmp(RR.trainMode,trainMode),'leg','auto')
                ylabel('R')
                title(sprintf('%s-%s-shared',signal,method))
                set(gca,'YLim',[0 0.55],'FontSize',8,'XTickLabel',condNames);
            case 'RcvRnc'
                % split Rcv versus Rnc
                figure(1)
                lineplot(RR.lambda(:,M1),RR.Rcv,'subset',...
                    strcmp(RR.method,method) & strcmp(RR.split,condType{condIdx}) & strcmp(RR.trainMode,trainMode),'leg','auto')
                hold on
                CAT.linecolor={'k'};
                CAT.linestyle={'-'};
                CAT.errorcolor={'k'};
                lineplot(RR.lambda(:,M1),RR.Rnc,'CAT',CAT,'subset',...
                    strcmp(RR.method,method)& strcmp(RR.split,condType{condIdx}) & strcmp(RR.trainMode,trainMode),'leg','auto')
                ylabel('R')
                xlabel('lambda')
                title(sprintf('%s-%s',signal,method))
                set(gca,'YLim',[.1 0.3],'FontSize',14,'ytick',[.1,.15,.2,.25,.3]);
        end
    case 'CONNECT:plotLambda'
        method=varargin{1}; % {'L2'}, {'nonNegative','L2} etc
        study=varargin{2}; % evaluating on which study, 1 or 2 ?
        signal=varargin{3}; % 'fitted'
        model=varargin{4}; % '162_tessellation_hem'
        trainMode=varargin{5}; % 'crossed'
        condType=varargin{6}; % 'unique'
        
        vararginoptions({varargin{7:end}},{'lambdaRange','lambdaIdx','CAT'}); % option if plotting individual map analysis.
        
        T=[];
        for m=1:length(method),
            Tt=load(fullfile(studyDir{2},connDir,'glm4','eval',sprintf('%s_%s_%s_sc%d_%s.mat',model,signal,method{m},study,condType)));
            if strcmp(method{m},'winnerTakeAll'),
                % figure out model numbers
                modelIdx=length(unique(T.modelIdx))+1;
                Tt.modelIdx=repmat(modelIdx,size(Tt.SN,1),1);
            end
            Tt.study=repmat(study,length(Tt.SN),1);
            T=addstruct(T,Tt);
        end
        
        CAT.errorwidth=1;
        xyplot(T.modelIdx,[T.Rcv],T.modelIdx,...
            'subset', strcmp(T.split,condType) & strcmp(T.trainMode,trainMode),'split',T.method,'CAT',CAT);
        hold on
        clear CAT
        
        % plot noise ceiling
        %         CAT.errorwidth=0;
        CAT.markertype='none';
        CAT.linestyle={'--'};
        CAT.linewidth={2};
        CAT.errorcolor='k';
        CAT.linecolor='k';
        lineplot(T.modelIdx,T.Ry,'subset', ...
            strcmp(T.split,condType) & strcmp(T.trainMode,trainMode) ,'CAT',CAT)
    case 'CONNECT:plotTask'
        method=varargin{1}; % {'L2'}, {'nonNegative','L2} etc
        study=varargin{2}; % evaluating on which study, 1 or 2 ?
        signal=varargin{3}; % 'fitted'
        model=varargin{4}; % '162_tessellation_hem'
        trainMode=varargin{5}; % 'crossed'
        condType=varargin{6}; % 'unique'
        
        vararginoptions({varargin{7:end}},{'lambdaRange','lambdaIdx','CAT'}); % option if plotting individual map analysis.
        
        T=[];
        Tt=load(fullfile(studyDir{2},connDir,'glm4','eval',sprintf('%s_%s_%s_sc%d_%s.mat',model,signal,method,study,condType)));
        Tt.study=repmat(study,length(Tt.SN),1);
        T=addstruct(T,Tt);
        
        CAT.errorwidth=1;
        
        % get split
        split=getrow(T,strcmp(T.split,condType) & strcmp(T.trainMode,trainMode) & T.study==study);
        taskSplits=unique(split.splitby);
        
        % plot splits
        barplot(T.splitby,T.Rcv,'subset',strcmp(T.split,condType) & strcmp(T.trainMode,trainMode),'CAT',CAT);
        
        varargout={taskSplits};
    case 'CONNECT:relMax'
        method=varargin{1}; % {'L2'}, {'nonNegative','L2} etc
        study=varargin{2}; % evaluating on which study, 1 or 2 ? [1,2]
        signal=varargin{3}; % 'fitted'
        model=varargin{4}; % '162_tessellation_hem'
        trainMode=varargin{5}; % 'crossed'
        condType=varargin{6}; % 'unique'
        
        vararginoptions({varargin{7:end}},{'CAT'}); % option if plotting individual map analysis.
        
        T=[];
        for i=study,
            for m=1:length(method),
                Tt=load(fullfile(studyDir{2},connDir,'glm4','eval',sprintf('%s_%s_%s_sc%d_%s.mat',model,signal,method{m},i,condType)));
                Tt.study=repmat(i,length(Tt.SN),1);
                T=addstruct(T,Tt);
            end
        end
        
        xyplot(T.relMaxRegm,T.Rcv,T.modelIdx,...
            'subset', strcmp(T.split,condType) & strcmp(T.trainMode,trainMode),'CAT',CAT);
    case 'CONNECT:bestLambda'
        method=varargin{1}; % 'L2' etc
        study=varargin{2}; % evaluating on which study, 1 or 2 ? [1,2]
        type=varargin{3}; % 'group' or 'indiv'
        
        signal='fitted';
        trainMode='crossed';
        condType='unique';
        lambdaRange=[0 10000];
        
        vararginoptions(varargin(4:end),{'newRegion','model','lambdaRange'}); % 'condType' is either 'unique or 'all'
        
        if exist('newRegion'),
            T=load(fullfile(studyDir{2},connDir,'glm4','eval',sprintf('%s_%s_%s_sc%d_%s.mat',newRegion,signal,method,study,condType)));
        else
            T=load(fullfile(studyDir{2},connDir,'glm4','eval',sprintf('%s_%s_%s_sc%d_%s.mat',model,signal,method,study,condType)));
        end
        
        switch type,
            case 'group'
                lambdaIdx=find(T.lambda(T.lambda>lambdaRange(1) & T.lambda<lambdaRange(2)));
                T.lambdaIdx=zeros(size(T.lambda,1),1);
                T.lambdaIdx(lambdaIdx)=1;
                A=getrow(T,strcmp(T.split,condType) & strcmp(T.trainMode,trainMode) & T.lambdaIdx==1);
                idx=find(A.Rcv==max(A.Rcv));
                bestLambda=A.lambda(idx,:);
                fprintf('bestLambda is %2.2f \n',bestLambda(2))
            case 'indiv'
                for s=1:length(returnSubjs),
                    lambdaIdx=find(T.lambda(T.lambda>lambdaRange(1) & T.lambda<lambdaRange(2)));
                    T.lambdaIdx=zeros(size(T.lambda,1),1);
                    T.lambdaIdx(lambdaIdx)=1;
                    A=getrow(T,T.SN==returnSubjs(s) & strcmp(T.split,condType) & strcmp(T.trainMode,trainMode) & T.lambdaIdx==1);
                    idx=find(A.Rcv==max(A.Rcv));
                    bestLambda(s,:)=A.lambda(idx(1),:);
                    fprintf('bestLambda for s%d is %2.2f \n',returnSubjs(s),bestLambda(2))
                end
        end
        
        varargout={bestLambda};
        
    case 'MAP:cortex'           % Map of where projections come from
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
    case 'MAP:surface_cortex'
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
    case 'MAP:surface_cereb'    % Plot any data on the cerebellar flatmap
        data = varargin{1};     % Data to plot
        volIndx = varargin{2};  % Indices into the volume (mask)
        V = varargin{3};        % Cerebellar suit volume
        
        cmap = [];
        cscale = [];
        type = 'label';     % func / label
        xlims = [-100  100];
        ylims = [-85 85];
        vararginoptions(varargin(4:end),{'cmap','type','cscale','xlims','ylims'});
        
        % map features on group
        V.dat=zeros([V.dim(1) V.dim(2) V.dim(3)]);
        V.dat(volIndx)=data;
        
        switch (type)
            case 'label'
                stats = 'mode';
                if (isempty(cmap))
                    cmap = colorcube(max(data));
                end;
                outputtype = 'paint';
            case {'func','rgb'}
                stats = 'nanmean';
                if (isempty(cmap))
                    cmap = hot;
                end;
                outputtype = 'metric';
        end;
        
        % map the data and display on flatmap
        data=suit_map2surf(V,'space','SUIT','stats',stats);
        suit_plotflatmap(data,'type',type,'cmap',cmap,...
            'cscale',cscale,'xlims',xlims,'ylims',ylims);
        
        
        P=caret_struct(outputtype,'data',data);
        varargout={P};
    case 'MAP:suit_reslice'     % Maps certain stats from an individual to SUIT and then to flat map
        data    = varargin{1};
        sn      = varargin{2};
        
        stats='nanmean'; % default
        vararginoptions(varargin(3:end),{'stats'}); % 'condType' is either 'unique or 'all'
        
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
        
        % resample the weights/labels into SUIT space
        for r=1:size(data,1)
            Vout=Vmask;
            Vout.dat = zeros(Vout.dim);
            Vout.dat(linIn3)=data(r,:);
            %             Vout.dat=data(r,:);
            Vout.dt = [64 0];
            Vout.pinfo = [1 0 0]';
            DataSUIT(r,:)=spm_sample_vol(Vout,i2,j2,k2,0);
            V.dat = zeros(V.dim);
            Vres(r)=V;
            Vres(r).dat(linIn1)=DataSUIT(r,:);  % Offset by one to account for 1 being medial wall
            Vres(r).fname=sprintf('data_%2.2d.nii',r);
            Vres(r).pinfo=[1 0 0]';
        end;
        
        % Now map to surface-based representation
        D = suit_map2surf(Vres,'stats',stats);
        varargout={D,Vres,linIn1,DataSUIT};
    case 'MAP:suit_reslice_MK'  % Maps weights from an individual to SUIT and then to flat map
        data=varargin{1};
        sn=varargin{2};
        
        
        load(fullfile(studyDir{1},regDir,'data',subj_name{sn},'regions_Cerebellum_grey.mat'));
        
        % Determine the voxels we want to resample in SUIT space
        V=spm_vol(fullfile(studyDir{1},suitDir,'anatomicals','cerebellarGreySUIT.nii'));
        X=spm_read_vols(V);
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
        job.interp=0;
        D=suit_reslice_dartel(job);
        % delete temporary files
        deleteFiles=dir(fullfile(studyDir{2},'connectivity','eval',subj_name{sn},'*temp*'));
        for b=1:length(deleteFiles),
            delete(char(fullfile(studyDir{2},'connectivity','eval',subj_name{sn},deleteFiles(b).name)));
        end
        
        V.dat = zeros(V.dim);
        Vres=V;
        %         Vres.dat(linIn1)=permute(D,[2 3 4 1]);
        Vres.dat=permute(D,[2 3 4 1]);
        Vres.fname=sprintf('data_%2.2d.nii',1);
        Vres.pinfo=[1 0 0]';
        
        % Now map to surface-based representation
        D = suit_map2surf(Vres,'stats',stats);
        varargout={D,Vres};
    case 'MAP:predictions'      % Plots weights or evaluations onto flatmap for subj or group
        method=varargin{1}; % 'L2' or 'nonNegative'
        study=varargin{2}; % evaluating on which study, 1 or 2 ? [1,2]
        type=varargin{3}; % get 'bestLambda' for 'group' or 'indiv'
        signal=varargin{4}; % 'fitted'
        data=varargin{5}; % 'Cerebellum_grey'
        trainMode=varargin{6}; % 'crossed'
        condType=varargin{7}; % 'unique'
        splitby=varargin{8}; % 'no'
        sn=varargin{9}; % returnSubjs
        toPlot=varargin{10}; % 'relMax' or 'R' or 'numReg'
        
        vararginoptions(varargin(11:end),{'lambdaRange','newRegion','model'}); % 'condType' is either 'unique or 'all'
        
        % load in eval results
        evalDir=fullfile(studyDir{2},connDir,'glm4','eval');dircheck(evalDir);
        
        % figure out if 'newRegion' like 'motorCortex' or existing 'model'
        % like '162_tessellation_hem'
        if exist('newRegion');
            whichModel='newRegion';
            model=newRegion;
        else
            whichModel='model';
        end
        
        % get best lambda
        bestLambda=sc1_sc2_CCC('CONNECT:bestLambda',method,study,type,whichModel,model);
        
        % do eval for voxels
        sc1_sc2_CCC('CONNECT:eval',returnSubjs,'L2',study,'l2',bestLambda(2),'voxels',1,whichModel,model)
        T=load(fullfile(evalDir,sprintf('%s_%s_%s_sc%d_%s_voxels.mat',model,signal,method,study,condType)));
        
        % make connectivity file
        connectMetName=fullfile(studyDir{2},caretDir,'suit_flat','glm4',sprintf('connectivity-predict-%s.metric',model));
        
        N=[];
        for s=1:length(sn),
            
            A=getrow(T,T.SN==sn(s) & strcmp(T.trainMode,trainMode) & strcmp(T.split,condType));
            
            taskSplits=unique(T.splitby);
            % figure out tasksplits
            switch splitby
                case 'no'
                    for t=1:length(taskSplits),
                        switch toPlot,
                            case 'R'
                                data_task(t,:)=A.Rcv{t};
                            case 'relMax'
                                data_task(t,:)=A.relMaxRegv{t};
                            case 'numReg'
                                data_task(t,:)=A.numRegv{t};
                                keyboard
                        end
                    end
                    data={nanmean(data_task,1)};
                    taskSplits=1;
                    clear data_task
                case 'yes'
                    switch toPlot,
                        case 'R'
                            data=A.Rcv;
                        case 'relMax'
                            data=A.relMaxRegv;
                        case 'numReg'
                            data=A.numRegv;
                    end
            end
            
            for sp=1:length(taskSplits),
                % get flatmap coordinates
                [D,Vres]=sc1_sc2_CCC('MAP:suit_reslice',data{sp},sn(s));
                %                 [D,Vres]=sc1_sc2_CCC('MAP:cerebellum',data{sp},sn(s));
                Nn.data=D';
                Nn.SN=sn(s);
                Nn.split=taskSplits(sp);
                N=addstruct(N,Nn);
            end
            fprintf('subj%d done\n',sn(s))
        end
        
        % plot to surface
        %         suit_plotflatmap(nanmean(N.data,1)')
        
        if exist(connectMetName),
            S=caret_load(connectMetName);
            for i=1:size(S.data,2),
                dataAll(:,i)=S.data(:,i);
                dataName{i}=S.column_name{i};
            end
            dataAll(:,i+1)=nanmean(N.data,1)';
            dataName{i+1}=toPlot;
        else
            dataAll(:,1)=nanmean(N.data,1)';
            dataName{1}=toPlot;
        end
        M=caret_struct('metric','data',dataAll,'column_name',dataName);
        caret_save(connectMetName,M);
    case 'MAP:weights'
        method=varargin{1}; % 'L2' or 'nonNegative' or 'winnerTakeAll'
        study=varargin{2}; % evaluating on which study, 1 or 2 ? [1,2]
        type=varargin{3}; % get bestLambda for 'group' or 'indiv'
        signal=varargin{4}; % 'fitted'
        data=varargin{5}; % 'Cerebellum_grey'
        trainMode=varargin{6}; % 'crossed'
        
        sn=returnSubjs;
        lambdaRange=[0 10000];
        l1=0;
        
        vararginoptions(varargin(7:end),{'lambdaRange','newRegion','model','sn'}); % 'condType' is either 'unique or 'all'
        
        % figure out if 'newRegion' like 'motorCortex' or existing 'model'
        % like '162_tessellation_hem'
        if exist('newRegion');
            whichModel='newRegion';
            model=newRegion;
        else
            whichModel='model';
        end
        
        weightDir=fullfile(studyDir{2},connDir,'glm4','weights',method);dircheck(weightDir);
        
        % get best lambda
        bestLambda=sc1_sc2_CCC('CONNECT:bestLambda',method,study,type,whichModel,model);
        
        % do eval for voxels
        weightFile=fullfile(weightDir,sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d.mat',model,signal,method,study,round(l1),round(bestLambda(2))));
        M=load(weightFile);
        
        % make connectivity file
        connectMetName=fullfile(studyDir{2},caretDir,'suit_flat','glm4',sprintf('connectivity-weights-%s.metric',model));
        
        N=[];
        for s=1:length(sn),
            
            A=getrow(M,M.SN==sn(s) & strcmp(M.trainMode,trainMode));
            
            % get flatmap coordinates
            data=nanmean(A.W{1},1);
            [D,Vres]=sc1_sc2_CCC('MAP:suit_reslice',data,sn(s));
            Nn.data=D';
            Nn.SN=sn(s);
            N=addstruct(N,Nn);
            fprintf('subj%d done\n',sn(s))
        end
        
        dataAll(:,1)=nanmean(N.data,1)';
        dataName{1}='weights';
        M=caret_struct('metric','data',dataAll,'column_name',dataName);
        caret_save(connectMetName,M);
    case 'MAP:localRand'
        type = 'voxelwise'; % 'voxelwise','searchlight', or 'searchlightvoxel'
        compare={'connectivity-WTA-yeo.nii','connectivity-WTA-buckner.nii'};
        
        split = [1 1 1 2 2 2];
        radius = 10;
        vararginoptions(varargin,{'radius','compare','type','split'});
        
        evalDir=fullfile(studyDir{2},connDir,'glm4','eval');
        numMaps = length(compare);
        load(fullfile(studyDir{2},encodeDir,'glm4','cereb_avrgDataStruct.mat'));  % Just to get V and volIndx
        clear T;
        
        for i=1:numMaps
            Vi=spm_vol(fullfile(evalDir,compare{i}));
            [i1,j1,k1]=ind2sub(V.dim,volIndx');
            [x,y,z]=spmj_affine_transform(i1,j1,k1,V.mat);
            [i2,j2,k2]=spmj_affine_transform(x,y,z,inv(Vi.mat));
            c(:,i)=spm_sample_vol(Vi,i2,j2,k2,0);
        end;
        ar=[];
        [x,y,z]=ind2sub(V.dim,volIndx);
        [x,y,z]=spmj_affine_transform(x,y,z,V.mat);
        xyz=[x;y;z];
        D=surfing_eucldist(xyz,xyz);        % Eucledian distance
        pair = [];
        for i=1:numMaps-1
            for j=i+1:numMaps
                ar(:,end+1)=RandIndexLocal(type,c(:,i),c(:,j),D,radius);
                if (~isempty(split))
                    if (split(i)==split(j) && split(i)==1)
                        pair(end+1)=1;
                    elseif (split(i)==split(j) && split(i)==2)
                        pair(end+1)=2;
                    elseif (split(i)~=split(j))
                        pair(end+1)=3;
                    else
                        pair(end+1)=0;
                    end;
                end;
            end;
        end;
        numPlots = length(unique(pair));
        for i=1:numPlots
            subplot(1,numPlots,i);
            sc1_sc2_CCC('MAP:surface_cereb',mean(ar(:,pair==i),2),volIndx,V,'type','func');
        end;
        Output.type= type;
        Output.compare= compare;
        Output.split= split;
        Output.ar= ar;
        Output.pair= pair;
        varargout = {Output};
    case 'MAP:WTA'
        method=varargin{1}; % 'L2' or 'winnerTakeAll
        signal=varargin{2}; % 'fitted'
        model=varargin{3}; % '162_tessellation' or 'yeo' or 'desikan'
        trainMode=varargin{4}; % 'crossed'
        study=varargin{5}; % 1 or 2 or [1,2]
        
        bestLambda=[0 0];
        stats='mode';
        
        vararginoptions(varargin(6:end),{'stats'}); % 'condType' is either 'unique or 'all'
        
        weightDir=fullfile(studyDir{2},connDir,'glm4','weights',method);dircheck(weightDir);
        evalDir=fullfile(studyDir{2},connDir,'glm4','eval');
        
        % load in colourMap for each model
        cmap=load(fullfile(studyDir{2},caretDir,'suit_flat','glm4',sprintf('%s.txt',model)));
        
        N=[];
        for ss=study,
            
            % load in the connectivity weights (X): option to load in multiple
            % models (each can have different methods, trainings etc)
            T=load(fullfile(weightDir,sprintf('%s_%s_%s_sc%d_l1_%d_l2_%d.mat',model,signal,method,ss,round(bestLambda(1)),round(bestLambda(2)))));
            
            subjs=unique(T.SN);
            for s=1:length(subjs),
                A=getrow(T,strcmp(T.trainMode,trainMode) & T.SN==subjs(s));
                
                % get winner label and weights
                [weights_winner,labels]=max(A.C{1}(1:end,:),[],1); % 1 is medial wall
                
                % get flatmap coordinates for WTA and weights
                [labels_suit,Vres,volIndx,DataSUIT]=sc1_sc2_CCC('MAP:suit_reslice',labels,subjs(s),'stats',stats);
                [weights_winner_suit]=sc1_sc2_CCC('MAP:suit_reslice',weights_winner,subjs(s),'stats',stats);
                
                Nn.labelSURF=labels_suit';
                Nn.weights_winner=weights_winner_suit';
                Nn.SN=subjs(s);
                Nn.study=ss;
                Nn.labelVOL=DataSUIT;
                N=addstruct(N,Nn);
            end
        end
        weightWTA=fullfile(studyDir{2},caretDir,'suit_flat','glm4',sprintf('connectivity-WTA-%s.metric',model));
        labelWTA=fullfile(studyDir{2},caretDir,'suit_flat','glm4',sprintf('connectivity-WTA-%s.paint',model));
        areaWTA=fullfile(studyDir{2},caretDir,'suit_flat','glm4',sprintf('connectivity-WTA-%s.areacolor',model));
        
        % save out metric, paint, and area, and volume
        W=caret_struct('metric','data',nanmean(N.weights_winner,1)');
        L=caret_struct('paint','data',mode(N.labelSURF,1)','column_name',{'column_01'});
        A=caret_struct('area','data',cmap(:,2:4),'column_name',L.paintnames,...
            'column_color_mapping',repmat([-5 5],L.num_paintnames),'paintnames',L.paintnames);
        
        Yy=zeros(1,Vres.dim(1)*Vres.dim(2)*Vres.dim(3));
        Yy(1,volIndx)=mode(N.labelVOL,1)';
        Yy=reshape(Yy,[Vres.dim(1),Vres.dim(2),Vres.dim(3)]);
        Yy(Yy==0)=NaN;
        Vres.fname=fullfile(evalDir,sprintf('connectivity-WTA-%s.nii',model));
        spm_write_vol(Vres,Yy);
        
        caret_save(weightWTA,W);
        caret_save(labelWTA,L);
        caret_save(areaWTA,A);
        
    case 'AXES:checks'
        method=varargin{1}; % 'L2' or 'L1' or 'winnerTakeAll'
        study=varargin{2}; % studyNum (1 or 2)
        whatToPlot=varargin{3}; % 'condType','crossUncross','taskConds','RcvRnc'
        
        signal='fitted';
        model='162_tessellation_hem';
        trainMode='crossed';
        condType={'unique','all','shared'};
        condIdx=1; % default is 'unique'
        
        sc1_sc2_CCC('CONNECT:checks',method,study,signal,model,trainMode,condType,condIdx,whatToPlot)
    case 'AXES:plotLambda'
        study=varargin{1}; % 1 or 2 or [1,2]
        model=varargin{2}; % {'162_tessellation_hem','desikan','yeo'}
        
        signal='fitted';
        trainMode='crossed';
        condType='unique';
        
        vararginoptions({varargin{3:end}},{'model','newRegion','condType'}); % fracture further into sessions [either 1 or 2]
        
        CAT.errorwidth=.5;
        CAT.markertype='none';
        CAT.linewidth=3;
        CAT.linestyle={'-','-','-','-','-','-'};
        CAT.linewidth={2, 2, 2, 2, 2, 2};
        errorcolor={'g','r','b','c','m','y','m'};
        linecolor={'g','r','b','c','m','y','m'};
        
        % loop over models for L2
        for m=1:length(model),
            CAT.errorcolor=errorcolor{m};
            CAT.linecolor=linecolor{m};
            sc1_sc2_CCC('CONNECT:plotLambda',{'L2'},study,signal,model{m},trainMode,condType,'CAT',CAT)
            hold on
        end
        
        % Labelling
        %         set(gca,'YLim',[.15 0.25],'FontSize',12,'XLim',[lambdaRange(1)*-2 lambdaRange(2)]);
        set(gca,'Ylim',[.09 0.3],'FontSize',12,'XLim',[0 42],'xtick',{})
        ylabel('predictive accuracy (R)')
        xlabel('L2 norm')
        %         set(gcf,'units','centimeters','position',[5,5,9,12])
        %         legend(plotName,'Location','NorthWest')
    case 'AXES:plotRelMax'
        study=varargin{1}; % 1 or 2 or [1,2]
        model=varargin{2}; % {'162_tessellation_hem','desikan','yeo'}
        
        signal='fitted';
        trainMode='crossed';
        condType='unique';
        
        CAT.errorwidth=.5;
        CAT.markertype='none';
        CAT.linewidth=3;
        CAT.linestyle={'-','-','-','-','-','-'};
        CAT.linewidth={2, 2, 2, 2, 2, 2};
        errorcolor={'g','r','b','c','c','y','m'};
        linecolor={'g','r','b','c','c','y','m'};
        
        % loop over models for L2
        for m=1:length(model),
            CAT.errorcolor=errorcolor{m};
            CAT.linecolor=linecolor{m};
            sc1_sc2_CCC('CONNECT:relMax',{'L2'},study,signal,model{m},trainMode,condType,'CAT',CAT)
            hold on
        end
        
        % Labelling
        %         set(gca,'YLim',[.15 0.25],'FontSize',12,'XLim',[lambdaRange(1)*-2 lambdaRange(2)]);
        set(gca,'Ylim',[.1 0.25],'FontSize',16)
        ylabel('predictive accuracy (R)')
        xlabel('relMax')
    case 'AXES:plotTask' % plot predictive accuracy per task
        study=varargin{1}; % 1 or 2
        
        signal='fitted';
        trainMode='crossed';
        condType='unique';
        
        vararginoptions({varargin{2:end}},{'model','newRegion','condType'}); % fracture further into sessions [either 1 or 2]
        
        CAT.errorwidth=.5;
        CAT.markertype='none';
        CAT.linewidth=3;
        CAT.linestyle={'-','-','-','-','-','-'};
        CAT.linewidth={2, 2, 2, 2, 2, 2};
        CAT.errorcolor={'k'};
        CAT.linecolor={'k'};
        
        if exist('model'),
            taskSplits=sc1_sc2_CCC('CONNECT:plotTask','L2',study,signal,model,trainMode,condType,'CAT',CAT);
        else
            taskSplits=sc1_sc2_CCC('CONNECT:plotTask','L2',study,signal,newRegion,trainMode,condType,'CAT',CAT);
        end
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        A=getrow(D,D.StudyNum==study);
        xNames=A.condNames(taskSplits);
        
        % Labelling
        %         set(gca,'YLim',[.15 0.25],'FontSize',12,'XLim',[lambdaRange(1)*-2 lambdaRange(2)]);
        set(gca,'Ylim',[0 0.4],'FontSize',10,'xticklabels',xNames)
        ylabel('Predictive accuracy (R)')
        xlabel('Tasks')
        %         set(gcf,'units','centimeters','position',[5,5,9,12])
        %         legend(plotName,'Location','NorthWest')
    case 'AXES:flatmap-predictions'
        method=varargin{1}; % 'L2' or 'winnerTakeAll'
        study=varargin{2}; % 1 or 2
        toPlot=varargin{3}; % 'relMax' or 'R' or 'numReg'
        
        % model - '162_tessellation_hem','yeo','desikan' etc (Regions that
        % are defined in 'ROI:define'
        % newRegion - 'frontalCortex','motorCortex' - regions created from
        % '162_tessellation_hem'
        
        vararginoptions({varargin{4:end}},{'model','newRegion','condType'}); % fracture further into sessions [either 1 or 2]
        
        % set up parameters
        signal='fitted';
        data='Cerebellum_grey';
        trainMode='crossed';
        condType='unique';
        splitby='no'; % get pred maps for each task ?
        sn=returnSubjs; % returnSubjs
        type='group'; % get 'bestLambda' for 'group' or 'indiv'
        bestLambda=[];
        lambdaRange=[];
        
        % either plot the weights (WTA) or the predictions (usually from L2
        % models)
        if exist('newRegion'),
            sc1_sc2_CCC('MAP:predictions',method,study,type,signal,data,...
                trainMode,condType,splitby,sn,toPlot,'lambdaRange',[0 10000],'newRegion',newRegion)
        else
            sc1_sc2_CCC('MAP:predictions',method,study,type,signal,data,...
                trainMode,condType,splitby,sn,toPlot,'lambdaRange',[0 10000],'model',model)
        end
    case 'AXES:flatmap-weights'
        method=varargin{1}; % 'L2' or 'winnerTakeAll'
        study=varargin{2}; % 1 or 2
        
        type='group';
        signal='fitted';
        data='Cerebellum_grey';
        trainMode='crossed';
        
        vararginoptions({varargin{3:end}},{'model','newRegion','type','signal','data','trainMode'}); % fracture further into sessions [either 1 or 2]
        
        if exist('newRegion'),
            sc1_sc2_CCC('MAP:weights',method,study,type,signal,data,trainMode,'newRegion',newRegion)
        else
            sc1_sc2_CCC('MAP:weights',method,study,type,signal,data,trainMode,'model',model)
        end
    case 'AXES:flatmap-WTA'
        study=varargin{1}; % 1 or 2 or [1,2]
        
        % model - '162_tessellation_hem','yeo','desikan' etc (Regions that
        % are defined in 'ROI:define'
        % newRegion - 'frontalCortex','motorCortex' - regions created from
        % '162_tessellation_hem'
        
        vararginoptions({varargin{2:end}},{'model','newRegion'}); % fracture further into sessions [either 1 or 2]
        
        if exist('newRegion'),
            model=newRegion;
        end
        
        % set up parameters
        method='winnerTakeAll';
        signal='fitted';
        trainMode='crossed';
        
        % either plot the weights (WTA) or the predictions (usually from L2
        % models)
        sc1_sc2_CCC('MAP:WTA',method,signal,model,trainMode,study)
        
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

