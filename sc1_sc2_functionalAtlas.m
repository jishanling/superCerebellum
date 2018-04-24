function varargout=sc1_sc2_functionalAtlas(what,varargin)

% Directories
baseDir          = '/Users/maedbhking/Documents/Cerebellum_Cognition';
baseDir            = '/Volumes/MotorControl/data/super_cerebellum_new';
% baseDir          = '/Users/jdiedrichsen/Data/super_cerebellum_new';

studyDir{1}     =fullfile(baseDir,'sc1');
studyDir{2}     =fullfile(baseDir,'sc2');
suitDir         ='/suit';
caretDir        ='/surfaceCaret';
regDir          ='/RegionOfInterest/';
encodeDir       ='/encoding';

run = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'};

subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11',...
    's12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24',...
    's25','s26','s27','s28','s29','s30','s31'};

returnSubjs=[2,3,4,6,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31];

switch what
    
    case 'ACTIVITY:get_data'   % get tasksConds (averaged across runs + sessions for each dataset)
        subjs=length(returnSubjs);
        
        F=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % load in allSubjs data struct from sc1 & sc2
        H=[];
        for study=1:2,
            load(fullfile(studyDir{study},encodeDir,'glm4','cereb_avrgDataStruct.mat'));
            T.studyNum=repmat(study,size(T.SN,1),1);
            H=addstruct(H,T);
            clear T
        end
        
        % get session average for each study separately
        for s=1:subjs,
            idx=1;
            for study=1:2,
                numConds=length(unique(H.cond(H.studyNum==study)));
                for c=1:numConds, % get average across sessions
                    indx = H.cond==c & H.studyNum==study & H.SN==(returnSubjs(s));
                    avrgData(idx,:)=nanmean(H.data(indx,:),1);
                    idx=idx+1;
                    clear indx
                end
                fprintf('subj%d averaged sessions for study%d \n',returnSubjs(s),study)
            end
            % zero to NaN
            avrgData(avrgData==0)=NaN;
            
            % subtract condition avrg baseline (across 61-1 conditions)
            all=[1:size(avrgData,1)];
            for c=1:size(avrgData,1),
                X=nanmean(avrgData(all(all~=c),:),1);
                data(c,:,s)=avrgData(c,:)-X;
            end
            
            clear avrgData
            fprintf('subj%d new baseline \n',returnSubjs(s))
        end
        varargout={data,V,volIndx};
    case 'ACTIVITY:make_model' % make X matrix (feature models)
        study=varargin{1}; % 1, 2, or [1,2]
        type=varargin{2};  % 'yes' or 'no' to modelling each run separately ?
        model=varargin{3}; % 'task' or 'motor' model ?
        
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
        
        switch model,
            case 'task'
                % make model output struct
                for m=1:length(featNames),
                    M.modelIdx{m,1}=m;
                    M.modelNames{m,1}=featNames{m};
                end
                M.modelIdx{m+1}=[1:numConds];
                M.modelNames{m+1,1}='task';
            case 'motor'
                % make model output struct
                for m=1:length(featNames),
                    M.modelIdx{m,1}=m;
                    M.modelNames{m,1}=featNames{m};
                end
        end
        
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
        
        % add rest to the end (better way of doing this !)
        rest=X.x(X.idx==numConds,:);
        X.x(X.idx==numConds,:)=[];
        X=[X.x; rest];
        
        varargout={X,M,numConds,F};
    case 'ACTIVITY:patterns'   % estimate each of the task conditions (motor features are subtracted)
        sharedTasks=varargin{1}; % 'average' or 'all'
        
        subjs=length(returnSubjs);
        lambda=.01;
        
        % load in task structure file
        F=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % load in activity patterns for sc1+sc2 (averaged across sessions,
        % within dataset)
        [data,V,volIndx]=sc1_sc2_functionalAtlas('ACTIVITY:get_data','average'); % sharedTasks are averaged
        
        % get motor feature model
        [X,M,numConds]=sc1_sc2_functionalAtlas('ACTIVITY:make_model',[1,2],'no','motor'); % load in model
        
        % set up volume info
        numFeat=size(X,2)-numConds;
        Yy=zeros(numConds+numFeat,subjs,V.dim(1)*V.dim(2)*V.dim(3));
        C{1}.dim=V.dim;
        C{1}.mat=V.mat;
        
        % NaN to zero
        data(isnan(data))=0;
        
        % do ridge regression
        for s=1:subjs,
            X1=bsxfun(@minus,X,mean(X));
            X1=bsxfun(@rdivide,X1,sum(X1.^2));
            B=(X1'*X1+eye(numConds+numFeat)*lambda)\(X1'*data(:,:,s));
            % make volume
            Yy(:,s,volIndx)=B;
            clear X1 B
            fprintf('ridge regress done for subj%d done \n',returnSubjs(s))
        end;
        clear data
        
        % if numSubjs > 1 get avg
        Yy=permute(Yy,[2 1 3]);
        indices=nanmean(Yy,1);
        indices=reshape(indices,[size(indices,2),size(indices,3)]);
        
        % map vol2surf
        indices=reshape(indices,[size(indices,1) V.dim(1),V.dim(2),V.dim(3)]);
        for i=1:size(indices,1),
            data=reshape(indices(i,:,:,:),[C{1}.dim]);
            C{i}.dat=data;
        end
        M=caret_suit_map2surf(C,'space','SUIT','stats','nanmean','column_names',M.modelNames);  % MK created caret_suit_map2surf to allow for output to be used as input to caret_save
        
        switch sharedTasks,
            case 'all'
                % do nothing
            case 'average'     
        end
        
        % save out metric
        caret_save(fullfile(studyDir{2},caretDir,'suit_flat','glm4','unCorr_taskConds.metric'),M);
    case 'ACTIVITY:reliability'
        glm=varargin{1};
        type=varargin{2}; % 'cerebellum' or 'cortex' 'basalGanglia'
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        load(fullfile(studyDir{2},encodeDir,sprintf('glm%d',glm),sprintf('allVox_sc1_sc2_sess_%s.mat',type)));
        Y=Yy;clear Yy;
        
        numSubj=length(Y);
        
        S=[];
        for subj=1:numSubj,
            R=[];
            idx=1;
            for study=1:2,
                D1=getrow(D,D.StudyNum==study);
                sharedConds=D1.condNumUni.*D1.overlap;
                if study==1,
                    condNames=D1.condNames(find(sharedConds));
                end
                % sharedConds=sharedConds(randperm(numel(sharedConds{2}))); % Shuffle
                cN=condNum{study}-1;  % Important: In the allVox file, instruction is still included!
                pN=partNum{study};    % Partition Numner
                sN=(pN>8)+1;          % Sessions Number
                for se=1:2,
                    X1=indicatorMatrix('identity_p',cN.*(sN==se));  % This one is the matrix that related trials-> condition numbers
                    X2=indicatorMatrix('identity_p',sharedConds); % THis goes from condNum to shared condNumUni
                    Yf(:,:,idx,subj)=pinv(X1*X2)*Y{subj}{study};
                    Yf(:,:,idx,subj)=bsxfun(@minus,Yf(:,:,idx,subj),mean(Yf(:,:,idx,subj)));
                    idx=idx+1;
                end;
            end;
            for c=1:size(Yf,1),
                CORR(c,:,:,subj)=interSubj_corr(Yf(c,:,:,subj));
                T.SN      = returnSubjs(subj);
                T.within1 = CORR(c,1,2,subj);
                T.within2 = CORR(c,3,4,subj);
                T.across  = mean(mean(CORR(c,1:2,3:4,subj)));
                T.condNum = c;
                T.condNames={condNames{c}};
                R=addstruct(R,T);
                clear T
            end;
            S=addstruct(S,R);
            clear R
        end;
        save(fullfile(studyDir{2},regDir,'glm4','patternReliability.mat'),'S','CORR')
    case 'PLOT:reliabilityA'
        % load relability
        load(fullfile(studyDir{2},regDir,'glm4','patternReliability_cerebellum.mat'));
        
        %         figure();lineplot(S.condNum,[S.within1,S.within2,S.across],'leg',{'within1','within2','across'})
        A=tapply(S,{'SN'},{'across'},{'within1'},{'within2'});
        
        % within & between-subj reliability
        myboxplot([],[A.within1 A.within2 A.across],'style_twoblock','plotall',1);
        drawline(0,'dir','horz');
        ttest(sqrt(A.within1.*A.within2),A.across,2,'paired');
        
    case 'PREDICTIONS:motorFeats'   % predict motor feature model + task model (cross-validated with fixed lambda)
        sn=varargin{1};
        study=varargin{2};
        
        lambdas=[.01:.1:.5];
        subjs=length(sn);
        l=1;
        
        % load X
        [X,M,numConds]=sc1_sc2_functionalAtlas('ACTIVITY:make_model',study,'yes','task');
        
        % loop over subjects
        Ys=[];
        for s=1:subjs,
            
            encodeSubjDir = fullfile(studyDir{study},encodeDir,'glm4',subj_name{sn(s)}); % set directory
            
            % load Y (per subj)
            load(fullfile(encodeSubjDir,'Y_info_glm4_grey_nan.mat'));
            Yp=getrow(Y,Y.cond~=0);
            
            % loop over models
            for m=1:length(M.modelIdx),
                
                % modify X
                Xm=X(:,M.modelIdx{m});
                
                % normalise
                N = (numConds)*numel(run);
                B = indicatorMatrix('identity',Yp.run);
                R  = eye(N)-B*pinv(B);
                Xm = R*Xm;            % Subtract block mean (from X)
                Xm=bsxfun(@rdivide,Xm,sqrt(sum(Xm.*Xm)/(size(Xm,1)-numel(run))));
                Yact = R*Yp.data;    % Subtract block mean (from Y)
                Yact=bsxfun(@rdivide,Yact,sqrt(sum(Yact.*Yact)/(size(Yact,1)-numel(run))));
                
                % run encoding model (+ ridge regress)
                [M.R2(m,1),~,M.R2_vox(m,:)]=encode_crossval(Yact,Xm,Yp.run,'ridgeFixed','lambda',lambdas(l));
                
                M.modelNum(m,1)=m; % modelNum
                M.lambdas(m,1)=l;
                fprintf('subj%d:model %s ...\n',sn(s),M.modelNames{m})
                clear Xm
            end
            M.SN=repmat(sn(s),m,1);
            M.idx=repmat(Yp.nonZeroInd(1,:),m,1);
            Ys=addstruct(Ys,M);
            clear Yp
        end
        
        % save out results
        save(fullfile(studyDir{study},encodeDir,'glm4','encode_models.mat'),'Ys','-v7.3');
    case 'PREDICTIONS:average_SC12' % get the average predictions across datasets
        
        F=dload(fullfile(baseDir,'motorFeats.txt')); % load in motor features
        
        YN=[];
        for i=1:2,
            numCond=length(F.condNum(F.studyNum==i));
            load(fullfile(studyDir{i},encodeDir,'glm4','encode_models.mat'))
            tmp=getrow(Ys,Ys.modelNum>numCond);
            tmp.studyNum=repmat(i,size(tmp.SN,1),1);
            YN=addstruct(YN,tmp);
        end
        
        % get average of models across datasets
        avrgStudy=tapply(YN,{'R2_vox','idx','SN','modelNames','modelIdx','modelNum','lambdas'},{'studyNum'}); % tapply can't deal with multi-cols
    case 'PREDICTIONS:vol2surf'     % visualise predictions for motor feats
        study=varargin{1}; % 1 or 2 or [1,2]
        metric=varargin{2}; % 'yes' or 'no'
        
        % are we visualising studies separately or combined ?
        if length(study)>1, % studies [1,2]
            study=2;
            load(fullfile(studyDir{study},encodeDir,'glm4','encode_models_SC12.mat'))
            outName='encode_models_SC12';
        else
            load(fullfile(studyDir{study},encodeDir,'glm4','encode_models.mat'))
            outName='encode_models';
        end
        
        numModels=length(unique(Ys.modelNum));
        SN=length(unique(Ys.SN));
        
        V=spm_vol(fullfile(studyDir{1},suitDir,'anatomicals','cerebellarGreySUIT.nii'));
        C{1}.dim=V.dim;
        C{1}.mat=V.mat;
        
        Yy=zeros(SN,V.dim(1)*V.dim(2)*V.dim(3));
        
        % loop over models
        for m=1:numModels,
            
            M=getrow(Ys,Ys.modelNum==m);
            
            % make vol
            Yy(:,M.idx(1,:))=M.R2_vox;
            
            % get avrg across subjs
            indices=nanmean(Yy,1);
            
            % map vol2surf
            data=reshape(indices,[V.dim(1),V.dim(2),V.dim(3)]);
            C{m}.dat=data;
            modelNames(m)=unique(M.modelNames);
        end
        
        M=caret_suit_map2surf(C,'space','SUIT','stats','nanmean','column_names',modelNames);  % MK created caret_suit_map2surf to allow for output to be used as input to caret_save
        
        % save out metric ?
        switch metric,
            case 'yes'
                caret_save(fullfile(studyDir{study},caretDir,'suit_flat','glm4',sprintf('%s.metric',outName)),M);
            case 'no'
                % visualise motor features
                for m=1:M.num_cols,
                    figure(m)
                    title(M.column_name{m})
                    suit_plotflatmap(M.data(:,m))
                end
        end
        
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
        
        D=[];
        sn=unique(T.SN);
        numSubj = length(sn);
        RR=[];
        for f=unique(T.freq)'
            for s = 1:numSubj
                for se=1:2
                    temp = T.data(T.SN==sn(s) & T.sess==se & T.freq==f,:);
                    temp = bsxfun(@minus,temp,mean(temp));
                    D(:,:,s+(se-1)*length(sn))=temp;
                end;
            end;
            C=intersubj_corr(D);
            R.sess = [ones(1,numSubj) ones(1,numSubj)*2]';
            R.subj = [1:numSubj 1:numSubj]';
            R.subj = sn(R.subj);
            R.freq = f*ones(numSubj*2,1);
            SameSess = bsxfun(@eq,R.sess',R.sess);
            SameSubj = bsxfun(@eq,R.subj',R.subj);
            for i=1:numSubj*2;
                R.withinSubj(i,:)=C(i,SameSubj(i,:) & ~SameSess(i,:));
                R.betweenSubj(i,:)=mean(C(i,~SameSubj(i,:)));
                R.totSS(i,1) = nansum(nansum(D(:,:,i).^2));
            end;
            RR=addstruct(RR,R);
        end;
        save(fullfile(studyDir{2},'encoding','glm4','cereb_spatialCorr_freq.mat'),'-struct','RR');
        varargout={RR};
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
        xlabels={'overall','0-0.5','0.5-1','1-1.5','1.5-2','>2'};
        
        T=tapply(T,{'subj','freq'},{'withinSubj'},{'betweenSubj'},{'totSS'},'subset',ismember(T.subj,returnSubjs));
        T.freqK = T.freq>0;
        lineplot([T.freqK T.freq],[T.withinSubj T.betweenSubj],'CAT',CAT);
        ylabel('Correlations');
        xlabel('Cycles/cm');
        drawline(0,'dir','horz');
        set(gca,'XTickLabel',xlabels);
        title('Within and Between-subject correlation');
        set(gcf,'PaperPosition',[2 2 4.5 5.5]);
        wysiwyg;
        
    case 'REPRESENTATION:get_distances'
        type=varargin{1}; % 'cerebellum'
        removeMotor=varargin{2}; % 'hands','saccades','all','none'
        
        load(fullfile(studyDir{2},regDir,'glm4',sprintf('G_hat_sc1_sc2_%s.mat',type)))
        subjs=size(G_hat,3);
        numDist=size(G_hat,1);
        
        % load in condName info
        T=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % Get RDM
        for s=1:subjs,
            H=eye(numDist)-ones(numDist)/numDist; % centering matrix
            G_hat(:,:,s)=H*G_hat(:,:,s)*H'; % subtract out mean pattern
            IPM=rsa_vectorizeIPM(G_hat(:,:,s));
            con = indicatorMatrix('allpairs',[1:numDist]);
            N = rsa_squareIPM(IPM);
            D = rsa.rdm.squareRDM(diag(con*N*con'));
            fullRDM(:,:,s) = D;
        end
        % load feature matrix
        F=dload(fullfile(baseDir,'motorFeats.txt')); % load in motor features
        
        switch removeMotor,
            case 'hands'
                X   = [F.lHand./F.duration F.rHand./F.duration];
                X   = [X eye(numDist)];
            case 'saccades'
                X   = [F.saccades./F.duration];
                X   = [X eye(numDist)];
            case 'all'
                X   = [F.lHand./F.duration F.rHand./F.duration F.saccades./F.duration];
                X   = [X eye(numDist)];
            case 'none'
                X   = [eye(numDist)];
        end
        X   = bsxfun(@minus,X,mean(X));
        X   = bsxfun(@rdivide,X,sqrt(sum(X.^2)));  % Normalize to unit length vectors
        
        varargout={fullRDM,T,X};
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
    case 'REPRESENTATION:eigDecomp'
        CAT=varargin{1};
        % load in condName info
        T=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        load(fullfile(studyDir{2},regDir,'glm4',sprintf('G_hat_sc1_sc2_%s.mat','cerebellum')))
        
        %                 R=[];
        %                 for s=1:2,
        %                     SC=getrow(T,T.StudyNum==s);
        %                     idx=SC.condNum;
        %                     numDist=length(idx);
        %                     avrgG_hat=nanmean(G_hat(idx,idx),3);
        %                     H=eye(numDist)-ones(numDist)/length(numDist);
        %                     avrgG_hat=(H*avrgG_hat*H);  % center G matrix
        %
        %                     [V,L]   = eig(avrgG_hat);
        %                     [l,i]   = sort(diag(L),1,'descend');
        %                     S.study=repmat(s,length(i),1);
        %                     S.eig=real(l);
        %                     R=addstruct(R,S);
        %                 end
        %                   lineplot([1:size(R.study,1)]',R.eig,'split',R.study,'CAT',CAT);
        idx=T.condNum;
        numDist=length(idx);
        avrgG_hat=nanmean(G_hat(idx,idx),3);
        H=eye(numDist)-ones(numDist)/length(numDist);
        avrgG_hat=(H*avrgG_hat*H);  % center G matrix
        
        [V,L]   = eig(avrgG_hat);
        [l,i]   = sort(diag(L),1,'descend');
        S.eig=real(l);
        S.idx=i;
        
        lineplot(S.idx,S.eig,'CAT',CAT)
    case 'REPRESENTATION:RDM'
        type=varargin{1}; % 'all' or 'none'
        
        % load in fullRDM
        [fullRDM,T,~]=sc1_sc2_functionalAtlas('REPRESENTATION:get_distances','cerebellum',type);
        
        % Plot RDM
        reOrder=[1,2,6,7,8,9,10,11,12,13,14,17,18,22,23,3,4,5,15,16,19,20,21,24,25,26,...
            27,28,29,58,59,60,43,44,49,48,36,34,35,55,56,57,61,30,31,32,33,37,38,39,40,41,42,45,46,47,50,51,52,53,54]'; % reorder
        
        % threshold RDM
        fullRDM_thresh=reshape(fullRDM,[size(fullRDM,1)*size(fullRDM,2)],[]);
        for dd=1:size(fullRDM_thresh,1),
            [t(dd),p(dd)]=ttest(fullRDM_thresh(dd,:),[],1,'onesample');
        end
        ut=t; % unthresholded distances
        
        % FDR correct p-vals
        t(p>.00001)=0; % thresholded distances
        
        % use spm_P_FDR or hand code the multiple comparisons
        fprintf('%2.2f%% of the pairwise distances are significantly different from zero \n',(length(t(t>0))/length(t))*100);
        
        % get group
        avrgFullRDM=ssqrt(nanmean(fullRDM,3));
        numDist=size(avrgFullRDM,1);
        
        % visualise thresholded RDM (t-values)
        squareT=tril(reshape(t,[size(fullRDM,1),size(fullRDM,2)]));
        squareUT=reshape(ut,[size(fullRDM,1),size(fullRDM,2)]);
        idxUT=find(triu(squareUT));
        squareT(idxUT)=squareUT(idxUT);
        squareT(isnan(squareT))=0;
        figure();imagesc_rectangle(squareT(reOrder,reOrder),'YDir','reverse')
        caxis([0 1]);
        g=set(gca,'Ytick',[1:numDist]','YTickLabel',T.condNames(reOrder)','FontSize',7);
        g.Color='white';
        colorbar
    case 'REPRESENTATION:MDS'
        type=varargin{1}; % 'all','hands','saccades','none'
        
        % colour
        colour={[1 0 0],[0 1 0],[0 0 1],[0.3 0.3 0.3],[1 0 1],[1 1 0],[0 1 1],...
            [0.5 0 0.5],[0.8 0.8 0.8],[.07 .48 .84],[.99 .76 .21],[.11 .7 .68],...
            [.39 .74 .52],[.21 .21 .62],[0.2 0.2 0.2],[.6 .6 .6],[.3 0 .8],[.8 0 .4],...
            [0 .9 .2],[.1 .3 0],[.2 .4 0],[.63 0 .25],[0 .43 .21],[.4 0 .8]};
        
        vararginoptions({varargin{2:end}},{'CAT','colour'}); % option if doing individual map analysis
        
        % load in fullRDM
        [fullRDM,T,X]=sc1_sc2_functionalAtlas('REPRESENTATION:get_distances','cerebellum',type); % remove all motorFeats
        
        avrgFullRDM=ssqrt(nanmean(fullRDM,3));
        
        vecRDM = rsa.rdm.vectorizeRDM(avrgFullRDM);
        [Y,~] = rsa_classicalMDS(vecRDM,'mode','RDM');
        B = (X'*X+eye(size(X,2))*0.0001)\(X'*Y); % ridge regression
        Yr    = Y  - X(:,1:3)*B(1:3,:); % remove motor features
        
        clustTree = linkage(Yr,'average');
        indx = cluster(clustTree,'cutoff',1);
        
        % define cluster colour
        [V,L]   = eig(Yr*Yr');
        [l,i]   = sort(diag(L),1,'descend'); % Sort the eigenvalues
        V       = V(:,i);
        X       = bsxfun(@times,V,sqrt(l'));
        X = real(X);
        
        condIndx=[1:length(T.condNum)]';
        condNames=T.condNames;
        CAT.markercolor= {colour{indx}};
        CAT.markerfill = {colour{indx}};
        CAT.labelcolor = {colour{indx}};
        
        X1=X(condIndx,condIndx);
        figure()
        scatterplot3(X1(:,1),X1(:,2),X1(:,3),'split',condIndx,'CAT',CAT,'label',condNames);
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
    case 'PLOT:reliabilityD'
        % load relability
        load(fullfile(studyDir{2},regDir,'glm4','distanceReliability_cerebellum.mat'));
        
        % within & between-subj reliability
        myboxplot([],[T.within1 T.within2 T.across],'style_twoblock','plotall',1);
        ttest(sqrt(T.within1.*T.within2),T.across,2,'paired');
        
    case 'CONVERT:mni2suit'    % converts from mni image to SUIT (there's also a cifti2nii if image is coming from HCP)
        inputImages={'N2C_subcortex_atlas_subcortexGSR.nii',...
            'Buckner2011_7Networks_MNI152_FreeSurferConformed1mm_TightMask.nii'};
        
        inputDir=which(inputImages{2});
        cd(fileparts(inputDir))
        suit_mni2suit(inputImages{2})
    case 'CONVERT:cifti2paint' % haven't finished this yet but want to make paint files from the cifti surface files
        inputImages={'CortexSubcortex_ColeAnticevic_NetPartition_netassignments_v1_R.dlabel.nii',...
            'CortexSubcortex_ColeAnticevic_NetPartition_parcels_v1_L.dlabel.nii'};
        
        inputImage=which(inputImages{1});
        cd(fileparts(inputImage))
        
        % Read CIFTI
        cii=ft_read_cifti(inputImage);
        
        LH=cii.x1(cii.brainstructure==1); % LH is 1
        RH=cii.x1(cii.brainstructure==2); % RH is 2
    case 'CONVERT:makelob10'   % this makes lob10 map out of lob26 map
        V=spm_vol(fullfile(studyDir{2},'encoding','glm4','groupEval_lob26','map.nii'));
        Vi=spm_read_vols(V);
        
        % consolidate lob26 into lob10 (clunky but works)
        labelsOld=[0:34];
        labelsNew=[0,1,1,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,0,0,0,0,0,0];
        
        tmp=Vi(:);
        for l=1:length(labelsOld),
            tmp(tmp==labelsOld(l))=labelsNew(l);
        end
        
        % reconstruct and write out to volume
        Vi=reshape(tmp,[V.dim]);
        outName=(fullfile(studyDir{2},'encoding','glm4','groupEval_lob10')); dircheck(outName)
        V.fname=fullfile(outName,'map.nii');
        spm_write_vol(V,Vi);
        
    case 'MAP:vol2surf'
        % this function takes any labelled volume (already in SUIT space)
        % and plots to the surface
        inputMap=varargin{1}; % some options are 'Buckner_7Networks','SC1_9cluster','lob10', 'Cole_10Networks', 'SC2_90cluster' etc
        metric=varargin{2}; % 'yes' or 'no'
        
        vararginoptions({varargin{3:end}},{'border','sn'}); % option if doing individual map analysis
        
        if exist('sn'),
            inputDir=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn});
            mapName=sprintf('map_%s.nii',inputMap);
        else
            inputDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',inputMap));
            mapName='map.nii';
        end
        cd(inputDir);
        
        Vo=spm_vol(fullfile(inputDir,mapName));
        Vi=spm_read_vols(Vo);
        Vv{1}.dat=Vi;
        Vv{1}.dim=Vo.dim;
        Vv{1}.mat=Vo.mat;
        
        M=caret_suit_map2surf(Vv,'space','SUIT','stats','mode');
        M.data=round(M.data);
        % figure out colourMap
        if exist(fullfile(inputDir,'colourMap.txt'),'file'),
            cmap=load('colourMap.txt');
            cmap=cmap(:,2:4);
            cmapF=cmap/255;
        else
            cmapF=colorcube(max(M.data));
        end
        
        switch metric,
            case 'yes'
                % make paint and area files
                M=caret_struct('paint','data',M.data);
                A=caret_struct('area','data',cmapF*255,'column_name',M.paintnames,...
                    'column_color_mapping',repmat([-5 5],M.num_paintnames),'paintnames',M.paintnames);
                
                % save out paint and area files
                caret_save(fullfile(studyDir{1},caretDir,'suit_flat',sprintf('%s.paint',inputMap)),M);
                caret_save(fullfile(studyDir{1},caretDir,'suit_flat',sprintf('%s.areacolor',inputMap)),A);
            case 'no'
                if exist('border'),
                    suit_plotflatmap(M.data,'type','label','cmap',cmapF,'border',[])
                else
                    suit_plotflatmap(M.data,'type','label','cmap',cmapF) % colorcube(max(M.data))
                end
        end
    case 'MAP:make_metric' % make metric files. ONE for all SNN maps. ONE for all subjs (for final multi-task task)
        type=varargin{1}; % 'subjs' or 'SNNMaps'
        
        inputMap='SC12_10cluster';
        
        if strcmp(type,'subjs'),
            mapType=returnSubjs;
        else
            clusters=[5:24];
            for m=1:length(clusters),
                mapType{m}=sprintf('SC12_%dcluster',clusters(m));
            end
        end
        
        % looping over subjects or maps ?
        for s=1:length(mapType),
            switch type,
                case 'subjs'
                    inputDir=fullfile(studyDir{2},encodeDir,'glm4',subj_name{mapType(s)});
                    mapName=sprintf('map_%s.nii',inputMap);
                    colNames=subj_name{returnSubjs(s)};
                    outName='Multi-Task-indivSubjs';
                case 'SNNMaps'
                    inputDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType{s}));
                    mapName='map.nii';
                    colNames=sprintf('cluster%d',clusters(s));
                    outName='SNNMaps';
            end
            
            Vo=spm_vol(fullfile(inputDir,mapName));
            Vi=spm_read_vols(Vo);
            Vv{s}.dat=Vi;
            Vv{s}.dim=Vo.dim;
            Vv{s}.mat=Vo.mat;
            column_names{s}=colNames;
        end
        
        M=caret_suit_map2surf(Vv,'space','SUIT','stats','mode','column_names',column_names);
        
        % make paint and area files
        M=caret_struct('paint','data',M.data,'column_name',M.column_name);
        
        % save out paint and area files
        caret_save(fullfile(studyDir{1},caretDir,'suit_flat',sprintf('%s.paint',outName)),M);
    case 'MAP:optimal'     % figure out optimal map for multiple clusters
        % example:sc1_sc2_functionalAtlas('MAP:optimal',<subjNums>,1,6,'group')
        sn=varargin{1};     % 'group' or <subjNum>
        study=varargin{2};  % 1 or 2 or [1,2]
        K=varargin{3};      % K=numClusters (i.e. 5);
        numCount=5;         % How often the "same" solution needs to be found
        
        tol_rand = 0.90;    % Tolerance on rand coefficient to call it the same solution
        maxIter=100; % if it's not finding a similar solution - force stop at 100 iters
        
        % Set the String correctly
        studyStr = sprintf('SC%d',study);
        if length(study)>1
            studyStr='SC12'; % both studies combined
        end
        
        % Set output filename: group or indiv ?
        if strcmp(sn,'group'), % group
            sn=returnSubjs;
            outDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s_%dcluster',studyStr,K));
            dircheck(outDir);
            outName=fullfile(outDir,'SNN.mat');
        else % indiv
            outDir=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn});
            outName=fullfile(outDir,sprintf('SNN_%s_%dcluster.mat',studyStr,K));
        end
        
        % get data
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',sn,study,'build');
        
        % Intialize iterations[G
        bestErr = inf;
        bestSol = ones(size(X_C,1),1);
        iter=1; % How many iterations
        count=0;
        while iter<maxIter,
            [F,G,Info,winner]=semiNonNegMatFac(X_C,K,'threshold',0.01); % get current error
            errors(iter)=Info.error;    % record error
            randInd(iter)=RandIndex(bestSol,winner); %
            
            % Check if we have a similar solution
            if randInd(iter)>tol_rand % Similar solution
                count=count+1;       % count up one
                if (Info.error<bestErr)  % If we got slightly better - update
                    bestErr = Info.error;
                    bestSol = winner;
                    bestG   = G;
                    bestF   = F;
                    bestInfo = Info;
                end;
            else                     % Different (enough) solution
                if (Info.error<bestErr) % Is this a better solution
                    bestErr = Info.error;
                    bestSol = winner;
                    bestG   = G;
                    bestF   = F;
                    bestInfo = Info;
                    count = 0;         % first time we found this solution: reset counter
                end;
            end;
            fprintf('Error: %2.2f Rand:%2.2f, Best:%2.2f currently found %d times\n',errors(iter),randInd(iter),bestErr,count);
            if count>=numCount || iter>=maxIter,
                fprintf('Existing loop....\n');
                break;
            end;
            iter=iter+1;
        end;
        save(outName,'bestG','bestF','bestInfo','errors','randInd','iter','count','volIndx','V');
    case 'MAP:computeR2'   % temporary function to compute SSE for SNN
        study=varargin{1}; % 1 or 2 or [1,2]
        K=varargin{2};     % K=numClusters (i.e. 5);
        
        % Set the String correctly
        studyStr = sprintf('SC%d',study);
        if length(study)>1
            studyStr='SC12'; % both studies combined
        end
        
        outDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s_%dcluster',studyStr,K));
        outName=fullfile(outDir,'SNN.mat');
        load(outName);
        
        % get X
        X = sc1_sc2_functionalAtlas('EVAL:get_data','group',study,'build');
        
        % Provide fitting info
        % R2
        X_hat=bestF*bestG';
        R=X-X_hat;
        n=size(R,1);
        k2=size(bestG,2);
        SSR=nansum(R.*R);
        SST=nansum(X.*X);
        SSRadj=SSR./k2; % this is not correct
        SSTadj=SST;
        
        % R
        SXP=nansum(nansum(X.*X_hat,1));
        SPP=nansum(nansum(X_hat.*X_hat));
        
        bestInfo.R2_vox  = 1-nansum(R.*R)./nansum(X.*X);  % MK: R2=1-SSR/SST
        bestInfo.R2      = 1-nansum(SSR)./nansum(SST);
        bestInfo.R       = SXP./sqrt(nansum(SST).*SPP);   % MK: R=covar(X,X_hat)/var(X)*var(X_hat)
        bestInfo.R2adj   = 1-nansum(SSRadj)./nansum(SSTadj);
        bestInfo.error   = nansum(nansum(R.*R));
        
        fprintf('bestInfo recomputed for %d clusters \n',K)
        save(outName,'bestG','bestF','bestInfo','errors','randInd','iter','count','volIndx','V');
    case 'MAP:ICA'
        % example:sc1_sc2_functionalAtlas('MAP:optimal',<subjNums>,1,6,'group')
        sn=varargin{1};     % subject numbers
        study=varargin{2};  % 1 or 2 or [1,2]
        threshold=varargin{3};      % K=thresh (i.e. 90)
        type=varargin{4};   % 'group' or 'indiv'
        
        %numCount=5;         % How often the "same" solution needs to be found
        %tol_rand = 0.90;    % Tolerance on rand coefficient to call it the same solution
        %maxIter=100; % if it's not finding a similar solution - force stop at 100 iters
        
        % Set the String correctly
        studyStr = sprintf('SC%d',study);
        if length(study)>1
            studyStr='SC12'; % both studies combined
        end
        
        % Set output File name
        switch type,
            case 'group'
                sn=returnSubjs;
                outDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s_%dPOV',studyStr,threshold));
                dircheck(outDir);
                outName=fullfile(outDir,'ICAs.mat');
            case 'indiv'
                outDir=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn});
                outName=fullfile(outDir,sprintf('ICAs_%s_%dPOV.mat',studyStr,threshold));
        end
        
        % get data
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',sn,study,'build');
        
        threshold=threshold/100;
        [A_PW,S_PW,W_PW,winner]=pca_ica(X_C,'threshold',threshold);
        
        save(outName,'A_PW','S_PW','W_PW','volIndx','V');
    case 'MAP:visualise'
        sn=varargin{1}; % [2] or 'group'
        study=varargin{2}; % 1 or 2 or [1,2]
        anaType=varargin{3}; % 'SNN' or 'ICAs'
        K=varargin{4}; % for SNN - number of clusters (i.e. 5), for ica - thresh (i.e. 90)
        
        % figure out if ICA or SNN
        if strcmp(anaType,'SNN'),
            anaName='cluster';
        else
            anaName='POV';
        end
        
        % Set the String correctly
        studyStr = sprintf('SC%d',study);
        if length(study)>1
            studyStr='SC12'; % both studies combined
        end
        
        % figure out if individual or group
        if strcmp(sn,'group')
            outName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s_%d%s',studyStr,K,anaName),'map.nii');
            load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s_%d%s',studyStr,K,anaName),sprintf('%s.mat',anaType)));
        else
            outName=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('map_%s_%d%s.nii',studyStr,K,anaName));
            load(fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('%s_%s_%d%s.mat',anaType,studyStr,K,anaName)));
        end
        
        % transpose matrix from ICA
        if strcmp(anaType,'ICAs'),
            bestG=S_PW';
        end
        [~,groupFeat]=max(bestG,[],2);
        
        % map features on group
        %V=spm_vol(which('Cerebellum-SUIT.nii'));
        Yy=zeros(1,V.dim(1)*V.dim(2)*V.dim(3));
        Yy(1,volIndx)=groupFeat;
        Yy=reshape(Yy,[V.dim(1),V.dim(2),V.dim(3)]);
        Yy(Yy==0)=NaN;
        Vv{1}.dat=Yy;
        Vv{1}.dim=V.dim;
        Vv{1}.mat=V.mat;
        
        % save out vol of SNN feats
        exampleVol=fullfile(studyDir{2},suitDir,'glm4','s02','wdbeta_0001.nii'); % must be better way ??
        X=spm_vol(exampleVol);
        X.fname=outName;
        X.private.dat.fname=V.fname;
        spm_write_vol(X,Vv{1}.dat);
    case 'MAP:compare'  % this function is will compute randindex between maps (pairwise if more than 2 maps are provided)
        maps=varargin{1}; % give cell with any number of maps {'Cole_10Networks','Buckner_7Networks'} etc
        % example
        % sc1_sc2_functionalAtlas('MAP:compare',{'Cole_10Networks','Buckner_7Networks,'SC1_10cluster'})
        numMaps=length(maps);
        
        % get data
        [~,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',2,1,'build'); % just to get volIndx + V
        
        % make sure all parcels are sampled into the same space
        for m=1:numMaps,
            mapName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',maps{m}),'map.nii');
            [i,j,k]=ind2sub(V.dim,volIndx);
            [x,y,z]=spmj_affine_transform(i,j,k,V.mat);
            VA=spm_vol(mapName);
            [i1,j1,k1]=spmj_affine_transform(x,y,z,inv(VA.mat));
            L(m,:)=spm_sample_vol(VA,i1,j1,k1,0);
        end
        
        for i=1:numMaps,
            for j=1:numMaps,
                RI(i,j)=RandIndex(L(i,:),L(j,:));
            end
        end
        
        % print out values for all pairwise RI
        idx=combnk(1:length(maps),2);
        for i=1:length(idx),
            fprintf('RI:%s with %s is %2.2f \n',char(maps(idx(i,1))),char(maps(idx(i,2))),RI(idx(i,1),idx(i,2)))
        end
        varargout={RI};
    case 'MAP:PLOT:ICA' % evaluate POV (ICA) as a function of clusters
        mapDir=varargin{1}; % {'SC1','SC2','SC12'}
        mapType=varargin{2};
        POV=varargin{3}; % [50:5:95,96:100]
        
        S=[];
        for s=1:length(mapDir),
            for m=1:length(mapType),
                load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s_%s',mapDir{s},mapType{m}),'ICAs.mat'));
                T.K(m,1)=size(S_PW,1);
                T.mapName{m,1}=mapType{m};
                T.m(m,1)=m;
                T.type{m,1}=mapDir{s};
            end
            T.POV=POV';
            S=addstruct(S,T);
        end
        figure()
        lineplot(S.POV,S.K,'split',S.type,'leg','auto','style_thickline2x3','linewidth',4)
        ylabel('clusters')
        xlabel('variance')
    case 'MAP:PLOT:SNN'
        mapDir=varargin{1}; % {'SC1','SC2','SC12'}
        mapType=varargin{2};
        K=varargin{3}; % [5:24]
        
        vararginoptions({varargin{4:end}},{'CAT'}); % option if doing individual map analysis
        
        S=[];
        for s=1:length(mapDir),
            for m=1:length(mapType),
                load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s_%s',mapDir{s},mapType{m}),'SNN.mat'));
                T.R2(m,1)=bestInfo.R2;
                T.R2adj(m,1)=bestInfo.R2adj;
                T.mapName{m,1}=mapType{m};
                T.m(m,1)=m;
                T.type{m,1}=mapDir{s};
            end
            T.K=K';
            S=addstruct(S,T);
        end
        figure()
        lineplot(S.K,S.R2,'split',S.type,'CAT',CAT)
        hold on
        CAT.linecolor={'k'};
        CAT.markercolor={'k'};
        CAT.markerfill={'k'};
        lineplot(S.K,S.R2adj,'split',S.type,'CAT',CAT)
        
    case 'EVAL:get_data'
        sn=varargin{1}; % Subj numbers to include
        study=varargin{2}; % 1 or 2 or [1,2]
        type=varargin{3}; % 'build' or 'eval'. For build - we get group data. For eval - we get indiv data
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % input could be 'group' or <subjNums>
        if strcmp(sn,'group'),
            sn=returnSubjs;
        end
        
        % load data
        UFullAvrgAll=[];
        for f=1:length(study),
            load(fullfile(studyDir{study(f)},encodeDir,'glm4','cereb_avrgDataStruct.mat'));
            
            D1=getrow(D,D.StudyNum==f);
            idx=D1.condNum(D1.overlap==1); % get index for unique tasks
            % format data
            for s=1:length(sn),
                indx=T.SN==sn(s);
                UFull=T.data(indx,:);
                [numConds,numVox]=size(UFull);
                numConds=numConds/2;
                % get average across sessions
                for c=1:numConds,
                    UFullAvrg(c,:,s)=nanmean(UFull([c,c+numConds],:),1);
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
        varargout={X_C,volIndx,V,sn};
    case 'EVAL:visualise_data' % this just visualises the volume (averaged betas for unique taskConds) for each dataset
        sn=varargin{1}; % take example subject
        study=varargin{2}; % take example study
        mapType=varargin{3}; % take example mapType
        
        % load map
        inDir=fullfile(studyDir{2},encodeDir,'glm4');
        mapName=fullfile(inDir,subj_name{sn},sprintf('map_%s.nii',mapType));
        
        % load in func data to evaluate
        load(fullfile(studyDir{study},'encoding','glm4','cereb_avrgDataStruct.mat'));
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        D1=getrow(D,D.StudyNum==study);
        
        idx=D1.condNum(D1.overlap==0); % get index for unique tasks
        
        % Now get the parcellation sampled into the same space
        [i,j,k]=ind2sub(V.dim,volIndx);
        [x,y,z]=spmj_affine_transform(i,j,k,V.mat);
        VA= spm_vol(mapName);
        [i1,j1,k1]=spmj_affine_transform(x,y,z,inv(VA.mat));
        Parcel = spm_sample_vol(VA,i1,j1,k1,0);
        voxIn = Parcel>0;
        X=spm_read_vols(VA);
        voxNum = find(X>0);
        clear i j k x y z i1 j1 k1; % Free memory
        
        for s=sn,
            for c=1:length(idx),
                i1(c) = find(T.SN==s & T.sess==1 & T.cond==idx(c));
                i2(c) = find(T.SN==s & T.sess==2 & T.cond==idx(c));
            end
            D=(T.data(i1,voxIn)+T.data(i2,voxIn))/2;
        end
        % average across all unique taskConds
        Vavrg=nanmean(D,1);
        
        % make volume
        Yy=zeros(1,VA.dim(1)*VA.dim(2)*VA.dim(3));
        Yy(1,voxNum)=Vavrg;
        Yy=reshape(Yy,[VA.dim(1),VA.dim(2),VA.dim(3)]);
        Yy(Yy==0)=NaN;
        
        % write out volume
        VA.fname=fullfile(inDir,subj_name{sn},sprintf('uniqueTasks_SC%d.nii',study));
        VA.private.dat.fname=VA.fname;
        spm_write_vol(VA,Yy);
    case 'EVAL:crossval'
        sn=varargin{1}; % 'group' or <subjNum>
        mapType=varargin{2}; % options are 'lob10','lob26','Buckner_17Networks','Buckner_7Networks', 'Cole_10Networks','SC<studyNum>_<num>cluster'
        data=varargin{3}; % evaluating data from study [1] or [2] ?
        condType=varargin{4}; % 'unique' or 'all'. Are we evaluating on all taskConds or just those unique to either sc1 or sc2 ?
        
        % set outName - group or indiv ?
        if strcmp(sn,'group'),
            % load in map
            mapName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'map.nii');
            outName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_%s.mat',data,condType));
        else
            mapName=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('map_%s.nii',mapType));
            outName=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('%s_spatialBoundfunc%d_%s.mat',mapType,data,condType));
        end
        
        % load in func data to test (e.g. if map is sc1; func data should
        % be sc2)
        load(fullfile(studyDir{data},'encoding','glm4','cereb_avrgDataStruct.mat'));
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        D1=getrow(D,D.StudyNum==data);
        switch condType,
            case 'unique'
                % if funcMap - only evaluate unique tasks in sc1 or sc2
                idx=D1.condNum(D1.overlap==0); % get index for unique tasks
            case 'all'
                idx=D1.condNum;
        end
        
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
        % Now calculate the uncrossvalidated and crossvalidated
        % Estimation of the correlation for each subject
        for s=unique(T.SN)';
            for c=1:length(idx),
                i1(c) = find(T.SN==s & T.sess==1 & T.cond==idx(c));
                i2(c) = find(T.SN==s & T.sess==2 & T.cond==idx(c));
            end
            D=(T.data(i1,voxIn)+T.data(i2,voxIn))/2;
            
            fprintf('%d cross\n',s);
            R.SN = ones(length(R.N),1)*s;
            R.corr = mva_spatialCorr(T.data([i1;i2],voxIn),BIN,...
                'CrossvalPart',T.sess([i1;i2],1),'excludeNegVoxels',1);
            R.crossval = ones(length(R.corr),1);
            RR = addstruct(RR,R);
            fprintf('%d correl\n',s);
            R.corr=mva_spatialCorr(D,BIN);
            R.crossval = zeros(length(R.corr),1);
            RR = addstruct(RR,R);
        end;
        save(outName,'-struct','RR');
    case 'EVAL:crossval_parcels'
        sn=varargin{1}; % 'group' or <subjNum>
        mapType=varargin{2}; % options are 'lob10','lob26','Buckner_17Networks','Buckner_7Networks', 'Cole_10Networks','SC<studyNum>_<num>cluster'
        data=varargin{3}; % evaluating data from study [1] or [2] ?
        condType=varargin{4}; % 'unique' or 'all'. Are we evaluating on all taskConds or just those unique to either sc1 or sc2 ?
        type=varargin{5}; % 'group' or 'indiv'
        
        switch type,
            case 'group'
                % load in map
                mapName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'map.nii');
                outName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_%s_parcels.mat',data,condType));
            case 'indiv'
                mapName=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('map_%s.nii',mapType));
                outName=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('%s_spatialBoundfunc%d_%s.mat',mapType,data,condType));
        end
        
        % load in func data to test (e.g. if map is sc1; func data should
        % be sc2)
        load(fullfile(studyDir{data},'encoding','glm4','cereb_avrgDataStruct.mat'));
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        D1=getrow(D,D.StudyNum==data);
        switch condType,
            case 'unique'
                % if funcMap - only evaluate unique tasks in sc1 or sc2
                idx=D1.condNum(D1.overlap==0); % get index for unique tasks
            case 'all'
                idx=D1.condNum;
        end
        
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
        R=[];
        [BIN,R]=mva_spatialCorrBin(XYZ(:,voxIn),'Parcel',Parcel(1,voxIn));
        
        clear XYZ i k l x y z i1 j1 k1 VA; % Free memory
        % Now calculate the uncrossvalidated
        % Estimation of the correlation for each subject
        for s=unique(T.SN)';
            for c=1:length(idx),
                i1(c) = find(T.SN==s & T.sess==1 & T.cond==idx(c));
                i2(c) = find(T.SN==s & T.sess==2 & T.cond==idx(c));
            end
            D=(T.data(i1,voxIn)+T.data(i2,voxIn))/2;
            
            fprintf('%d cross\n',s);
            [~,~,binCorrVox]=mva_spatialCorr(D,BIN,'saveVox','yes'); % this function was changed to add binCorrVox as output
            
            R.corr=binCorrVox';
            R.SN=repmat(s,size(R.corr,1),1);
            R.crossval=repmat(0,size(R.corr,1),1); % no across-sess crossval (was taking too long !)
            RR=addstruct(RR,R);
        end;
        voxIn=volIndx(voxIn);
        save(outName,'RR','voxIn');
    case 'EVAL:averageK' % make new 'spatialBoundfunc4.mat' struct. [4] - average eval corr across studies
        mapType=varargin{1}; % ex. '6cluster' (for snn) or '90POV' (for ica)
        condType=varargin{2}; % evaluating on 'unique' or 'all' taskConds ?
        type=varargin{3}; % 'group' or 'indiv' ?
        
        studyType=[1,2]; % test on one dataset
        evalType=[2,1]; % evaluate on the other
        R=[];
        
        vararginoptions({varargin{4:end}},{'sn'}); % option if doing individual map analysis
        
        switch type,
            case 'group'
                for i=1:2,
                    T=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_SC%d_%s',studyType(i),mapType),sprintf('spatialBoundfunc%d_%s.mat',evalType(i),condType)));
                    T.studyNum=repmat([i],length(T.SN),1);
                    R=addstruct(R,T);
                end
                outDir=fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_SC2_%s',mapType),sprintf('spatialBoundfunc4_%s.mat',condType));
            case 'indiv'
                for i=1:2,
                    T=load(fullfile(studyDir{2},'encoding','glm4',subj_name{sn},sprintf('SC%d_%s_spatialBoundfunc%d_%s.mat',studyType(i),mapType,evalType(i),condType)));
                    T.studyNum=repmat([i],length(T.SN),1);
                    R=addstruct(R,T);
                end
                outDir=fullfile(studyDir{2},'encoding','glm4',subj_name{sn},sprintf('SC2_%s_spatialBoundfunc4_%s.mat',mapType,condType));
        end
        R=rmfield(R,{'distmin','distmax','N'});
        % get average of both structures here
        A=tapply(R,{'bin','SN','bwParcel','crossval'},{'corr'});
        
        % distances are diff across evals so need to get dist per bin:
        for b=1:length(unique(R.bin)),
            dist=mode(round(R.dist(R.bin==b)));
            idx=find(A.bin==b);
            A.dist(idx,1)=dist;
        end
        
        save(outDir,'-struct','A');
    case 'EVAL:averageOther'
        mapType=varargin{1}; % options are 'lob10','Buckner_7Networks','Buckner_17Networks','Cole_10Networks','SC12_<K>cluster' or 'SC12_<thresh>POV' etc'
        condType=varargin{2}; % evaluating on 'unique' or 'all' taskConds ?
        
        R=[];
        for i=1:2,
            T=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_%s.mat',i,condType)));
            T.studyNum=repmat([i],length(T.SN),1);
            R=addstruct(R,T);
        end
        R=rmfield(R,{'distmin','distmax','N'});
        % get average of both structures here
        A=tapply(R,{'bin','SN','bwParcel','crossval','dist'},{'corr'});
        
        % distances are diff across evals so need to get dist per bin:
        for b=1:length(unique(R.bin)),
            dist=mode(round(R.dist(R.bin==b)));
            idx=find(A.bin==b);
            A.dist(idx,1)=dist;
        end
        
        outDir=fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc4_%s.mat',condType));
        save(outDir,'-struct','A');
    case 'EVAL:PLOT:CURVES'
        mapType=varargin{1}; % options are 'lob10','lob26','bucknerRest','SC<studyNum>_<num>cluster', or 'SC<studyNum>_POV<num>'
        data=varargin{2}; % evaluating data from study [1] or [2], both [3] or average of [1] and [2] after eval [4]
        type=varargin{3}; % 'group' or 'leaveOneOut' or 'indiv'
        crossval=varargin{4}; % [0] - no crossval; [1] - crossval
        condType=varargin{5}; % evaluating on 'all' or 'unique' taskConds ??
        
        vararginoptions({varargin{6:end}},{'CAT','sn'}); % option if doing individual map analysis
        
        switch type,
            case 'group'
                T=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_%s.mat',data,condType)));
            case 'indiv'
                T=[];
                for s=1:length(sn)
                    tmp=load(fullfile(studyDir{2},'encoding','glm4',subj_name{sn(s)},sprintf('%s_spatialBoundfunc%d_%s.mat',mapType,data,condType)));
                    S=getrow(tmp,tmp.SN~=sn(s)); % only predicting other subjs
                    T=addstruct(T,S);
                end
            otherwise
                fprintf('no such case')
        end
        
        xyplot(T.bin,T.corr,T.bin,'split',T.bwParcel,'subset',T.crossval==crossval & T.dist<=35,'CAT',CAT,'leg',{'within','between'});
        
        % do stats (over all bins)
        C=getrow(T,T.crossval==1 & T.dist<=35); % only crossval and dist<35
        S=tapply(C,{'bwParcel','SN'},{'corr'});
        fprintf('overall \n')
        ttest(S.corr(S.bwParcel==1), S.corr(S.bwParcel==0),2,'independent');
        
        % do stats (per bin)
        bins=unique(C.bin);
        for b=1:length(bins),
            fprintf('bin %d \n',b)
            ttest(C.corr(C.bwParcel==1 & C.bin==b), C.corr(C.bwParcel==0 & C.bin==b),2,'independent');
        end
        
        % summary stats
        x1=nanmean(S.corr(S.bwParcel==0));x2=nanmean(S.corr(S.bwParcel==1));
        SEM1=std(S.corr(S.bwParcel==0))/sqrt(length(returnSubjs));SEM2=std(S.bwParcel==1)/sqrt(length(returnSubjs));
        fprintf('average within corr is %2.2f; CI:%2.2f-%2.2f \n average between corr is %2.2f; CI:%2.2f-%2.2f \n',...
            nanmean(S.corr(S.bwParcel==0)),x1-(1.96*SEM1),x1+(1.96*SEM1),nanmean(S.corr(S.bwParcel==1)),...
            x2-(1.96*SEM2),x2+(1.96*SEM2));
    case 'EVAL:PLOT:DIFF'
        mapType=varargin{1}; % {'lob10','bucknerRest','atlasFinal9'}
        data=varargin{2}; % evaluating data from study [1] or [2] or [3]?
        type=varargin{3}; % 'group' or 'indiv'
        crossval=varargin{4}; % [0]-no crossval; [1]-crossval
        condType=varargin{5}; % evaluating on 'unique' or 'all' taskConds ??
        
        vararginoptions({varargin{6:end}},{'CAT','sn'}); % option if plotting individual map analysis
        
        P=[];
        for m=1:length(mapType),
            switch type,
                case 'group'
                    T=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType{m}),sprintf('spatialBoundfunc%d_%s.mat',data(m),condType)));
                case 'indiv'
                    T=[];
                    for s=1:length(sn)
                        tmp=load(fullfile(studyDir{2},'encoding','glm4',subj_name{sn(s)},sprintf('%s_spatialBoundfunc%d_%s.mat',mapType{m},data(m),condType)));
                        S=getrow(tmp,tmp.SN==sn(s)); % only predicting the same subject !
                        T=addstruct(T,S);
                    end
            end
            A=getrow(T,T.crossval==crossval);
            A.type=repmat({sprintf('%d.%s',m,mapType{m})},length(A.bin),1);
            A.m=repmat(m,length(A.bin),1);
            P=addstruct(P,A);
            clear A
        end
        
        % plot boxplot of different clusters
        W=getrow(P,P.bwParcel==0); % within
        B=getrow(P,P.bwParcel==1); % between
        W.diff=W.corr-B.corr;
        
        lineplot(W.bin,W.diff,'split',W.type,'subset',W.dist<=35,'CAT',CAT,'leg','auto');
        
        % do stats (integrate over spatial bins)
        C=getrow(W,W.dist<=35);
        S=tapply(C,{'m','SN','type'},{'corr'});
        
        for m=1:length(mapType),
            y(m,:)=S.corr(S.m==m);
        end
        
        % do F test first - are there any differences across maps ?
        
        % do individual t-tests
        %         idx=combnk([1:length(mapType)],2)';
        %         for t=1:length(idx),
        %             fprintf('%s:%s \n',mapType{idx(1,t)},mapType{idx(2,t)});
        %             ttest(y(idx(1,t),:),y(idx(2,t),:),2,'independent')
        %         end
    case 'EVAL:PLOT:Parcels'
        mapType=varargin{1}; % options are 'lob10','lob26','bucknerRest','SC<studyNum>_<num>cluster', or 'SC<studyNum>_POV<num>'
        data=varargin{2}; % evaluating data from study [1] or [2]
        type=varargin{3}; % 'group' or 'leaveOneOut' or 'indiv'
        condType=varargin{4}; % evaluating on 'all' or 'unique' taskConds ??
        
        bins=[1:7]; % we only want bins 1-7 (<35 mm)
        vararginoptions({varargin{6:end}},{'sn'}); % option if doing individual map analysis
        
        [~,~,V]=sc1_sc2_functionalAtlas('EVAL:get_data',2,data,'build');  % put any subjNum (doesn't matter here)
        
        switch type,
            case 'group'
                load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_%s_parcels.mat',data,condType)));
            case 'indiv'
                load(fullfile(studyDir{2},'encoding','glm4',subj_name{sn},sprintf('%s_spatialBoundfunc%d_%s.mat',mapType,data,condType)));
            otherwise
                fprintf('no such case')
        end
        
        for b=1:length(bins),
            % compute diff for diff bins
            within=getrow(RR,RR.bwParcel==0 & RR.bin==bins(b));
            between=getrow(RR,RR.bwParcel==1 & RR.bin==bins(b));
            diff=within.corr-between.corr;
            
            columnName{b}=sprintf('bin%d-%dmm',bins(b),unique(round(within.dist)));
            
            % map features on group
            Yy=zeros(1,V.dim(1)*V.dim(2)*V.dim(3));
            Yy(1,voxIn)=nanmean(diff,1);
            Yy=reshape(Yy,[V.dim(1),V.dim(2),V.dim(3)]);
            Yy(Yy==0)=NaN;
            Vv{b}.dat=Yy;
            Vv{b}.dim=V.dim;
            Vv{b}.mat=V.mat;
            clear Yy diff
        end
        
        % map vol2surf
        M=caret_suit_map2surf(Vv,'space','SUIT','stats','nanmean','column_names',columnName);
        
        caret_save(fullfile(studyDir{2},caretDir,'suit_flat','glm4',sprintf('diff_%s.metric',mapType)),M);
    case 'EVAL:visualiseBounds' % Generates metric file and border file for the pacellation
        mapType = varargin{1}; % 'lob10','Buckner_7Networks','SC12_10cluster' etc
        
        EvalDir = fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType));
        SurfDir = fullfile(studyDir{1},'surfaceCaret','suit_flat');
        load(fullfile(EvalDir,'boundaries.mat'));
        
        
        % Map the clusters
        V.dat=zeros([V.dim(1) V.dim(2) V.dim(3)]);
        V.dat(volIndx)=Cluster;
        Mcl=suit_map2surf(V,'space','SUIT','stats','mode','stats',@mode);
        
        % Map the parcel
        V.dat(volIndx)=Parcel;
        Mpa=suit_map2surf(V,'space','SUIT','stats','mode','stats',@mode);
        
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
        
        % Evaluate the strength of each border
        T=load(fullfile(EvalDir,'BoundariesFunc3_all.mat'));
        for  i=1:size(Edge,1)
            % Make sure that the bin is calcualted both for within and
            % between
            A=pivottable(T.bin,T.bwParcel,T.corr,'nanmean','subset',T.crossval==1 & ...
                T.Edge(:,1)==Edge(i,1) & T.Edge(:,2)==Edge(i,2));
            EdgeWeight(i,1)=nanmean(A(:,1)-A(:,2));  % Difference within - between
        end;
        
        % check if a color map is provided
        if exist(fullfile(EvalDir,'colourMap.txt'),'file'),
            cmap=load('colourMap.txt');
            cmap=cmap(:,2:4);
            cmapF=cmap/255;
        else
            cmapF=colorcube(max(Parcel));
        end
        
        % Make the plot
        suit_plotflatmap(Mpa,'type','label','border',[],'cmap',cmapF);
        hold on;
        LineWeight=EdgeWeight*200;
        % LineWeight(LineWeight>20)=20;
        for  b=1:length(Border)
            if (Border(b).numpoints>0 & EdgeWeight(b)>0)
                p=plot(Border(b).data(:,1),Border(b).data(:,2),'k.');
                set(p,'MarkerSize',LineWeight(b));
            end;
        end;
        hold off;
        
        keyboard;
        
    case 'STRENGTH:get_bound'
        % This goes from a group parcellation map and generates a
        % structure of clusters and boundaries from the volume
        mapType = varargin{1};
        
        bDir    = fullfile(studyDir{2},'encoding','glm4');
        EvalDir = fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType));
        load(fullfile(bDir,'cereb_avrgDataStruct.mat'));
        mapName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'map.nii');
        
        % Get the parcellation data
        [i,j,k]=ind2sub(V.dim,volIndx);
        [x,y,z]=spmj_affine_transform(i,j,k,V.mat);
        VA= spm_vol(mapName);
        [i1,j1,k1]=spmj_affine_transform(x,y,z,inv(VA.mat));
        Parcel = spm_sample_vol(VA,i1,j1,k1,0);
        
        % Find unique clusters for each parcel
        numParcel = max(Parcel); % was 28 (from Joern's code)
        Cluster = nan(size(Parcel));
        n=1;
        coords=[i;j;k];     % Voxel coordinates in original image
        for p=1:numParcel
            indx = find(Parcel==p);
            A= spm_clusters(coords(:,indx));
            numCluster=max(A);
            for c=1:numCluster
                clInd = (A==c);
                N=sum(clInd);
                if N>=5  % ignore clusters of smaller than 5
                    Cluster(indx(clInd))=n;
                    n=n+1;
                end;
            end;
        end;
        
        % Check how assignment went
        pivottable(Cluster',Parcel',Cluster','length','subset',~isnan(Cluster'));
        
        % Now detect boundaries between adjacent parcels
        numCluster = max(Cluster);
        n=1;
        for i=1:numCluster
            for j=i+1:numCluster
                D=surfing_eucldist(coords(:,Cluster==i),coords(:,Cluster==j));
                if min(D(:))< 1.4 % direct connectivity scheme
                    Edge(n,:) = [i j];
                    n=n+1;
                end;
            end;
        end;
        
        % Visualise using graph toolbox
        G = graph(Edge(:,1),Edge(:,2));
        plot(G)
        save(fullfile(EvalDir,'boundaries.mat'),'V','volIndx','Parcel','Cluster','coords','Edge');
    case 'STRENGTH:fullEval'
        mapType = varargin{1};
        EvalDir = fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType));
        
        T1=sc1_sc2_functionalAtlas('STRENGTH:eval_bound',mapType,1,'unique');
        T2=sc1_sc2_functionalAtlas('STRENGTH:eval_bound',mapType,2,'unique');
        T1.study = ones(length(T1.SN),1)*1;
        T2.study = ones(length(T1.SN),1)*2;
        T=addstruct(T1,T2);
        save(fullfile(EvalDir,'BoundariesFunc3_all.mat'),'-struct','T');
        varargout={T};
    case 'STRENGTH:eval_bound'
        mapType = varargin{1}; % 'SC12_10cluster','Buckner_7Networks'
        study   = varargin{2}; % evaluating data from study [1] or [2] ?
        condType = varargin{3}; % 'unique' or 'all'
        
        spatialBins = [0:3:35];
        
        bDir    = fullfile(studyDir{2},'encoding','glm4');
        EvalDir = fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType));
        load(fullfile(EvalDir,'boundaries.mat'));
        numBins = length(spatialBins)-1;
        
        % Get the condition numbers
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        D1=getrow(D,D.StudyNum==study);
        switch condType,
            case 'unique'
                % if funcMap - only evaluate unique tasks in sc1 or sc2
                idx=D1.condNum(D1.overlap==0); % get index for unique tasks
            case 'all'
                idx=D1.condNum;
        end;
        
        % Load activity data
        load(fullfile(studyDir{study},encodeDir,'glm4','cereb_avrgDataStruct.mat'));
        RR=[];    % Output structure.
        
        % Now build a structure all boundaries
        for i=1:size(Edge,1)
            indx = Cluster==Edge(i,1) | Cluster==Edge(i,2);  % Get all the involved voxels
            fprintf('Edge %d %d\n',Edge(i,1),Edge(i,2));
            
            % Determine the spatial bins for this pair of regions
            [BIN,R]=mva_spatialCorrBin(coords(:,indx),'Parcel',Cluster(indx)','spatialBins',spatialBins);
            N=length(R.N);
            
            % Now determine the correlation for each subject
            for s=unique(T.SN)';
                fprintf('Subj:%d\n',s);
                for c=1:length(idx),
                    i1(c) = find(T.SN==s & T.sess==1 & T.cond==idx(c));
                    i2(c) = find(T.SN==s & T.sess==2 & T.cond==idx(c));
                end
                D=(T.data(i1,indx)+T.data(i2,indx))/2; % average data
                fprintf('%d cross\n',s);
                R.Edge = repmat(Edge(i,:),N,1);
                R.SN = ones(N,1)*s;
                R.corr = mva_spatialCorr(T.data([i1;i2],indx),BIN,...
                    'CrossvalPart',T.sess([i1;i2],1),'excludeNegVoxels',1,'numBins',N);
                R.crossval = ones(N,1);
                if (length(R.corr)~=N)
                    keyboard;
                end;
                RR = addstruct(RR,R);
                fprintf('%d correl\n',s);
                R.corr=mva_spatialCorr(D,BIN,'numBins',N);
                R.crossval = zeros(N,1);
                if (length(R.corr)~=N)
                    keyboard;
                end;
                RR = addstruct(RR,R);
            end;
        end;
        varargout={RR};
    case 'STRENGTH:visualise_bound'
        mapType = varargin{1};
        
        EvalDir = fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType));
        SurfDir = fullfile(studyDir{1},'surfaceCaret','suit_flat');
        load(fullfile(EvalDir,'boundaries.mat'));
        
        
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
        
        % Evaluate the strength of each border
        T=load(fullfile(EvalDir,'BoundariesFunc3_all.mat'));
        for  i=1:size(Edge,1)
            % Make sure that the bin is calcualted both for within and
            % between
            A=pivottable(T.bin,T.bwParcel,T.corr,'nanmean','subset',T.crossval==1 & ...
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
        
        % Make the plot
        suit_plotflatmap(Mpa,'type','label','border',[],'cmap',cmapF);
        hold on;
        LineWeight=EdgeWeight*200;
        LineWeight(LineWeight>20)=20;
        for  b=1:length(Border)
            if (Border(b).numpoints>0 & EdgeWeight(b)>0)
                p=plot(Border(b).data(:,1),Border(b).data(:,2),'k.');
                set(p,'MarkerSize',LineWeight(b));
                weights(b)=EdgeWeight(b);
            end;
        end;
        
        %         % plot percent (%) of variance explained by each functional boundary ?
        %         for b=1:length(weights),
        %             weightPer=(100/sum(weights))*weights(b);
        %             if (Border(b).numpoints>0 & EdgeWeight(b)>0)
        %                 p=text(double(Border(b).data(1,1)),double(Border(b).data(1,2)),sprintf('%s%%',num2str(round(weightPer))));
        %                 set(p,'FontSize',20);
        %             end
        %         end
        hold off;
        
    case 'ENCODE:get_features'
        mapType=varargin{1};
        
        D=dload(fullfile(baseDir,'featureTable_jd.txt')); % Read feature table
        S=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt')); % List of task conditions
        
        load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'SNN.mat'));
        W=pivottablerow(S.condNumUni,bestF,'mean(x,1)'); % we only want unique taskConds
        
        % get new condNames (unique only)
        condNames=[S.condNames(S.StudyNum==1);S.condNames(S.StudyNum==2 & S.overlap==0)];
        
        % Make the feature matrix
        D.LeftHand    = D.leftHandPresses ./D.duration;
        D.RightHand   = D.rightHandPresses ./D.duration;
        D.Saccade    = D.saccades./D.duration;
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
        %         A = -eye(numFeat);
        %         b = zeros(numFeat,1);
        for p=1:numClusters,
            % U(:,p) =
            % cplexqp(XX+lambda(2)*eye(numFeat),ones(numFeat,1)*lambda(1)-XY(:,p),A,b);
            % % I don't have the IBM cplexqp routine
            U(:,p) = (XX + eye(numFeat)*lambda(1))\(XY(:,p)); % just normal ridge regression
        end;
        
        % Get corr between feature weights
        C=corr(F,W);
        
        % Present the list of the highest three correlation for each
        % cluster
        for i=1:numClusters,
            [a,b]=sort(C(:,i),'descend');
            B.clusters(i,1)=i;
            % get 3 highest corrs
            for f=1:3,
                B.featNames{i,f}=FeatureNames{b(f)};
                B.featIdx(i,f)=b(f);
                B.featCorrs(i,f)=a(f);
            end
        end;
        
        varargout={B,F,W,C,condNames,FeatureNames};
        
        %                 % Make word lists for word map
        %                 figure(2);
        %                 DD = U;WORD=FeatureNames;
        %                 DD(DD<0)=0;
        %                 DD=DD./max(DD(:));
        %                 for j=1:numClusters,
        %                     subplot(3,4,j);
        %                     set(gca,'XLim',[0 2.5],'YLim',[-0.2 1.2]);
        %                     title(sprintf('Cluster%d',j))
        %                     for i=1:numFeat
        %                         if (DD(i,j)>0)
        %                             siz=ceil(DD(i,j)*20);
        %                             text(unifrnd(0,1,1),unifrnd(0,1,1),WORD{i},'FontSize',siz);
        %                         end;
        %                     end;
        %                 end;
        %                 set(gcf,'PaperPosition',[1 1 60 30]);
        %                 wysiwyg;
    case 'ENCODE:project_featSpace'
        mapType=varargin{1};
        toPlot=varargin{2}; % 'winner' or 'all' or 'featMatrix'
        
        % get features
        [B,F,~,C,condNames,FeatureNames]=sc1_sc2_functionalAtlas('ENCODE:get_features',mapType);
        
        switch toPlot,
            case 'all'
                imagesc_rectangle(C);
                caxis([0 1]);
                t=set(gca,'Ytick',[1:length(FeatureNames)]','YTickLabel',FeatureNames');
                t.Color='white';
                colorbar
            case 'winner'
                % get cluster colours
                cmap=load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'colourMap.txt'));
                cmap=cmap/255;
                
                % figure out where to position features
                for i=1:size(B.featNames,1),
                    subplot(3,4,i);
                    set(gca,'XLim',[0 2.5],'YLim',[-0.2 1.2]);
                    text(0,1.2,sprintf('Network %d',i),'FontSize',20,'Color',cmap(i,2:4));
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
        mapType=varargin{1};
        toPlot=varargin{2}; % 1 is [4,10] - left & right hand tasks etc
        
        vararginoptions({varargin{3:end}},{'CAT'});
        
        % get features
        [B,F,W,~,condNames]=sc1_sc2_functionalAtlas('ENCODE:get_features',mapType);
        
        D=dload(fullfile(baseDir,'motorFeats.txt')); % Read feature table
        S=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt')); % List of task conditions
        
        % project back to task-space
        load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'SNN.mat'));
        %         W=bestF;
        L=W'*W;
        I=diag(diag(sqrt(L))); % diag not sum
        X=W/I;
        
        if strcmp(toPlot,'all'),
            imagesc_rectangle(X);
            caxis([0 1]);
            t=set(gca,'Ytick',[1:length(condNames)]','YTickLabel',condNames');
            t.Color='white';
            colorbar
        else
            % load in colourMap
            cmap=load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'colourMap.txt'));
            cmap=cmap(:,2:4)/255;
            
            % feature loadings on network
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
            CAT.markercolor=colourIdx;
            CAT.markerfill=colourIdx;
            CAT.labelcolor=colourIdx;
            
            scatterplot(X(:,toPlot(1)),X(:,toPlot(2)),'label',condNames,'regression','linear','intercept',0,'draworig','CAT',CAT);
            xlabel(sprintf('Network%d',toPlot(1)));ylabel(sprintf('Network%d',toPlot(2)));
        end
        
    case 'AXES:eigenValues'
        % Aesthetics
        CAT.markersize=8;
        CAT.markertype='^';
        CAT.linewidth=4;
        CAT.linecolor={'r','k'};
        CAT.markercolor={'r','k'};
        CAT.markerfill={'r','k'};
        
        sc1_sc2_functionalAtlas('REPRESENTATION:eigDecomp',CAT)
        set(gcf,'units','points','position',[5,5,1000,1000])
        set(gca,'FontSize',14);
        ylabel('Eigenvalues')
        xlabel('# of dimensions')
    case 'AXES:reliability'
        toPlot=varargin{1}; % 'D' (for distance) or 'A' (for activity)
        toPlotName=varargin{2}; % 'distance reliability' or 'activity reliability'
        
        sc1_sc2_functionalAtlas(sprintf('PLOT:reliability%s',toPlot))
        
        % Labelling
        set(gca,'FontSize',18);
        set(gca,'XTickLabel',{'SC1','SC2','across'})
        ylabel('Activity Correlation (R)');
        title(toPlotName)
        set(gcf,'units','points','position',[5,5,1000,1000])
    case 'AXES:RDM'
        type=varargin{1}; % 'all' or 'none' - removeMotor
        
        % aesthetics
        sc1_sc2_functionalAtlas('REPRESENTATION:RDM',type)
        set(gcf,'units','points','position',[5,5,1000,1000])
    case 'AXES:MDS' % plots MDS plot
        type=varargin{1}; % 'all' or 'none' - removeMotor
        
        % aesthetics
        colour={[1 0 0],[0 1 0],[0 0 1],[0.3 0.3 0.3],[1 0 1],[1 1 0],[0 1 1],...
            [0.5 0 0.5],[0.8 0.8 0.8],[.07 .48 .84],[.99 .76 .21],[.11 .7 .68],...
            [.39 .74 .52],[.21 .21 .62],[0.2 0.2 0.2],[.6 .6 .6],[.3 0 .8],[.8 0 .4],...
            [0 .9 .2],[.1 .3 0],[.2 .4 0],[.63 0 .25],[0 .43 .21],[.4 0 .8]};
        CAT.markersize=10;
        CAT.markertype='o';
        CAT.labelsize=14;
        
        sc1_sc2_functionalAtlas('REPRESENTATION:MDS',type,'CAT',CAT,'colour',colour)
        set(gcf,'units','points','position',[5,5,1000,1000])
        axis equal;
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'Box','on');
        xlabel('PC1');ylabel('PC2');zlabel('PC3')
        view([81 9]);
    case 'AXES:map' % takes any volume and plots on surface
        toPlot=varargin{1}; % 'Cole_10Networks', 'Buckner_7Networks','SC12_10cluster'
        toPlotName=varargin{2}; % 'Cole10', 'Buckner7', ...
        
        vararginoptions({varargin{3:end}},{'border','sn'}); % option if doing individual map analysis
        
        % plot borders ?
        if exist('border'),
            sc1_sc2_functionalAtlas('MAP:vol2surf',toPlot,'no','border',[])
        else
            sc1_sc2_functionalAtlas('MAP:vol2surf',toPlot,'no')
        end
        
        % one subject or group ?
        if exist('sn'),
            sc1_sc2_functionalAtlas('MAP:vol2surf',toPlot,'no','sn',sn)
        end
        
        title(toPlotName);
        set(gcf,'units','points','position',[5,5,1000,1000])
        set(gca,'Xticklabel',[],'Yticklabel',[])
        axis off
    case 'AXES:allSNN_similarity' % matrix of rand indices between all SNN maps (SC12)
        K=[5:25];
        for k=1:length(K),
            toPlot{k}=sprintf('SC12_%dcluster',K(k));
            toPlotNames{k}=sprintf('%d',K(k));
        end
        
        [RI]=sc1_sc2_functionalAtlas('MAP:compare',toPlot);
        
        figure()
        imagesc_rectangle(RI,'YDir','reverse');
        caxis([0 1]);
        t=set(gca,'Ytick',[1:length(toPlot)]','YTickLabel',toPlotNames',...
            'FontSize',10,'Xtick',[1:length(toPlot)]','XTickLabel',toPlotNames','FontSize',14);
        t.Color='white';
        colorbar
        set(gcf,'units','points','position',[5,5,1000,1000])
    case 'AXES:allSNN_fit' % plot SSE as a function of clusters
        % OR plot upper and lower bounds of multi-task map
        K=[5:24,25];
        for p=1:length(K),
            toPlot{p}=sprintf('%dcluster',K(p));
        end
        
        % Aesthetics
        CAT.markersize=8;
        CAT.markertype='^';
        CAT.linewidth=4;
        CAT.linecolor={'r','k'};
        CAT.markercolor={'r','k'};
        CAT.markerfill={'r','k'};
        
        sc1_sc2_functionalAtlas('MAP:PLOT:SNN',{'SC12'},toPlot,K,'CAT',CAT);
        set(gcf,'units','points','position',[5,5,1000,1000])
        set(gca,'FontSize',18);
        ylabel('R2')
        xlabel('clusters')
        legend('R2','R2adj','Location','SouthEast')
    case 'AXES:map_similarity' % RS: resting state maps (Buckner, Cole)
        toPlot=varargin{1}; % {'Buckner_7Networks','Buckner_17Networks','Cole_10Networks'}
        toPlotNames=varargin{2}; % {'Buckner7',...}
        
        [RI]=sc1_sc2_functionalAtlas('MAP:compare',toPlot);
        
        imagesc_rectangle(RI,'YDir','reverse');
        caxis([0 1]);
        t=set(gca,'Ytick',[1:length(toPlot)]','YTickLabel',toPlotNames',...
            'FontSize',10,'Xtick',[1:length(toPlot)]','XTickLabel',toPlotNames','FontSize',18,'Color','white');
        colorbar
        set(gcf,'units','points','position',[5,5,1000,1000])
    case 'AXES:group_curves' % make separate graphs for 'lob10','Buckner_7Networks','Buckner_17Networks','Cole_10Networks','SC12_10cluster'
        toPlot=varargin{1}; % 'SC12_10cluster'
        plotName=varargin{2}; % 'Multi-Task:Upper'
        
        % Aesthetics
        CAT.markersize=8;
        CAT.markertype='^';
        CAT.linewidth=4;
        CAT.linecolor={'r','k'};
        CAT.markercolor={'r','k'};
        CAT.markerfill={'r','k'};
        CAT.errorcolor={'r','k'};
        
        sc1_sc2_functionalAtlas('EVAL:PLOT:CURVES',toPlot,4,'group',1,'unique','CAT',CAT);
        
        % Labelling
        set(gca,'YLim',[0 0.6],'XLim',[0 8],'FontSize',18,'xtick',[0 8],'XTickLabel',{'0','35'});
        xlabel('Spatial Distances (mm)');
        ylabel('Activity Correlation (R)');
        title(plotName);
        set(gcf,'units','points','position',[5,5,1000,1000])
        %         legend({'within parcels', 'between parcels'},'Location','SouthWest')
    case 'AXES:MT_UpperLower' % makes graph for upper and lower bounds of multi-task plot
        toPlot={'SC12_10cluster','SC2_10cluster'};
        toPlotName='Multi-Task Parcellation';
        
        % aesthetics
        CAT.markersize=8;
        CAT.markertype='^';
        CAT.linewidth=4;
        CAT.linecolor={'r','k'};
        CAT.markercolor={'r','k'};
        CAT.markerfill={'r','k'};
        CAT.errorcolor={'r','k'};
        
        % plot upper and lower within and between curves
        sc1_sc2_functionalAtlas('EVAL:PLOT:CURVES',toPlot{1},4,'group',1,'unique','CAT',CAT);
        hold on
        CAT.linestyle='--';
        CAT.linewidth={2,2};
        sc1_sc2_functionalAtlas('EVAL:PLOT:CURVES',toPlot{2},4,'group',1,'unique','CAT',CAT);
        hold off
        
        % Labelling
        set(gca,'YLim',[0 0.6],'XLim',[0 8],'FontSize',18,'xtick',[0 8],'XTickLabel',{'0','35'});
        xlabel('Spatial Distances (mm)');
        ylabel('Activity Correlation (R)');
        title(toPlotName)
        set(gcf,'units','points','position',[5,5,1000,1000])
        %         legend('within parcels', 'between parcels','Location','SouthWest')
    case 'AXES:diff_allMaps'  % make summary graph for diff curves for all maps
        toPlot=varargin{1}; % {'lob10','Buckner_7Networks','Cole_10Networks','SC12_10cluster','SC2_10cluster'}
        plotName=varargin{2}; % {'Lobular','Buckner7','Cole10','Multi-Task:Upper','Multi-Task:Lower'}
        evalNums=varargin{3}; % [4 4 4 4 4]
        
        % aesthetics
        CAT.markertype='none';
        CAT.errorwidth=1.5; 
        CAT.linestyle={'-','-','-','-','-'};
        CAT.errorcolor={'k','g','b','r','r'};
        CAT.linecolor={'k','g','b','r','r'};
        CAT.linewidth={3, 3, 3, 3, 3};
        
        sc1_sc2_functionalAtlas('EVAL:PLOT:DIFF',toPlot,evalNums,'group',1,'unique','CAT',CAT); % always take crossval + unique
        
        % Labelling
        set(gca,'YLim',[0 0.18],'XLim',[0 8],'FontSize',18,'xtick',[0 8],'XTickLabel',{'0','35'});
        xlabel('Spatial Distances (mm)');
        ylabel('Difference');
        set(gcf,'units','points','position',[5,5,1000,1000])
        legend(plotName,'Location','NorthWest')
    case 'AXES:diff_SNNMaps' % make summary graph for diff curves for all SNN maps
        toPlot=varargin{1};% {'SC12_5cluster','SC12_7cluster','SC12_9cluster','SC12_11cluster','SC12_13cluster','SC12_15cluster','SC12_17cluster'}
        plotName=varargin{2}; % {'5','7','9','11','13','15','17};
        evalNums=varargin{3}; % repmat([4],length(plotName),1)
        
        % aesthetics
        CAT.markersize=8;
        CAT.markertype={'^','v','s'};
        CAT.linewidth=4;
        CAT.errorcolor={'g','y','c','m'};
        CAT.linecolor={'k','k','k','m'};
        CAT.markercolor={'g','y','c','m'};
        CAT.markerfill={'g','y','c','m'};
        CAT.linestyle={'-','-','-','-'};
        CAT.linewidth={4, 4, 4, 2};
        
        sc1_sc2_functionalAtlas('EVAL:PLOT:DIFF',toPlot,evalNums,'group',1,'unique','CAT',CAT); % always take crossval + unique
        
        % Labelling
        set(gca,'YLim',[0 0.18],'XLim',[0 8],'FontSize',18,'xtick',[0 8],'XTickLabel',{'0','35'});;
        xlabel('Spatial Distances (mm)');
        ylabel('Difference');
        set(gcf,'units','points','position',[5,5,1000,1000])
        legend(plotName,'Location','NorthWest')
    case 'AXES:MT_indiv_curves' % make average curve graph for individual subjects (or average of individual subjects)
        
        % aesthetics
        CAT.markersize=8;
        CAT.markertype='^';
        CAT.linewidth=4;
        CAT.linecolor={'r','k'};
        CAT.markercolor={'r','k'};
        CAT.markerfill={'r','k'};
        
        sc1_sc2_functionalAtlas('EVAL:PLOT:CURVES','SC12_10cluster',4,'indiv',1,'unique','CAT',CAT,'sn',returnSubjs); % always take crossval + unique
        
        % Labelling
        set(gca,'YLim',[0 0.6],'XLim',[0 8],'FontSize',18);
        xlabel('Spatial Bins');
        ylabel('Activity Correlation (R)');
        title('Individual')
        set(gcf,'units','points','position',[5,5,1000,1000])
        legend('within parcels', 'between parcels','Location','SouthEast')
    case 'AXES:MT_group_indiv' % group versus indiv for multi-task map
        % aesthetics
        CAT.markersize=8;
        CAT.markertype='^';
        CAT.linewidth={2 1};
        CAT.errorcolor={'r','r'};
        CAT.linecolor={'r','r'};
        CAT.markercolor={'r','r'};
        CAT.markerfill={'r','r'};
        CAT.linestyle={'-','--'};
        
        sc1_sc2_functionalAtlas('EVAL:PLOT:DIFF',{'SC12_10cluster','SC2_10cluster'},[4 4],'group',1,'unique','CAT',CAT); % always take crossval + unique
        hold on
        CAT.errorcolor={'k','k'};
        CAT.linecolor={'k','k'};
        CAT.markercolor={'k','k'};
        CAT.markerfill={'k','k'};
        CAT.linestyle={'-','--'};
        CAT.linewidth={2 1};
        sc1_sc2_functionalAtlas('EVAL:PLOT:DIFF',{'SC12_10cluster','SC2_10cluster'},[4 4],'indiv',1,'unique','CAT',CAT,'sn',returnSubjs);
        
        % Labelling
        set(gca,'FontSize',18);
        xlabel('Spatial Bins');
        ylabel('Difference');
        set(gcf,'units','points','position',[5,5,1000,1000])
        %         legend('group','individual','Location','NorthWest')
    case 'AXES:boundary_strength' % makes separate graphs for 'lob10','Buckner_7Networks','Buckner_17Networks','Cole_10Networks','SC12_10cluster'
        toPlot=varargin{1}; % 'lob10','Buckner_7Networks' etc
        toPlotName=varargin{2}; % 'Lobular', 'Buckner7' etc
        
        sc1_sc2_functionalAtlas('STRENGTH:visualise_bound',toPlot)
        %         title(toPlotName);
        set(gca,'FontSize',18);
        set(gca,'Xticklabel',[],'Yticklabel',[])
        axis off
        set(gcf,'units','points','position',[5,5,1000,1000])
    case 'AXES:featSpace'  % graph the feature loadings for each network
        toPlot=varargin{1}; % 'winner' or 'all'
        
        sc1_sc2_functionalAtlas('ENCODE:project_featSpace','SC12_10cluster',toPlot)
        
        set(gca,'FontSize',9);
        set(gcf,'units','points','position',[5,5,1000,1000])
    case 'AXES:taskSpace'  % graph the task loadings for each network
        toPlot=varargin{1}; % either [4,10] or 'all' (if you want all taskLoadings
        
        % Aesthetics
        CAT.markersize=12;
        CAT.labelsize=18;
        sc1_sc2_functionalAtlas('ENCODE:project_taskSpace','SC12_10cluster',toPlot,'CAT',CAT)
        
        set(gcf,'units','points','position',[5,5,1000,1000])
        set(gca,'FontSize',9);
    case 'AXES:featMatrix' % graph the feature matrix
        sc1_sc2_functionalAtlas('ENCODE:project_featSpace','SC12_10cluster','featMatrix')
        
        set(gcf,'units','points','position',[5,5,1000,1000])
        set(gca,'FontSize',9);
    case 'AXES:varianceFreq'
        % Aesthetics
        CAT.markersize=8;
        CAT.markertype='^';
        CAT.linewidth=4;
        CAT.linecolor={'r','k'};
        CAT.markercolor={'r','k'};
        CAT.markerfill={'r','k'};
        CAT.errorcolor={'r','k'};
        
        sc1_sc2_functionalAtlas('PLOT:spatialFreqCorr',CAT)
        set(gca,'FontSize',18);
        set(gcf,'units','points','position',[5,5,1000,1000])
    case 'AXES:interSubjCorr'
        % Aesthetics
        CAT.markersize=8;
        CAT.markertype='^';
        CAT.linewidth=4;
        CAT.linecolor={'r','k'};
        CAT.markercolor={'r','k'};
        CAT.markerfill={'r','k'};
        CAT.errorcolor={'r','k'};
        
        sc1_sc2_functionalAtlas('PLOT:interSubjCorr',CAT)
        set(gca,'FontSize',18);
        set(gcf,'units','points','position',[5,5,1000,1000])
        legend({'within','between'},'Location','SouthWest')
        
    case 'FIGURE1' % Motor Feature Model
    case 'FIGURE2' % Representational Structure
        sc1_sc2_functionalAtlas('AXES:MDS','all')
    case 'FIGURE3' % Lobular versus Functional
        subplot(2,2,1)
        sc1_sc2_functionalAtlas('AXES:group_curves','lob10','Lobular Parcellation')
        subplot(2,2,2)
        sc1_sc2_functionalAtlas('AXES:MT_UpperLower','Multi-Task Parcellation')
        subplot(2,2,3)
        sc1_sc2_functionalAtlas('AXES:boundary_strength','lob10','Lobular Parcellation')
        subplot(2,2,4)
        sc1_sc2_functionalAtlas('AXES:boundary_strength','SC12_10cluster','Multi-Task Parcellation')
    case 'FIGURE5' % Resting State Parcellations
        subplot(2,3,1)
        sc1_sc2_functionalAtlas('AXES:group_curves','Buckner_7Networks','Buckner7')
        subplot(2,3,2)
        sc1_sc2_functionalAtlas('AXES:group_curves','Buckner_17Networks','Buckner17')
        subplot(2,3,3)
        sc1_sc2_functionalAtlas('AXES:group_curves','Cole_10Networks','Cole10')
        subplot(2,3,4)
        sc1_sc2_functionalAtlas('AXES:boundary_strength','Buckner_7Networks','Buckner7')
        subplot(2,3,5)
        sc1_sc2_functionalAtlas('AXES:boundary_strength','Buckner_17Networks','Buckner17')
        subplot(2,3,6)
        sc1_sc2_functionalAtlas('AXES:boundary_strength','Cole_10Networks','Cole10')
    case 'FIGURE6' % Summary Graph
        sc1_sc2_functionalAtlas('AXES:diff_allMaps',{'lob10','Cole_10Networks','Buckner_7Networks','SC12_10cluster','SC2_10cluster'},...
            {'Lobular','Cole10','Buckner7','Multi-Task:Upper','Multi-Task:Lower'},...
            [4,4,4,4,4])
    case 'FIGURE7' % Multiple Multi-Task Parcellations
        subplot(2,2,[3 4])
        sc1_sc2_functionalAtlas('AXES:diff_SNNMaps',{'SC12_7cluster','SC12_10cluster','SC12_17cluster','lob10'},{'Multi-Task:7','Multi-Task:10','Multi-Task:17','Lobular'},repmat([4],4,1))
        %         subplot(2,2,2)
        %         sc1_sc2_functionalAtlas('AXES:boundary_strength','SC12_17cluster','Multi-Task17')
        %         subplot(2,2,1)
        %         sc1_sc2_functionalAtlas('AXES:boundary_strength','SC12_7cluster','Multi-Task7')
    case 'FIGURE8a'% Multi-Task Map
        sc1_sc2_functionalAtlas('AXES:map','SC12_10cluster','Multi-Task')
    case 'FIGURE8b'% Graph the "winner" features for each cluster
        sc1_sc2_functionalAtlas('AXES:featSpace','winner')
    case 'FIGURE8c'% Assigning Semantic Labels
        figure(1)
        sc1_sc2_functionalAtlas('AXES:taskSpace',[4,10])
        figure(2)
        sc1_sc2_functionalAtlas('AXES:taskSpace',[3,6])
    case 'FIGURE9' % Group Versus Individual
        %         subplot(2,4,1)
        %         sc1_sc2_functionalAtlas('AXES:MT_indiv_curves',2,'S02')
        %         subplot(2,4,3)
        %         sc1_sc2_functionalAtlas('AXES:MT_indiv_curves',4,'S04')
        %         subplot(2,4,[5:8])
        %         sc1_sc2_functionalAtlas('AXES:MT_group_indiv')
        %         subplot(2,4,2)
        %         sc1_sc2_functionalAtlas('AXES:map','SC12_10cluster','S02','sn',2)
        %         subplot(2,4,4)
        %         sc1_sc2_functionalAtlas('AXES:map','SC12_10cluster','S04','sn',4)
        subplot(2,2,[1,2])
        sc1_sc2_functionalAtlas('AXES:interSubjCorr')
        subplot(2,2,[3,4])
        sc1_sc2_functionalAtlas('AXES:MT_group_indiv')
        
    case 'SUPP1'   % Map Similarity
        sc1_sc2_functionalAtlas('AXES:map_similarity',{'lob10','SC12_10cluster','Buckner_7Networks','Buckner_17Networks','Cole_10Networks'},...
            {'Lobular','Multi-Task','Buckner7','Buckner17','Cole10'});
    case 'SUPP2'   % RDM
        sc1_sc2_functionalAtlas('AXES:RDM','all')
    case 'SUPP3'   % Reliability of activity
        subplot(2,1,1)
        sc1_sc2_functionalAtlas('AXES:reliability','A','Activity Reliability')
        subplot(2,1,2)
        sc1_sc2_functionalAtlas('AXES:reliability','D','Distance Reliability')
    case 'SUPP4'   % Feature & Task Loading Matrices
        subplot(2,2,1)
        sc1_sc2_functionalAtlas('AXES:featSpace','all')
        subplot(2,2,2)
        sc1_sc2_functionalAtlas('AXES:taskSpace','all')
        subplot(2,2,[3,4])
        sc1_sc2_functionalAtlas('AXES:featMatrix')
end
% Local functions
function dircheck(dir)
if ~exist(dir,'dir');
    warning('%s doesn''t exist. Creating one now. You''re welcome! \n',dir);
    mkdir(dir);
end

% InterSubj Corr
function C=intersubj_corr(Y)
numSubj=size(Y,3);
for i=1:numSubj
    for j=1:numSubj
        C(i,j)=nansum(nansum(Y(:,:,i).*Y(:,:,j)))/...
            sqrt(nansum(nansum(Y(:,:,i).*Y(:,:,i)))*...
            nansum(nansum(Y(:,:,j).*Y(:,:,j))));
    end;
end

function E = expectedObs(pivotTable)
%Number of observations
N = sum(nansum(pivotTable));

%Number of observations marginals
X = nansum(pivotTable,1);
Y = nansum(pivotTable,2);

%Propotions
Px = X/N;
Py = Y/N;

%Expected number of observations
E = Py * Px * N;

function G = Gtest(map1,map2, validIdx)
%Argument controller
if nargin == 2
    %Valid vertices index vector
    validIdx = ones(length(map1),1);
end
%Observations
O = pivottable(map1, map2, map1, 'length', 'subset', validIdx);
E = expectedObs(O);
%G-test
G = 2 * sum(nansum( O .* (log(O) - log(E)) ));






