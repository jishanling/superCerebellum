function varargout=sc1_sc2_functionalAtlas(what,varargin)

% Directories
baseDir          = '/Users/maedbhking/Documents/Cerebellum_Cognition';
% baseDir            = '/Volumes/MotorControl/data/super_cerebellum_new';
% baseDir          = '/Users/jdiedrichsen/Data/super_cerebellum_new';

studyDir{1}     =fullfile(baseDir,'sc1');
studyDir{2}     =fullfile(baseDir,'sc2');
suitDir         ='/suit';
caretDir        ='/surfaceCaret';
regDir          ='/RegionOfInterest/';
encodeDir       ='/encoding';

subj_name = {'s01','s02','s03','s04','s05','s06','s08','s09','s10','s11',...
    's12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24',...
    's25','s26','s27','s28','s29','s30','s31'};

returnSubjs=[2,3,4,6,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31];

switch what
    
    case 'TASKSPACE:overlap'
        
    case 'ACTIVITY:get_data'
        % load in allSubjs data struct from sc1 & sc2
        H=[];
        for study=1:2,
            load(fullfile(studyDir{study},encodeDir,'glm4','cereb_avrgDataStruct.mat'));
            T.studyNum=repmat(study,size(T.SN,1),1);
            H=addstruct(H,T);
            clear T
        end
        
        sn=unique(H.SN);
        
        % get session average for each study separately
        for s=1:length(sn),
            idx=1;
            for study=1:2,
                numConds=length(unique(H.cond(H.studyNum==study)));
                for c=1:numConds, % get average across sessions
                    indx = H.cond==c & H.studyNum==study & H.SN==(sn(s));
                    avrgData(idx,:)=nanmean(H.data(indx,:),1);
                    idx=idx+1;
                    clear indx
                end
                fprintf('subj%d averaged sessions for study%d \n',sn(s),study)
            end
            % subtract condition avrg baseline (across 61-1 conditions)
            all=[1:size(avrgData,1)];
            for c=1:size(avrgData,1),
                X=nanmean(avrgData(all(all~=c),:),1);
                data(c,:,s)=avrgData(c,:)-X;
            end
            clear avrgData
            fprintf('subj%d new baseline \n',sn(s))
        end
        
        varargout={data,sn,V,volIndx};
    case 'ACTIVITY:get_motorFeats'
        D.saccades=[28;28;98;31;31;37;37;29;29;46;38;38;36;41;38;38;20;20;25;25;50;22;22;26;26;32;32;32;45;...
            60;15;15;15;26;26;50;16;30;35;71;60;58;38;38;43;43;43;25;25;30;30;23;23;50;40;40;40;98;21;21;45];
        D.duration = [15;15;30;15;15;15;15;15;15;30;15;15;30;30;15;15;15;15;15;15;30;15;15;15;15;10;10;10;30;...
            30;10;10;10;15;15;30;10;10;10;30;30;30;15;15;10;10;10;15;15;15;15;10;10;10;10;10;10;30;15;15;30];
        D.lHand = [0;15;2;0;0;8;7;0;0;0;0;0;0;0;12;12;0;7;0;0;0;3.5;4;0;0;5;5;5;0;...
            2;2;2;2;0;0;0;1;1;1;0;0;0;12;12;0;0;0;0;0;0;0;1;1;1;5;5;5;2;0;0;0];
        D.rHand = [0;0;0;0;0;0;0;5;5;0;8;7;15;0;12;12;0;0;0;7;0;3.5;4;0;0;0;0;0;0;...
            2;0;0;0;0;0;0;1;1;1;0;0;0;12;12;3;3;3;0;7;5;5;1;1;1;0;0;0;0;0;0;0];
        featNames={'lHand','rHand','saccades'}; 
        varargout={D,featNames};
    case 'ACTIVITY:taskConds'
        stats=varargin{1}; % 'yes' or 'no'
        
        pThresh=.05;
        
        [data,sn,V,volIndx]=sc1_sc2_functionalAtlas('ACTIVITY:get_data');
        numTasks=size(data,1); 
        
        [D,featNames]=sc1_sc2_functionalAtlas('ACTIVITY:get_motorFeats');

        T=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        condNames=[T.condNames;featNames']; 
        
        X=[eye(numTasks) D.lHand./D.duration D.rHand./D.duration D.saccades./D.duration];
        numTasks=size(X,2); 
        
        % set up volume info
        Yy=zeros(numTasks,length(sn),V.dim(1)*V.dim(2)*V.dim(3));
        C{1}.dim=V.dim;
        C{1}.mat=V.mat;
        
        % do ridge regression
        for s=1:length(sn),
            X1=bsxfun(@minus,X,mean(X));
            X1=bsxfun(@rdivide,X1,sum(X1.^2));
            B=(X1'*X1+eye(numTasks)*.5)\(X1'*data(:,:,s));
            % make volume
            Yy(:,s,volIndx)=B;
            clear X1 B
            fprintf('subj%d done \n',sn(s))
        end;
        
        switch stats,
            case 'yes'
                Yy(Yy==0)=nan;
                for ii=1:numTasks,
                    data=ssqrt(Yy(ii,:,:));
                    data=reshape(data,[size(data,2),size(data,3)]);
                    cStats{ii}= caret_getcSPM('onesample_t','data',data(1:length(sn),:)');
                    P=spm_P_FDR(cStats{ii}.con.Z,cStats{ii}.con.df,'Z',1,sort(cStats{ii}.con.Z_P,'ascend'));
                    %                     P=spm_P_Bonf(cStats{ii}.con.Z,cStats{ii}.con.df,'Z',size(cStats{ii}.data,1),1);
                    c=cStats{ii}.con.Z;
                    c(P>pThresh)=nan;
                    indices(ii,:)=c;
                    clear c
                end
            case 'no'
                % if numSubjs > 1 get avg
                Yy=permute(Yy,[2 1 3]);
                indices=nanmean(Yy,1);
                indices=reshape(indices,[size(indices,2),size(indices,3)]);
        end
        clear data
        
        % map vol2surf
        indices=reshape(indices,[size(indices,1) V.dim(1),V.dim(2),V.dim(3)]);
        for i=1:size(indices,1),
            data=reshape(indices(i,:,:,:),[C{1}.dim]);
            C{i}.dat=data;
        end
        M=caret_suit_map2surf(C,'space','SUIT','stats','nanmean','column_names',condNames);  % MK created caret_suit_map2surf to allow for output to be used as input to caret_save
        
        % save out metric file
        if strcmp(stats,'yes'),
            save(fullfile(studyDir{2},caretDir,'suit_flat','glm4','FDRCorr_taskConds.metric'),M);
        else
            caret_save(fullfile(studyDir{2},caretDir,'suit_flat','glm4','unCorr_taskConds.metric'),M);
        end
    case 'ACTIVITY:motorFeats'
        stats=varargin{1}; % 'yes' or 'no'
        
        pThresh=.05;
        
        [data,sn,V,volIndx]=sc1_sc2_functionalAtlas('ACTIVITY:get_data');
        
        [D,featNames]=sc1_sc2_functionalAtlas('ACTIVITY:get_motorFeats');
        
        X=[D.lHand./D.duration D.rHand./D.duration D.saccades./D.duration];
        
        % set up volume info
        numFeat=size(X,2);
        Yy=zeros(numFeat,length(sn),V.dim(1)*V.dim(2)*V.dim(3));
        C{1}.dim=V.dim;
        C{1}.mat=V.mat;
        
        % do ridge regression
        for s=1:length(sn),
            X1=bsxfun(@minus,X,mean(X));
            X1=bsxfun(@rdivide,X1,sum(X1.^2));
            B=(X1'*X1+eye(numFeat)*.5)\(X1'*data(:,:,s));
            % make volume
            Yy(:,s,volIndx)=B;
            clear X1 B
            fprintf('subj%d done \n',sn(s))
        end;
        
        switch stats,
            case 'yes'
                Yy(Yy==0)=nan;
                for ii=1:numFeat,
                    data=ssqrt(Yy(ii,:,:));
                    data=reshape(data,[size(data,2),size(data,3)]);
                    cStats{ii}= caret_getcSPM('onesample_t','data',data(1:length(sn),:)');
                    P=spm_P_FDR(cStats{ii}.con.Z,cStats{ii}.con.df,'Z',1,sort(cStats{ii}.con.Z_P,'ascend'));
                    %                     P=spm_P_Bonf(cStats{ii}.con.Z,cStats{ii}.con.df,'Z',size(cStats{ii}.data,1),1);
                    c=cStats{ii}.con.Z;
                    c(P>pThresh)=nan;
                    indices(ii,:)=c;
                    clear c
                end
            case 'no'
                % if numSubjs > 1 get avg
                Yy=permute(Yy,[2 1 3]);
                indices=nanmean(Yy,1);
                indices=reshape(indices,[size(indices,2),size(indices,3)]);
        end
        clear data
        
        % map vol2surf
        indices=reshape(indices,[size(indices,1) V.dim(1),V.dim(2),V.dim(3)]);
        for i=1:size(indices,1),
            data=reshape(indices(i,:,:,:),[C{1}.dim]);
            C{i}.dat=data;
        end
        M=caret_suit_map2surf(C,'space','SUIT','stats','nanmean','column_names',featNames');  % MK created caret_suit_map2surf to allow for output to be used as input to caret_save
        
        % visualise motor features
        for m=1:numFeat,
            figure()
            suit_plotflatmap(M.data(:,m))
        end
        
        % save out metric file
        if strcmp(stats,'yes'),
            save(fullfile(studyDir{2},caretDir,'suit_flat','glm4','FDRCorr_motorFeats.metric'),M);
        else
            caret_save(fullfile(studyDir{2},caretDir,'suit_flat','glm4','unCorr_motorFeats.metric'),M);
        end
    case 'ACTIVITY:reliability'
        glm=varargin{1};
        type=varargin{2}; % 'cerebellum' or 'cortex'
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        load(fullfile(studyDir{2},encodeDir,sprintf('glm%d',glm),sprintf('allVox_sc1_sc2_sess_%s.mat',type)));
        Y=Yy;clear Yy;
        
        numSubj=length(Y);
        
        for subj=1:numSubj,
            idx=1;
            for study=1:2,
                
                D1=getrow(D,D.StudyNum==study);
                
                cN=condNum{study}-1;  % Important: In the allVox file, instruction is still included!
                pN=partNum{study};    % Partition Numner
                sN=(pN>8)+1;          % Sessions Number
                for se=1:2,
                    X1=indicatorMatrix('identity_p',cN.*(sN==se));  % This one is the matrix that related trials-> condition numbers
                    X2=indicatorMatrix('identity_p',D1.condNumUni.*D1.overlap); % THis goes from condNum to shared condNumUni
                    Yf(:,:,idx,subj)=pinv(X1*X2)*Y{subj}{study};
                    Yf(:,:,idx,subj)=bsxfun(@minus,Yf(:,:,idx,subj),mean(Yf(:,:,idx,subj)));
                    idx=idx+1;
                end;
            end;
            CORR(:,:,subj)=interSess_corr(Yf(:,:,:,subj));
            T.SN(subj,1)      = subj;
            T.within1(subj,1) = CORR(1,2,subj);
            T.within2(subj,1) = CORR(3,4,subj);
            T.across(subj,1)  = mean(mean(CORR(1:2,3:4,subj)));
        end
        save(fullfile(studyDir{2},regDir,'glm4','patternReliability.mat'),'T','CORR')
    case 'FIGURES:reliabilityA'
        % load relability 
        load(fullfile(studyDir{2},regDir,'glm4','patternReliability.mat')); 
        
        % within & between-subj reliability
        myboxplot([],[T.within1 T.within2 T.across],'style_tukey');
        drawline(0,'dir','horz');
        set(gca,'XTickLabel',{'SC1','SC2','across'});
        set(gcf,'PaperPosition',[2 2 4 3]);
        wysiwyg;
        ttest(sqrt(T.within1.*T.within2),T.across,2,'paired');
        
        % group reliability
        X=nanmean(CORR,3);
        fprintf('group reliability for study1 is %2.2f \n',X(1,2));
        fprintf('group reliability for study2 is %2.2f \n',X(3,4));
        fprintf('group reliability for shared tasks is %2.2f \n',mean(mean(X(1:2,3:4))));   

    case 'REPRESENTATION:get_distances'
        type=varargin{1}; % 'cortex' or 'cerebellum'
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
        F=sc1_sc2_functionalAtlas('ACTIVITY:get_motorFeats');
        
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
        type=varargin{2}; % 'cerebellum' or 'cortex'
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
            for j=1:numSess
                dist(:,j  ,i)  = ssqrt(diag(C*G_hat_sc1(i1,i1,j,i)*C'));
                dist(:,j+2,i)  = ssqrt(diag(C*G_hat_sc2(i2,i2,j,i)*C'));
            end;
            CORR(:,:,i)    = corr(dist(:,:,i));
            T.SN(i,1)      = i;
            T.within1(i,1) = CORR(1,2,i);
            T.within2(i,1) = CORR(3,4,i);
            T.across(i,1)  = mean(mean(CORR(1:2,3:4,i)));
        end;

        save(fullfile(studyDir{2},regDir,'glm4','distanceReliability.mat'),'T','dist')
    case 'FIGURES:RDM'
        % load in fullRDM
        [fullRDM,T,~]=sc1_sc2_functionalAtlas('REPRESENTATION:get_distances','cerebellum','all'); 

        % Plot RDM
        reOrder=[1,2,6,7,8,9,10,11,12,13,14,17,18,22,23,3,4,5,15,16,19,20,21,24,25,26,...
            27,28,29,58,59,60,43,44,49,48,36,34,35,55,56,57,61,30,31,32,33,37,38,39,40,41,42,45,46,47,50,51,52,53,54]'; % reorder
        figure()
        avrgFullRDM=ssqrt(nanmean(fullRDM,3));
        numDist=size(avrgFullRDM,1);
        imagesc_rectangle(avrgFullRDM(reOrder,reOrder),'YDir','reverse');
        caxis([0 1]);
        t=set(gca,'Ytick',[1:numDist]','YTickLabel',T.condNames(reOrder)','FontSize',12,'FontWeight','bold');
        t.Color='white';
        colorbar   
    case 'FIGURES:MDS'
         % load in fullRDM
        [fullRDM,T,X]=sc1_sc2_functionalAtlas('REPRESENTATION:get_distances','cerebellum','none'); 
        
        avrgFullRDM=ssqrt(nanmean(fullRDM,3));

        vecRDM = rsa.rdm.vectorizeRDM(avrgFullRDM);
        [Y,~] = rsa_classicalMDS(vecRDM,'mode','RDM');
        B = (X'*X+eye(size(X,2))*0.0001)\(X'*Y); % ridge regression
        Yr    = Y  - X(:,1:3)*B(1:3,:); % remove motor features
        
        clustTree = linkage(Yr,'average');
        indx = cluster(clustTree,'cutoff',1);
        
        % define cluster colour
        colour={[1 0 0],[0 1 0],[0 0 1],[0.3 0.3 0.3],[1 0 1],[1 1 0],[0 1 1],[0.5 0 0.5],[0.8 0.8 0.8],[.07 .48 .84],[.99 .76 .21],[.11 .7 .68],[.39 .74 .52],[.21 .21 .62],[0.2 0.2 0.2],[.6 .6 .6],[.3 0 .8],[.8 0 .4],[0 .9 .2],[.1 .3 0],[.2 .4 0],[.63 0 .25],[0 .43 .21],[.4 0 .8]};
        [V,L]   = eig(Yr*Yr');
        [l,i]   = sort(diag(L),1,'descend'); % Sort the eigenvalues
        V       = V(:,i);
        X       = bsxfun(@times,V,sqrt(l'));
        X = real(X);
        
        for study=1:2,
            % which condNames and numbers are we taking ?
            if study==2,
                condIndx=T.condNum(T.StudyNum==study)+length(T.condNum(T.StudyNum==1)); 
            else
                condIndx=T.condNum(T.StudyNum==study); 
            end
            indxShort=indx(T.condNum(T.StudyNum==study)); 
            CAT.markercolor= {colour{indxShort}};
            CAT.markerfill = {colour{indxShort}};
            CAT.markersize = 10;
            
            X1=X(condIndx,condIndx);
            figure(study)
            scatterplot3(X1(:,1),X1(:,2),X1(:,3),'split',condIndx,'CAT',CAT,'label',T.condNames(T.StudyNum==study),'markertype','o');
            set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'Box','on');
            hold on;
            plot3(0,0,0,'+');
            % Draw connecting lines
%             for i=1:15,
%                 ind=clustTree(i,1:2);
%                 X(end+1,:)=(X(ind(1),:)+X(ind(2),:))/2;
%                 line(X(ind,1),X(ind,2),X(ind,3));
%             end;
            hold off;
            set(gcf,'PaperPosition',[2 2 8 8]);
            axis equal;
            set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'Box','on');
            wysiwyg;
            view([81 9]);
            clear X1  indxShort
        end     
    case 'FIGURES:reliabilityD'
        % load relability 
        load(fullfile(studyDir{2},regDir,'glm4','distanceReliability.mat')); 
        
        % within & between-subj reliability
        myboxplot([],[T.within1 T.within2 T.across],'style_tukey');
        set(gca,'XTickLabel',{'SC1','SC2','across'});
        set(gcf,'PaperPosition',[2 2 4 3]);
        wysiwyg;
        ttest(sqrt(T.within1.*T.within2),T.across,2,'paired');
        
        % group reliability
        X=nanmean(dist,3);
        groupCorr=corr(X);
        fprintf('group reliability for study1 is %2.2f \n',groupCorr(1,2));
        fprintf('group reliability for study2 is %2.2f \n',groupCorr(3,4));
        fprintf('group reliability for shared tasks is %2.2f \n',mean(mean(groupCorr(1:2,3:4))));    
    
    case 'ATLAS:finalMap'
        % example: 'sc1_sc2_functionalAtlas('ATLAS:finalMap',[2],1,13,'leaveOneOut')
        K=varargin{1}; % number of clusters
        metric=varargin{2}; % 'yes' or 'no'
        
        sn=returnSubjs;
        dircheck(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_atlasFinal%d',K)));
        outName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_atlasFinal%d',K),'map.nii');
        
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',sn,[1,2],'build');
        [F,G,Err1]=semiNonNegMatFac(X_C,K,'threshold',0.01);
        save(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_atlasFinal%d',K),'parcels.mat'),'F','G','volIndx','V');
        fprintf('final atlas created for %s \n',K)
        
        [~,groupFeat]=max(G,[],2);
        
        % map features on group
        Yy=zeros(1,V.dim(1)*V.dim(2)*V.dim(3));
        Yy(1,volIndx)=groupFeat;
        Yy=reshape(Yy,[V.dim(1),V.dim(2),V.dim(3)]);
        Yy(Yy==0)=NaN;
        Vv{1}.dat=Yy;
        Vv{1}.dim=V.dim;
        Vv{1}.mat=V.mat;
        
        % save out vol of feats
        exampleVol=fullfile(studyDir{2},suitDir,'glm4','s02','wdbeta_0001.nii'); % must be better way ??
        X=spm_vol(exampleVol);
        X.fname=outName;
        X.private.dat.fname=V.fname;
        spm_write_vol(X,Vv{1}.dat);
        
        % save out metric file
        % map vol 2 surf
        M=caret_suit_map2surf(Vv,'space','SUIT','stats','mode','type','paint');
        
        switch metric,
            case 'yes'
                M.column_name={'atlasFinal'};
                caret_save(fullfile(studyDir{2},caretDir,'suit_flat','glm4','atlasFinal.paint'),M);
                % make areacolour
                cmap=load(fullfile(studyDir{2},caretDir,'suit_flat','glm4','features.txt'));
                M.column_name=M.paintnames;
                M.column_color_mapping=repmat([-5 5],length(M.column_name),1);
                M.data=cmap(1:length(M.column_name),2:4);
                caret_save(fullfile(studyDir{2},caretDir,'suit_flat','glm4','atlasFinal.areacolor'),M);
            case 'no'
                figure()
                title('final atlas')
                M=caret_suit_map2surf(Vv,'space','SUIT','stats','mode');
                suit_plotflatmap(M.data,'type','label','cmap',colorcube(K))
        end
        
    case 'EVAL:get_data'
        sn=varargin{1}; % Subj numbers to include
        study=varargin{2}; % 1 or 2 or [1,2]
        type=varargin{3}; % 'build' or 'eval'. For build - we get group data. For eval - we get indiv data
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
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
    case 'EVAL:unCrossval:GROUP'
        % code is written now so that func map is built on sc1+sc2 (allConds) and
        % evaluated on sc1+sc2 (allConds)
        mapType=varargin{1}; % options are 'lob10','lob26','bucknerRest', or 'atlasFinal<num>'
        eval=varargin{2}; % [1] [2] or [1,2]
        
        if size(eval,2)==2,
            outName='spatialBoundfunc3.mat';
            outDir=fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),outName);
        else
            outName=sprintf('spatialBoundfunc%d.mat',eval);
            outDir=fullfile(studyDir{eval},'encoding','glm4',sprintf('groupEval_%s',mapType),outName);
        end
        
        % load in mapOR
        mapName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'map.nii');
        
        % load in func data (sc1+sc2) to test
        [X_C,volIndx,V,sn] = sc1_sc2_functionalAtlas('EVAL:get_data',returnSubjs,eval,'eval');
        
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
        % Now calculate the uncrossvalidated estimation of the correlation for each subject
        for s=1:length(sn),
            B=X_C(:,:,s);
            R.SN = ones(length(R.N),1)*sn(s);
            fprintf('%d correl \n',sn(s));
            R.corr=mva_spatialCorr(B,BIN);
            R.crossval = zeros(length(R.corr),1);
            RR = addstruct(RR,R);
            fprintf('uncrossval eval done for subj %d \n',sn(s));
        end;
        save(outDir,'-struct','RR');
    case 'EVAL:crossval:MAKE'
        % example: 'sc1_sc2_functionalAtlas('EVAL:make_SNN_map',[2],1,13,'leaveOneOut')
        sn=varargin{1}; % subjNum or 'group'
        study=varargin{2}; % 1 or 2
        K=varargin{3}; % number of clusters
        type=varargin{4}; % 'group','indiv','leaveOneOut'
        
        % figure out if individual or group or leave one out
        switch type,
            case 'group'
                sn=returnSubjs;
                outDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K));dircheck(outDir);
                outName=fullfile(outDir,'SNN.mat');
            case 'indiv'
                outDir=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn});
                outName=fullfile(outDir,sprintf('SNN_SC%d_%dcluster.mat',study,K));
            case 'leaveOneOut'
                outDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K)); dircheck(outDir);
                outName=fullfile(outDir,sprintf('SNN_leaveOut_%s.mat',subj_name{sn}));
                sn=returnSubjs(returnSubjs~=sn);
        end
        
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',sn,study,'build');
        [F,G,Err1]=semiNonNegMatFac(X_C,K,'threshold',0.01);
        fprintf('SNN map (%d clusters) created for %s \n',K,type)
        save(fullfile(outName),'F','G','volIndx','V');
    case 'EVAL:crossval:VISUALISE'
        sn=varargin{1}; % [2] or 'group'
        study=varargin{2}; % 1 or 2 or [1,2]
        K=varargin{3}; % number of clusters
        type=varargin{4}; % 'group','indiv',or 'leaveOneOut'
        
        % figure out if individual or group
        switch type,
            case 'group'
                outName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K),'map.nii');
                load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K),'SNN.mat'));
            case 'indiv'
                outName=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('map_SC%d_%dcluster.nii',study,K));
                load(fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('SNN_SC%d_%dcluster.mat',study,K)));
            case 'leaveOneOut'
                outName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K),sprintf('map_leaveOut_%s.nii',subj_name{sn}));
                load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K),sprintf('SNN_leaveOut_%s.mat',subj_name{sn})));
        end
        
        [~,groupFeat]=max(G,[],2);
        
        % map features on group
        Yy=zeros(1,V.dim(1)*V.dim(2)*V.dim(3));
        Yy(1,volIndx)=groupFeat;
        Yy=reshape(Yy,[V.dim(1),V.dim(2),V.dim(3)]);
        Yy(Yy==0)=NaN;
        Vv{1}.dat=Yy;
        Vv{1}.dim=V.dim;
        Vv{1}.mat=V.mat;
        
        % save out vol of SNN feats
        exampleVol=fullfile(studyDir{study},suitDir,'glm4','s02','wdbeta_0001.nii'); % must be better way ??
        X=spm_vol(exampleVol);
        X.fname=outName;
        X.private.dat.fname=V.fname;
        spm_write_vol(X,Vv{1}.dat);
        
        figure()
        title(sprintf('%d clusters',K))
        M=caret_suit_map2surf(Vv,'space','SUIT','stats','mode');
        suit_plotflatmap(M.data,'type','label','cmap',colorcube(K))
    case 'EVAL:crossval:GROUP'
        % should be that func map from sc1 is always evaluated
        % on sc2 data (and vice versa)
        mapType=varargin{1}; % options are 'lob10','lob26','bucknerRest','SC<studyNum>_<num>cluster', or 'SC<studyNum>_POV<num>'
        data=varargin{2}; % evaluating data from study [1] or [2] ?
        condType=varargin{3}; % 'unique' or 'all'. Are we evaluating on all taskConds or just those unique to either sc1 or sc2 ?
        
        % load in map
        mapName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'map.nii');
        
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
        for sn=unique(T.SN)';
            for c=1:length(idx),
                i1(c) = find(T.SN==sn & T.sess==1 & T.cond==idx(c));
                i2(c) = find(T.SN==sn & T.sess==2 & T.cond==idx(c));
            end
            B=(T.data(i1,voxIn)+T.data(i2,voxIn))/2;
            %             fprintf('%d cross\n',sn);
            R.SN = ones(length(R.N),1)*sn;
            R.corr = mva_spatialCorr(T.data([i1;i2],voxIn),BIN,...
                'CrossvalPart',T.sess([i1;i2],1),'excludeNegVoxels',1);
            R.crossval = ones(length(R.corr),1);
            RR = addstruct(RR,R);
            fprintf('%d correl\n',sn);
            R.corr=mva_spatialCorr(B,BIN);
            R.crossval = zeros(length(R.corr),1);
            RR = addstruct(RR,R);
        end;
        save(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_%s.mat',data,condType)),'-struct','RR');
    case 'EVAL:crossval:INDIV' % NEED TO UPDATE FOR SNN !!
        sn=varargin{1}; % [2]
        study=varargin{2}; % 1 or 2
        var=varargin{3}; % .95
        data=varargin{4}; % 'func1' or 'func2'
        
        outName=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('SC%d_%dPOV_spatialBound%s.mat',study,var*100,data));
        mapName=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('map_SC%d_%dPOV.nii',study,var*100));
        
        % load in func data to test
        switch data
            case 'func1'
                load(fullfile(studyDir{1},'encoding','glm4','cereb_avrgDataStruct.mat'));
            case 'func2'
                load(fullfile(studyDir{2},'encoding','glm4','cereb_avrgDataStruct.mat'));
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
        % Now calucalte the uncrossvalidated and crossvalidated
        % Estimation of the correlation for each subject
        for s=unique(T.SN)';
            i1 = find(T.SN==s & T.sess==1);
            i2 = find(T.SN==s & T.sess==2);
            D=(T.data(i1,voxIn)+T.data(i2,voxIn))/2; % why divide by 2 ?
            
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
    case 'EVAL:crossval:LEAVEONEOUT' % NEED TO UPDATE FOR SNN !!
        sn=varargin{1}; % [2]
        study=varargin{2}; % 1 or 2
        var=varargin{3}; % .95
        data=varargin{4}; % 'func1' or 'func2
        
        outName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100),sprintf('spatialBound%s_eval_%s.mat',data,subj_name{sn}));
        mapName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100),sprintf('map_leaveOut_%s.nii',subj_name{sn}));
        
        % load in func data to test
        switch data
            case 'func1'
                load(fullfile(studyDir{1},'encoding','glm4','cereb_avrgDataStruct.mat'));
            case 'func2'
                load(fullfile(studyDir{2},'encoding','glm4','cereb_avrgDataStruct.mat'));
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
        % Now calucalte the uncrossvalidated and crossvalidated
        % Estimation of the correlation for each subject
        i1 = find(T.SN==sn & T.sess==1);
        i2 = find(T.SN==sn & T.sess==2);
        D=(T.data(i1,voxIn)+T.data(i2,voxIn))/2; % why divide by 2 ?
        
        fprintf('%d cross\n',sn);
        R.SN = ones(length(R.N),1)*sn;
        R.corr = mva_spatialCorr(T.data([i1;i2],voxIn),BIN,...
            'CrossvalPart',T.sess([i1;i2],1),'excludeNegVoxels',1);
        R.crossval = ones(length(R.corr),1);
        RR = addstruct(RR,R);
        fprintf('%d correl\n',sn);
        R.corr=mva_spatialCorr(D,BIN);
        R.crossval = zeros(length(R.corr),1);
        RR = addstruct(RR,R);
        save(outName,'-struct','RR');
    case 'EVAL:averageClust' % make new 'spatialBoundfunc4.mat' struct. [4] - average eval corr across studies
        clusterNum=varargin{1}; % [5:2:23]
        condType=varargin{2}; % evaluating on 'unique' or 'all' taskConds ?
        
        mapType=[1,2]; % test on one dataset
        eval=[2,1]; % evaluate on the other
        R=[];
        for i=1:2,
            T=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_SC%d_%dcluster',mapType(i),clusterNum),sprintf('spatialBoundfunc%d_%s.mat',eval(i),condType)));
            T.studyNum=repmat([i],length(T.SN),1);
            R=addstruct(R,T);
        end
        R=rmfield(R,{'distmin','distmax','N'});
        % get average of both structures here
        A=tapply(R,{'bin','SN','bwParcel','crossval','dist'},{'corr'});
        
        outDir=fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_SC2_%dcluster',clusterNum),sprintf('spatialBoundfunc4_%s.mat',condType));
        save(outDir,'-struct','A');
    case 'EVAL:averageOther'
        mapType=varargin{1}; % options are 'lob10','bucknerRest','atlasFinal<num>' etc'
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
        
        outDir=fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc4_%s.mat',condType));
        save(outDir,'-struct','A');
        
    case 'EVAL:PLOT:CURVES'
        mapType=varargin{1}; % options are 'lob10','lob26','bucknerRest','SC<studyNum>_<num>cluster', or 'SC<studyNum>_POV<num>'
        data=varargin{2}; % evaluating data from study [1] or [2], both [3] or average of [1] and [2] after eval [4]
        type=varargin{3}; % 'group' or 'leaveOneOut'
        crossval=varargin{4}; % [0] - no crossval; [1] - crossval
        condType=varargin{5}; % evaluating on 'all' or 'unique' taskConds ??
        
        switch type,
            case 'group'
                T=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_%s.mat',data,condType)));
            case 'leaveOneOut'
                T=[];
                for s=1:length(returnSubjs),
                    P=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_eval_%s.mat',data,subj_name{returnSubjs(s)})));
                    T=addstruct(T,P);
                end
            otherwise
                fprintf('no such case')
        end
        
        xyplot(T.dist,T.corr,T.bin,'split',T.bwParcel,'leg',{'within Parcel','between Parcels'},...
            'subset',T.crossval==crossval & T.dist<=35,'style_symbols4*2','markersize',8,...
            'markertype','^','linewidth',2,'linecolor',{'r','k'},'markercolor',{'r','k'},...
            'markerfill',{'r','k'});
        set(gca,'YLim',[0 0.6],'XLim',[0 35]);
        xlabel('Spatial Distance (mm)');
        ylabel('Activity correlation');
        title(sprintf('%s-func%d-%dcrossval',mapType,data,crossval)); %0-no crossval; %1-crossval
    case 'EVAL:PLOT:DIFF'
        mapType=varargin{1}; % {'lob10','bucknerRest','atlasFinal9'}
        data=varargin{2}; % evaluating data from study [1] or [2] or [3]?
        crossval=varargin{3}; % [0]-no crossval; [1]-crossval
        condType=varargin{4}; % evaluating on 'unique' or 'all' taskConds ??
        
        P=[];
        for m=1:length(mapType),
            T=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType{m}),sprintf('spatialBoundfunc%d_%s.mat',data,condType)));
            
            A=getrow(T,T.crossval==crossval);
            A.type=repmat({mapType{m}},length(A.bin),1);
            A.m=repmat(m,length(A.bin),1);
            P=addstruct(P,A);
            clear A
        end
        %         pivottable(T.SN,[T.bin T.bwParcel],T.corr,'length','subset',T.crossval==crossval)
        
        % plot boxplot of different clusters
        W=getrow(P,P.bwParcel==0); % within
        B=getrow(P,P.bwParcel==1); % between
        W.diff=W.corr-B.corr;
        
        myboxplot(W.m,W.diff,'subset',W.dist<=35,'style_twoblock','plotall',0) % within-between diff
        ylabel('Within/Between Difference');
        
        figure()
        lineplot(W.m,W.diff,'subset',W.dist<=35); 
    case 'EVAL:PLOT:(UN)CROSSVAL'
        mapType=varargin{1}; % {'lob10','bucknerRest','atlasFinal9'}
        data=varargin{2}; % evaluating data from study [1] or [2] or [4]? (not [3] - because there is no crossval)
        condType=varargin{3}; % evaluating on 'unique' or 'all' ?
        P=[];
        for m=1:length(mapType),
            T=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType{m}),sprintf('spatialBoundfunc%d_%s.mat',data,condType)));
            T.type=repmat({mapType{m}},length(T.bin),1);
            T.m=repmat(m,length(T.bin),1); 
            P=addstruct(P,T);
            clear T
        end
        %         pivottable(T.SN,[T.bin T.bwParcel],T.corr,'length','subset',T.crossval==crossval)
        
        % plot boxplot of different clusters
        W=getrow(P,P.bwParcel==0); % within
        B=getrow(P,P.bwParcel==1); % between
        W.diff=W.corr-B.corr;
        
        lineplot([W.type],W.diff,'split',W.crossval,'leg',{'uncrossval','crossval'},'subset',W.dist<=30,'style_shade')
        ylabel('Within/Between Diff')
    case 'EVAL:PLOT:INDIV' % NEED TO UPDATE FOR SNN !!
        study=varargin{1};
        var=varargin{2};
        data=varargin{3};
        type=varargin{4}; % 'leaveOneOut' or 'Indiv'
        
        switch type,
            case 'indiv'
                T=load(fullfile(studyDir{study},'encoding','glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100),sprintf('spatialBound%s.mat',data)));
            case 'leaveOneOut'
                T=[];
                for s=1:length(returnSubjs),
                    P=load(fullfile(studyDir{study},'encoding','glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100),sprintf('spatialBound%s_eval_%s.mat',data,subj_name{returnSubjs(s)})));
                    T=addstruct(T,P);
                end
            otherwise
                fprintf('no such case')
        end
        
        figure('Name',type)
        for s=1:length(unique(T.SN)),
            %         subplot(2,1,1);
            %         xyplot(T.dist,T.corr,T.bin,'split',T.bwParcel,'leg',{'within Parcel','between Parcels'},'subset',T.crossval==0 ,'style_thickline');
            %         set(gca,'YLim',[0 0.5],'XLim',[0 76]);
            %         xlabel('Spatial Distance (mm)');
            %         ylabel('Activity correlation');
            %         title(sprintf('SC%d(%dPOV)-%s-without crossval',study,var*100,data));
            %         subplot(2,1,2);
            boxplot(T.SN,T.corr,'split',T.bwParcel,'leg',{'within Parcel','between Parcels'},'subset',T.crossval==1 & T.SN==returnSubjs(s),'style_thickline');
        end
        set(gca,'YLim',[0 0.5],'XLim',[0 76]);
        xlabel('Spatial Distance (mm)');
        ylabel('Activity correlation');
        title(sprintf('SC%d(%dPOV)-%s-with crossval',study,var*100,data));
        set(gcf,'PaperPosition',[2 4 10 12]);
        wysiwyg;
    case 'FIGURES:CORR' % Makes the summary figure of within / between correlations
        toPlot = {'SC2_5cluster','SC2_7cluster','SC2_9cluster','SC2_11cluster','SC2_13cluster','SC2_15cluster','SC2_17cluster'};
        crossval=varargin{1}; % 0-uncrossval; 1-crossval
        evalNum=varargin{2}; % SC1-[1],SC2-[2],concat SC1+SC2 [3],average of SC1+SC2 eval [4]
        condType=varargin{3}; % evaluating on 'unique' or 'all' taskConds ?? 
        numPlots = numel(toPlot);
        for i=1:numPlots
            subplot(1,numPlots,i);
            sc1_sc2_functionalAtlas('EVAL:PLOT:CURVES',toPlot{i},evalNum,'group',crossval,condType);
        end;
        %         set(gcf,'PaperPosition',[2 4 12 3]);
        %         wysiwyg;
    case 'FIGURES:DIFF' % Makes the summary figure of within / between diff
        toPlot = {'SC2_5cluster','SC2_7cluster','SC2_9cluster','SC2_11cluster','SC2_13cluster','SC2_15cluster','SC2_17cluster','SC2_19cluster','SC2_21cluster','SC2_23cluster'}; 
        crossval=varargin{1}; % 0-uncrossval; 1-crossval
        evalNum=varargin{2}; % SC1-[1],SC2-[2],concat SC1+SC2 [3],average of SC1+SC2 eval [4]
        condType=varargin{3}; % evaluating on 'unique' or 'all' taskConds ?? 
        sc1_sc2_functionalAtlas('EVAL:PLOT:DIFF',toPlot,evalNum,crossval,condType);
        set(gcf,'PaperPosition',[2 4 12 3]);
        %         wysiwyg;    
        
    case 'ENCODE:project_taskSpace'
        K=varargin{1}; % how many clusters ?
        
        % load raw Y (cereb)
        [Y,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',returnSubjs,[1,2],'no');
        
        % load clusters for cereb
        load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_atlasFinal%d',K),'parcels.mat'));
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % project back to task-space
        U=Y*G; % get U * S
        L=U'*U;
        S=diag(diag(sqrt(L)));
        X=U/S;
        
        % visualise clusters in task-space
        P=caret_load(fullfile(caretDir,'suit_flat','glm4','atlasFinal.paint'));
        figure();imagesc_rectangle(X,'YDir','reverse');
        caxis([0 1]);
        t=set(gca,'Ytick', 1:length(D.condNames),'YTickLabel',D.condNames,'FontSize',12,'FontWeight','bold',...
            'Xtick',1:size(X,2),'XTickLabel',P.paintnames);
        t.Color='white';
        
        T.taskWeights=X;
        T.taskNames=D.condNames;
        T.featNames=P.paintnames;
        
        save(fullfile(studyDir{2},encodeDir,'glm4','atlasFinal_condLoads.mat'),'T')
        
        % Plot left/right hands
        
        %         % scatterplots
        %         CAT.markertype='o';
        %         CAT.markersize=8;
        %
        %         % plot left/right hands
        %         scatterplot(X(:,14),X(:,10),'label',allCondNames,'CAT',CAT,'regression','linear','intercept',0,'printcorr','draworig')
        %         xlabel('right hand tasks');ylabel('left hand tasks');
    case 'ENCODE:features'
        K=varargin{1}; % how many clusters in atlasFinal ?
        
        D=dload(fullfile(baseDir,'featureTable_jd.txt'));        % Read feature table
        S=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));      % List of task conditions
        load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%datlasFinal',K),'condLoads.mat'));
        W=pivottablerow(S.condNumUni,T.taskWeights,'mean(x,1)');
        
        figure(1);
        subplot(2,2,1);
        imagesc_rectangle(W);
        
        set(gca,'YTick',[1:47],'YTickLabel',D.conditionName,'XTickLabel',T.featNames,'FontSize',6);
        
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
        numClusters = length(T.featNames)-1;
        
        subplot(2,2,2);
        C=corr(F,W);
        imagesc_rectangle(C);
        set(gca,'YTick',[1:numFeat],'YTickLabel',FeatureNames,'XTickLabel',T.featNames,'FontSize',6);
        
        subplot(2,2,[3:4]);
        imagesc_rectangle(F);
        set(gca,'YTick',[1:numCond],'YTickLabel',D.conditionName,'XTick',[1:numFeat],'XTickLabel',FeatureNames,'FontSize',6);
        set(gca,'XTickLabelRotation',60);
        
        % Run a multiple regression analysis on clusters onto features
        lambda = [0.001 0.001];
        X=bsxfun(@minus,F,mean(F,1));
        Y=bsxfun(@minus,W,mean(W,1));
        X=bsxfun(@rdivide,X,sqrt(mean(X.^2)));
        Y=bsxfun(@rdivide,Y,sqrt(mean(Y.^2)));
        XX=X'*X;
        XY=X'*Y;
        A = -eye(numFeat);
        b = zeros(numFeat,1);
        for p=1:numClusters,
            % U(:,p) =
            % cplexqp(XX+lambda(2)*eye(numFeat),ones(numFeat,1)*lambda(1)-XY(:,p),A,b);
            % % I don't have the IBM cplexqp routine (will download it
            % later)
            U(:,p) = (XX + eye(numFeat)*lambda(2))\(XY(:,p)); % just normal ridge regression
        end;
        
        % Present the list of the highest three correlation for each
        % cluster
        for i=1:numClusters,
            [a,b]=sort(C(:,i),'descend');
            fprintf('%s: %s(%2.2f) %s(%2.2f) %s(%2.2f)\n',T.featNames{i},...
                FeatureNames{b(1)},a(1),...
                FeatureNames{b(2)},a(2),...
                FeatureNames{b(3)},a(3));
        end;
        
        % Make word lists for word map
        figure(2);
        DD = U;WORD=FeatureNames;
        DD(DD<0)=0;
        DD=DD./max(DD(:));
        for j=1:size(DD,2)
            subplot(3,4,j);
            set(gca,'XLim',[0 2.5],'YLim',[-0.2 1.2]);
            title(T.featNames{j})
            for i=1:size(DD,1)
                if (DD(i,j)>0)
                    siz=ceil(DD(i,j)*20);
                    text(unifrnd(0,1,1),unifrnd(0,1,1),WORD{i},'FontSize',siz);
                end;
            end;
        end;
        set(gcf,'PaperPosition',[1 1 60 30]);wysiwyg;
        

end
% Local functions
function dircheck(dir)
if ~exist(dir,'dir');
    warning('%s doesn''t exist. Creating one now. You''re welcome! \n',dir);
    mkdir(dir);
end

% InterSubj Corr
function C=interSess_corr(Y)
numSess=size(Y,3);
for i=1:numSess
    for j=1:numSess
        C(i,j)=nansum(nansum(Y(:,:,i).*Y(:,:,j)))/...
            sqrt(nansum(nansum(Y(:,:,i).*Y(:,:,i)))*...
            nansum(nansum(Y(:,:,j).*Y(:,:,j))));
    end;
end




