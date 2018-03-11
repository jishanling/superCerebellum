function varargout=sc1_sc2_functionalAtlas(what,varargin)

% Directories
baseDir          = '/Users/maedbhking/Documents/Cerebellum_Cognition';
%baseDir            = '/Volumes/MotorControl/data/super_cerebellum_new';
% baseDir          = '/Users/jdiedrichsen/Data/super_cerebellum_new';

studyDir{1}     =fullfile(baseDir,'sc1');
studyDir{2}     =fullfile(baseDir,'sc2');
suitDir         ='/suit';
caretDir        ='/surfaceCaret';
regDir          ='/RegionOfInterest/';
encodeDir       ='/encoding';

subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11',...
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
                CORR(c,:,:,subj)=interSess_corr(Yf(c,:,:,subj));
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
    case 'FIGURES:reliabilityA'
        % load relability
        load(fullfile(studyDir{2},regDir,'glm4','patternReliability.mat'));
        
        figure();lineplot(S.condNum,[S.within1,S.within2,S.across],'leg',{'within1','within2','across'})
        A=tapply(S,{'SN'},{'across'},{'within1'},{'within2'});
        
        % within & between-subj reliability
        myboxplot([],[A.within1 A.within2 A.across],'style_tukey');
        drawline(0,'dir','horz');
        set(gca,'XTickLabel',{'SC1','SC2','across'});
        set(gcf,'PaperPosition',[2 2 4 3]);
        wysiwyg;
        ttest(sqrt(A.within1.*A.within2),A.across,2,'paired');
        
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
        
        save(fullfile(studyDir{2},regDir,'glm4','distanceReliability.mat'),'T','dist')
    case 'FIGURES:RDM'
        % load in fullRDM
        [fullRDM,T,~]=sc1_sc2_functionalAtlas('REPRESENTATION:get_distances','cerebellum','none');
        
        % Plot RDM
        reOrder=[1,2,6,7,8,9,10,11,12,13,14,17,18,22,23,3,4,5,15,16,19,20,21,24,25,26,...
            27,28,29,58,59,60,43,44,49,48,36,34,35,55,56,57,61,30,31,32,33,37,38,39,40,41,42,45,46,47,50,51,52,53,54]'; % reorder
        
        % threshold RDM
        fullRDM_thresh=reshape(fullRDM,[size(fullRDM,1)*size(fullRDM,2)],[]);
        for dd=1:size(fullRDM_thresh,1),
            [t(dd),p(dd)]=ttest(fullRDM_thresh(dd,:),[],1,'onesample');
        end
        ut=t; % unthresholded distances
        t(p>.0001)=0; % thresholded distances
        
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
        t=set(gca,'Ytick',[1:numDist]','YTickLabel',T.condNames(reOrder)','FontSize',7);
        t.Color='white';
        colorbar
        
        % visualise RDM (distances - unthresh)
        figure()
        imagesc_rectangle(avrgFullRDM(reOrder,reOrder),'YDir','reverse');
        caxis([0 1]);
        %         t=set(gca,'Ytick',[1:numDist]','YTickLabel',T.condNames(reOrder)','FontSize',12,'FontWeight','bold');
        t=set(gca,'Ytick',[1:numDist]','YTickLabel',T.condNames(reOrder)','FontSize',7);
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
        
        condIndx=[1:length(T.condNum)]';
        indxShort=indx;
        condNames=T.condNames;
        CAT.markercolor= {colour{indxShort}};
        CAT.markerfill = {colour{indxShort}};
        CAT.markersize = 10;
        
        X1=X(condIndx,condIndx);
        figure()
        scatterplot3(X1(:,1),X1(:,2),X1(:,3),'split',condIndx,'CAT',CAT,'label',condNames,'markertype','o');
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'Box','on');
        hold on;
        plot3(0,0,0,'+');
        % Draw connecting lines
        for i=1:15,
            ind=clustTree(i,1:2);
            X(end+1,:)=(X(ind(1),:)+X(ind(2),:))/2;
            line(X(ind,1),X(ind,2),X(ind,3));
        end;
        hold off;
        set(gcf,'PaperPosition',[2 2 8 8]);
        axis equal;
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'Box','on');
        wysiwyg;
        view([81 9]);
        clear X1  indxShort
    case 'FIGURES:reliabilityD'
        % load relability
        load(fullfile(studyDir{2},regDir,'glm4','distanceReliability.mat'));
        
        % within & between-subj reliability
        myboxplot([],[T.within1 T.within2 T.across],'style_tukey');
        set(gca,'XTickLabel',{'SC1','SC2','across'});
        set(gcf,'PaperPosition',[2 2 4 3]);
        wysiwyg;
        ttest(sqrt(T.within1.*T.within2),T.across,2,'paired');
        
    case 'CONVERT:mni2suit'
        inputImages={'N2C_subcortex_atlas_subcortexGSR.nii',...
            'Buckner2011_7Networks_MNI152_FreeSurferConformed1mm_TightMask.nii'};
        
        inputDir=which(inputImages{2});
        cd(fileparts(inputDir))
        suit_mni2suit(inputImages{2})
        
    case 'MAP:vol2surf'
        % this function takes any labelled volume (already in SUIT space)
        % and plots to the surface
        inputMap=varargin{1}; % some options are 'Buckner_7Networks','SC1_9cluster','lob10', 'Cole_10Networks', 'SC2_90cluster' etc
        metric=varargin{2}; % 'yes' or 'no'
        
        inputDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',inputMap));
        cd(inputDir);
        
        Vo=spm_vol(fullfile(inputDir,'map.nii'));
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
                suit_plotflatmap(M.data,'type','label','cmap',cmapF) % colorcube(max(M.data))
        end
    case 'MAP:optimal'  % figure out optimal map for multiple clusters
        % example:sc1_sc2_functionalAtlas('MAP:optimal',<subjNums>,1,6,'group')
        sn=varargin{1};     % subject numbers
        study=varargin{2};  % 1 or 2 or [1,2]
        K=varargin{3};      % K=numClusters (i.e. 5);
        type=varargin{4};   % 'group' or 'indiv'
        numCount=5;         % How often the "same" solution needs to be found
        
        tol_rand = 0.90;    % Tolerance on rand coefficient to call it the same solution
        maxIter=100; % if it's not finding a similar solution - force stop at 100 iters
        
        % Set the String correctly
        studyStr = sprintf('SC%d',study);
        if length(study)>1
            studyStr='AtlasFinal';
        end
        
        % Set output File name
        switch type,
            case 'group'
                sn=returnSubjs;
                outDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s_%dcluster',studyStr,K));
                dircheck(outDir);
                outName=fullfile(outDir,'SNN.mat');
            case 'indiv'
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
                end;
            else                     % Different (enough) solution
                if (Info.error<bestErr) % Is this a better solution
                    bestErr = Info.error;
                    bestSol = winner;
                    bestG   = G;
                    bestF   = F;
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
        save(outName,'bestG','bestF','errors','randInd','iter','count','volIndx','V');
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
            studyStr='AtlasFinal';
        end
        
        % Set output File name
        switch type,
            case 'group'
                sn=returnSubjs;
                outDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s_%dPOV',studyStr,K));
                dircheck(outDir);
                outName=fullfile(outDir,'ICAs.mat');
            case 'indiv'
                outDir=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn});
                outName=fullfile(outDir,sprintf('ICAs_%s_%dPOV.mat',studyStr,K));
        end
        
        % get data
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',sn,study,'build');
        
        threshold=threshold/100;
        [A_PW,S_PW,W_PW,winner]=pca_ica(X_C,'threshold',threshold);
        
        save(outName,'A_PW','S_PW','W_PW','volIndx','V');
        
        % Intialize iterations[G
        %         bestSol = ones(size(X_C,1),1);
        %         iter=1; % How many iterations
        %         count=0;
        %         K=K/100;
        %         while iter<maxIter,
        %             [F,G,~,winner]=pca_ica(X_C,'K',K); % .9 instead of 90
        %             randInd(iter)=RandIndex(bestSol,winner);
        %
        %             % Check if we have a similar solution
        %             if randInd(iter)>tol_rand % Similar solution
        %                 count=count+1;       % count up one
        %                 bestSol = winner;
        %                 bestG   = G;
        %                 bestF   = F;
        %                 fprintf('RandIndex is %2.2f \n',randInd);
        %             else                     % Different (enough) solution
        %                 bestSol = winner;
        %                 bestG   = G;
        %                 bestF   = F;
        %                 count = 0;         % first time we found this solution: reset counter
        %                 fprintf('RandIndex is %2.2f \n',randInd);
        %             end;
        %             if count>=numCount || iter>=maxIter,
        %                 fprintf('Existing loop....\n');
        %                 break;
        %             end;
        %             iter=iter+1;
        %         end;
        %         save(outName,'bestG','bestF','randInd','iter','count','volIndx','V');
    case 'MAP:visualise'
        sn=varargin{1}; % [2] or 'group'
        study=varargin{2}; % 1 or 2 or [1,2]
        anaType=varargin{3}; % 'SNN' or 'ICAs'
        K=varargin{4}; % for SNN - number of clusters (i.e. 5), for ica - thresh (i.e. 90)
        type=varargin{5}; % 'group', or 'indiv'
        
        % figure out if ICA or SNN
        if strcmp(anaType,'SNN'),
            anaName='cluster';
        else
            anaName='POV';
        end
        
        % figure out if individual or group
        switch type,
            case 'group'
                outName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC%d_%d%s',study,K,anaName),'map.nii');
                load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC%d_%d%s',study,K,anaName),sprintf('%s.mat',anaType)));
            case 'indiv'
                outName=fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('map_SC%d_%d%s.nii',study,K,anaName));
                load(fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('%s_SC%d_%d%s.mat',anaType,study,K,anaName)));
        end
        
        % transpose matrix from ICA
        if size(S_PW,1)<1000,
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
        exampleVol=fullfile(studyDir{study},suitDir,'glm4','s02','wdbeta_0001.nii'); % must be better way ??
        X=spm_vol(exampleVol);
        X.fname=outName;
        X.private.dat.fname=V.fname;
        spm_write_vol(X,Vv{1}.dat);
    case 'MAP:compare' % this function is will compute randindex between maps (pairwise if more than 2 maps are provided)
        maps=varargin{1}; % give cell with any number of maps {'Cole_10Networks','Buckner_7Networks'} etc
        % example
        % sc1_sc2_functionalAtlas('MAP:compare',{'Cole_10Networks','Buckner_7Networks,'SC1_10cluster'})
        numMaps=length(maps);
        
        % get data
        [~,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',2,1,'build'); % just to get volIndx + V
        
        % make sure all parcels are sampled into the same space
        for m=1:numMaps,
            mapName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',maps{m}),'map_old.nii');
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
        varargout={RI};
    case 'MAP:final'
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
                data=suit_map2surf(Vv,'space','SUIT','stats','mode');
                suit_plotflatmap(data,'type','label','cmap',colorcube(K))
        end
    case 'MAP:plotMapErr'
        sn=varargin{1};     % 'group' or subjNum
        study=varargin{2};  % 1 or 2
        anaType=varargin{3}; % 'SNN' or 'cluster'
        K=varargin{3};      % Number of clusters (i.e. 5) or ICA thresh (i.e. 90 etc)
        type=varargin{4};   % 'group' or 'indiv'
        
        % figure out if ICA or SNN
        if strcmp(anaType,'SNN'),
            anaName='cluster';
        else
            anaName='POV';
        end
        
        % figure out if individual or group
        switch type,
            case 'group'
                load(fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_SC%d_%d%s',study,K,anaName),sprintf('%s.mat',anaType)));
            case 'indiv'
                load(fullfile(studyDir{2},encodeDir,'glm4',subj_name{sn},sprintf('%s_SC%d_%d%s.mat',anaType,study,K,anaName)));
        end;
        
        % plot bestF
        subplot(2,2,[1,2])
        title(sprintf('Map with %d clusters',K))
        imagesc_rectangle(bestF,'YDir','reverse')
        caxis([0 1]);
        t=set(gca,'Xtick',[1:size(bestF,2)]','Ytick',[1:size(bestF,1)]');
        xlabel('clusters')
        ylabel('taskConds')
        
        % plot errors
        if strcmp(anaType{1},'SNN'),
            subplot(2,2,3)
            lineplot([1:length(errors)]',errors')
            ylabel('errors')
        end
        
        % plot randInd
        subplot(2,2,4)
        lineplot([1:length(randInd)]',randInd')
        ylabel('Rand Index')
    case 'FIGURES:similarityMatrix' % this function will plot a matrix of map similarities (randIndex is the metric of similarity)
        anaType=varargin{1}; % 'cluster' or 'POV'
        toPlotNames={'1.10Cluster','2.Cole','3.7Cluster','4.Buckner7','5.17Cluster','6.Buckner17'};
        
        K=[50:5:95];
        idx=1;
        for s=1:2,
            for k=1:length(K),
                toPlot{idx}=sprintf('SC%d_%d%s',s,K(k),anaType);
                toPlotNames{idx}=sprintf('%d.SC%d-%d',idx,s,K(k));
                idx=idx+1;
            end
            idx=k+1;
        end
        
        
        [RI]=sc1_sc2_functionalAtlas('MAP:compare',toPlot);
        
        figure()
        imagesc_rectangle(RI,'YDir','reverse');
        caxis([0 1]);
        t=set(gca,'Ytick',[1:length(toPlot)]','YTickLabel',toPlotNames','FontSize',7,'Xtick',[1:length(toPlot)]','FontSize',10);
        t.Color='white';
        colorbar
        
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
        mapType=varargin{1}; % options are 'lob10','lob26','Buckner_7Networks','Buckner_17Networks','Cole_12Networks', or 'atlasFinal<num>'
        evalType=varargin{2}; % [1] [2] or [1,2]
        
        if size(evalType,2)==2,
            outName='spatialBoundfunc3_all.mat';
            outDir=fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),outName);
        else
            outName=sprintf('spatialBoundfunc%d.mat',evalType);
            outDir=fullfile(studyDir{evalType},'encoding','glm4',sprintf('groupEval_%s',mapType),outName);
        end
        
        % load in map
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
    case 'EVAL:crossval'
        sn=varargin{1}; % 'group' or <subjNum>
        mapType=varargin{2}; % options are 'lob10','lob26','Buckner_17Networks','Buckner_7Networks', 'Cole_10Networks','SC<studyNum>_<num>cluster'
        data=varargin{3}; % evaluating data from study [1] or [2] ?
        condType=varargin{4}; % 'unique' or 'all'. Are we evaluating on all taskConds or just those unique to either sc1 or sc2 ?
        type=varargin{5}; % 'group' or 'indiv'
        
        switch type,
            case 'group'
                % load in map
                mapName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'map.nii');
                outName=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_%s.mat',data,condType));
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
        [BIN,R]=mva_spatialCorrBin(XYZ(:,voxIn),'Parcel',Parcel(1,voxIn));
        clear XYZ i k l x y z i1 j1 k1 VA Parcel; % Free memory
        % Now calculate the uncrossvalidated and crossvalidated
        % Estimation of the correlation for each subject
        for s=unique(T.SN)';
            for c=1:length(idx),
                i1(c) = find(T.SN==s & T.sess==1 & T.cond==idx(c));
                i2(c) = find(T.SN==s & T.sess==2 & T.cond==idx(c));
            end
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
                    T=load(fullfile(studyDir{2},'encoding','glm4',subj_name{sn},sprintf('SC%d_%s_spatialBoundfunc%d_%s.mat',studyType(i),mapType,clusterNum,evalType(i),condType)));
                    T.studyNum=repmat([i],length(T.SN),1);
                    R=addstruct(R,T);
                end
                outDir=fullfile(studyDir{2},'encoding','glm4',subj_name{sn},sprintf('%dcluster_spatialBoundfunc4_%s.mat',clusterNum,condType));
        end
        R=rmfield(R,{'distmin','distmax','N'});
        % get average of both structures here
        A=tapply(R,{'bin','SN','bwParcel','crossval','dist'},{'corr'});
        
        save(outDir,'-struct','A');
    case 'EVAL:averageOther'
        mapType=varargin{1}; % options are 'lob10','Buckner_7Networks','Buckner_17Networks','Cole_10Networks','atlasFinal<num>' etc'
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
        type=varargin{3}; % 'group' or 'leaveOneOut' or 'indiv'
        crossval=varargin{4}; % [0] - no crossval; [1] - crossval
        condType=varargin{5}; % evaluating on 'all' or 'unique' taskConds ??
        
        vararginoptions({varargin{6:end}},{'sn'}); % option if doing individual map analysis
        
        switch type,
            case 'group'
                T=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_%s.mat',data,condType)));
            case 'leaveOneOut'
                T=[];
                for s=1:length(returnSubjs),
                    P=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_eval_%s.mat',data,subj_name{returnSubjs(s)})));
                    T=addstruct(T,P);
                end
            case 'indiv'
                T=load(fullfile(studyDir{2},'encoding','glm4',subj_name{sn},sprintf('%s_spatialBoundfunc%d_%s.mat',mapType,data,condType)));
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
        
        vararginoptions({varargin{5:end}},{'K'}); % option if plotting clusters
        
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
        
        figure()
        if exist('K'),
            lineplot(W.m,W.diff,'subset',W.dist<=35,'style_shade');
            t=set(gca,'XTickLabel',K');
        else
            lineplot(W.m,W.diff,'subset',W.dist<=35,'style_shade');
        end
        ylabel('Within/Between Diff')
        xlabel('Maps')
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
        
        lineplot([W.type],W.diff,'split',W.crossval,'leg',{'uncrossval','crossval'},'subset',W.dist<=35,'style_shade')
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
        
    case 'FIGURES:CORR:GROUP' % Makes the summary figure of within / between correlations for GROUP
        toPlot = {'lob10','Buckner_7Networks','Buckner_17Networks','Cole_10Networks'};
        crossval=varargin{1}; % 0-uncrossval; 1-crossval
        evalNum=varargin{2}; % SC1-[1],SC2-[2],concat SC1+SC2 [3],average of SC1+SC2 eval [4]
        condType=varargin{3}; % evaluating on 'unique' or 'all' taskConds ??
        
        numPlots = numel(toPlot);
        figure()
        for i=1:numPlots
            subplot(1,numPlots,i);
            sc1_sc2_functionalAtlas('EVAL:PLOT:CURVES',toPlot{i},evalNum,'group',crossval,condType);
        end;
        
        %                 set(gcf,'PaperPosition',[2 4 12 3]);
        %                 fig.PaperPositionMode = 'auto';
        %                 wysiwyg;
    case 'FIGURES:CORR:INDIV'
        sn=varargin{1}; % <subjNum>
        crossval=varargin{2}; % 0-uncrossval; 1-crossval
        evalNum=varargin{3}; % SC1-[1],SC2-[2],concat SC1+SC2 [3],average of SC1+SC2 eval [4]
        condType=varargin{4}; % evaluating on 'unique' or 'all' taskConds ??
        
        numPlots = numel(sn);
        figure()
        for i=1:numPlots
            subplot(1,numPlots,i);
            sc1_sc2_functionalAtlas('EVAL:PLOT:CURVES','9cluster',evalNum,'indiv',crossval,condType,'sn',sn(i));
        end;
        %         set(gcf,'PaperPosition',[2 4 12 3]);
        %         wysiwyg;
    case 'FIGURES:DIFF:GROUP' % Makes the summary figure of within / between diff
        % example: toPlot= {'lob10','Buckner_7Networks','Buckner_17Networks','Cole_10Networks','SC2_10cluster'};
        crossval=varargin{1}; % 0-uncrossval; 1-crossval
        evalNum=varargin{2}; % SC1-[1],SC2-[2],concat SC1+SC2 [3],average of SC1+SC2 eval [4]
        condType=varargin{3}; % evaluating on 'unique' or 'all' taskConds ??
        studyNum=varargin{4}; % if evaluating on 1,4: studyNum=2; if evaluating on 2,studyNum=1
        
        %K=[50:5:95,5:24]; % whatever is hardcoded here will determine ICA or SNN
        % what are we plotting ? SNN or ICA maps ?
        %         for k=1:length(K),
        %             if max(K)<50,
        %                 toPlot{k}=sprintf('SC%d_%dcluster',studyNum,K(k));
        %             else
        %                 toPlot{k}=sprintf('SC%d_%dPOV',studyNum,K(k));
        %             end
        %         end
        toPlot={'lob10','Buckner_7Networks','Buckner_17Networks','Cole_10Networks','SC2_9cluster','SC2_90POV'};
        
        
        if exist('K'),
            sc1_sc2_functionalAtlas('EVAL:PLOT:DIFF',toPlot,evalNum,crossval,condType,'K',K);
        else
            sc1_sc2_functionalAtlas('EVAL:PLOT:DIFF',toPlot,evalNum,crossval,condType);
        end
        set(gcf,'PaperPosition',[2 4 12 3]);
        %         wysiwyg;
    case 'FIGURES:DIFF:INDIV'
        sn=varargin{1}; % <subjNum>
        crossval=varargin{2}; % [0]-no crossval; [1]-crossval
        condType=varargin{3}; % evaluating on 'unique' or 'all' taskConds ??
        
        for s=1:length(sn),
            P=[];
            
            T=load(fullfile(studyDir{2},'encoding','glm4',subj_name{sn(s)},sprintf('9cluster_spatialBoundfunc4_%s.mat',condType)));
            A=getrow(T,T.crossval==crossval);
            P=addstruct(P,A);
            clear A
            %         pivottable(T.SN,[T.bin T.bwParcel],T.corr,'length','subset',T.crossval==crossval)
            
            % plot boxplot of different clusters
            W=getrow(P,P.bwParcel==0); % within
            B=getrow(P,P.bwParcel==1); % between
            W.diff=W.corr-B.corr;
            
            subplot(1,length(sn),s);
            lineplot(W.SN,W.diff,'subset',W.dist<=35,'style_shade','errorfcn','stderr(fisherz(x))');
            ylabel('Within/Between Diff')
            title(sprintf('%s',subj_name{sn(s)}));
            clear W
        end
    case 'FIGURES:GROUP:INDIV'
        sn=varargin{1}; % <subjNum>
        crossval=varargin{2}; % [0]-no crossval; [1]-crossval
        condType=varargin{3}; % evaluating on 'unique' or 'all' taskConds ??
        
        % load in group map
        T=load(fullfile(studyDir{2},'encoding','glm4','groupEval_SC2_9cluster',sprintf('spatialBoundfunc4_%s.mat',condType)));
        T.mapType=repmat({'group'},length(T.SN),1);
        T.type=repmat([1],length(T.SN),1); % 1 is group
        G=getrow(T,T.crossval==crossval);
        
        A=[];
        for s=1:length(sn),
            P=load(fullfile(studyDir{2},'encoding','glm4',subj_name{sn(s)},sprintf('9cluster_spatialBoundfunc4_%s.mat',condType)));
            P.mapType=repmat({'indiv'},length(P.SN),1);
            P.type=repmat([2],length(T.SN),1); % 2 is indiv
            I=getrow(P,P.crossval==crossval & P.SN==sn(s));
            A=addstruct(A,I);
            clear I P
        end
        G=addstruct(G,A);
        
        W=getrow(G,G.bwParcel==0); % within
        B=getrow(G,G.bwParcel==1); % between
        W.diff=W.corr-B.corr;
        
        % integrate over dist (large error bars !!)
        D=getrow(W,W.dist<=35);
        Dd=tapply(D,{'SN','mapType','type'},{'diff'});
        
        lineplot(Dd.SN,Dd.diff,'split',Dd.mapType,'leg','auto','style_thickline2x3')
        
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






