function varargout=sc1_sc2_functionalAtlas(what,varargin)

% Directories
baseDir          = '/Users/maedbhking/Documents/Cerebellum_Cognition';
% baseDir            = '/Volumes/MotorControl/data/super_cerebellum';
baseDir          = '/Users/jdiedrichsen/Data/super_cerebellum_new';

studyDir{1}     =fullfile(baseDir,'sc1');
studyDir{2}     =fullfile(baseDir,'sc2');
suitDir         ='/suit';
caretDir        ='/surfaceCaret';
regDir          ='/RegionOfInterest/';
encodeDir       ='/encoding';

subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11',...
    's12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24',...
    's25','s26','s27','s28','s29','s30','s31'};

returnSubjs=[2,3,4,6,7,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31];

% Remove feats for functional maps
removeFeats={[],... % SC1_Sess1
    [],... % SC1_Sess2
    [],... % SC2_Sess1
    []};   % SC2_Sess2

switch what
    
    case 'TASKSPACE:overlap'
    case 'TASKSPACE:contrasts'
        
    case 'EVAL:runGroup'
        study=varargin{1}; % 1 or 2 or [1,2]
        var=varargin{2};% examples: .95 or .90
        type=varargin{3}; % 'generalisability' or 'reliability'
        
        sc1_sc2_functionalAtlas('ICA:make_map','group',study,var*100,'group')
        sc1_sc2_functionalAtlas('ICA:visualise_map','group',study,var*100,'group')
        switch type
            case 'generalisability'
                % Generalisability
                sc1_sc2_functionalAtlas('ICA:funcBound:GROUP','SC1',var*100,'func2');
                sc1_sc2_functionalAtlas('ICA:funcBound:GROUP','SC2',var*100,'func1');
            case 'reliability'
                % Reliability
                sc1_sc2_functionalAtlas('ICA:funcBound:GROUP','SC1',var*100,'func1');
                sc1_sc2_functionalAtlas('ICA:funcBound:GROUP','SC2',var*100,'func2');
        end
    case 'EVAL:runIndiv'
        var={.80,.95}; % 80% for SC1 and 95% for SC2
        
        for s=1:length(returnSubjs),
            for study=1:2,
                sc1_sc2_functionalAtlas('ICA:make_map',returnSubjs(s),study,var{study},'indiv')
                sc1_sc2_functionalAtlas('ICA:visualise_map',returnSubjs(s),study,var{study},'indiv')
            end
            sc1_sc2_functionalAtlas('ICA:funcBound:INDIV',returnSubjs(s),1,var{1},'func2')
            sc1_sc2_functionalAtlas('ICA:funcBound:INDIV',returnSubjs(s),2,var{2},'func1')
        end
    case 'EVAL:runLeaveOneOut'
        var={.80,.95};
        
        for s=1:length(returnSubjs),
            for study=1:2,
                sc1_sc2_functionalAtlas('ICA:make_map',returnSubjs(s),study,var{study},'leaveOneOut')
                sc1_sc2_functionalAtlas('ICA:visualise_map',returnSubjs(s),study,var{study},'leaveOneOut')
            end
            sc1_sc2_functionalAtlas('ICA:funcBound:LEAVEONEOUT',returnSubjs(s),1,var{1},'func2')
            sc1_sc2_functionalAtlas('ICA:funcBound:LEAVEONEOUT',returnSubjs(s),2,var{2},'func1')
        end
    case 'EVAL:get_data'
        sn=varargin{1}; % Subje numbers to include
        study=varargin{2}; % 1 or 2 or [1,2]

        % load data
        UFullAvrgAll=[];
        for f=1:length(study),
            load(fullfile(studyDir{study(f)},encodeDir,'glm4','cereb_avrgDataStruct.mat')); % JD Bug here: f instead of study(f)
            
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
            
            % if func1+func2 - concatenate
            if length(study)>1,
                UFullAvrgAll=[UFullAvrgAll;UFullAvrg];
            else
                UFullAvrgAll=UFullAvrg;
            end
        end
        
        % if group - get mean
        UFull=nanmean(UFullAvrgAll,3);
        
        % center the data
        X_C=bsxfun(@minus,UFull,mean(UFull));
        varargout={X_C,volIndx,V};
    case 'EVAL:make_ICA_map'
        % example: 'sc1_sc2_functionalAtlas('EVAL:make_ICA_map',[2],1,.95)
        sn=varargin{1}; % subjNum or 'group'
        study=varargin{2}; % 1 or 2 or [1,2]
        var=varargin{3}; % var explained. (.95,.90,.80 etc)
        type=varargin{4}; % 'group','indiv','leaveOneOut'
        
        vararginoptions({varargin{5:end}},{'var'}); % var explained. (.95,.90,.80 etc)
        
        % figure out if individual or group or leave one out
        switch type,
            case 'group'
                sn=returnSubjs;
                outDir=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100));
                outName=fullfile(outDir,'ICAs.mat');
            case 'indiv'
                outDir=fullfile(studyDir{study},encodeDir,'glm4',subj_name{sn});
                outName=fullfile(outDir,sprintf('ICAs_SC%d_%dPOV.mat',study,var*100));
            case 'leaveOneOut'
                outDir=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100));
                outName=fullfile(outDir,sprintf('ICAs_leaveOut_%s.mat',subj_name{sn}));
                sn=returnSubjs(returnSubjs~=sn);
        end
        
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('ICA:get_data',sn,study);
        % calculate the pca for data
        [E,D] = pcamat(X_C,1,size(X_C,1),'off','off'); % JD: Where is pcamat?
        
        % get PCs explaining <var> of POV
        [EV,I]=sort(diag(D),'descend');
        for i=1:length(EV),
            POV(i)=EV(i)/trace(D);
        end
        EVindx=find(cumsum(POV)<=var);
        
        % reconstruct E and D
        E=E(:,I);
        D=D(I,I);
        E=E(:,EVindx);
        D=D(EVindx,EVindx);
        
        % prewhiten the data
        whiteningMatrix = inv(sqrt (D)) * E';
        dewhiteningMatrix = E * sqrt (D);
        X_PW =  whiteningMatrix * X_C;
        
        % run ICA faster on centered data
        [S_PW,A_PW,W_PW] = fastica(X_C,'pcaE',E,'pcaD',D,'whiteSig',X_PW,'whiteMat',whiteningMatrix,'dewhiteMat',dewhiteningMatrix);
        
        dircheck(fullfile(outDir));
        save(fullfile(outName),'S_PW','A_PW','W_PW','X_PW','whiteningMatrix','dewhiteningMatrix','volIndx','V');   
    case 'EVAL:make_SNN_map'
       % example: 'sc1_sc2_functionalAtlas('EVAL:make_SNN_map',[2],1,13,'leaveOneOut')
        sn=varargin{1}; % subjNum or 'group'
        study=varargin{2}; % 1 or 2 or [1,2]
        K=varargin{3}; % number of clusters
        type=varargin{4}; % 'group','indiv','leaveOneOut'
        
        vararginoptions({varargin{5:end}},{'var'}); % var explained. (.95,.90,.80 etc)
        
        % figure out if individual or group or leave one out
        switch type,
            case 'group'
                sn=returnSubjs;
                outDir=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K));
                outName=fullfile(outDir,'SNN.mat');
            case 'indiv'
                outDir=fullfile(studyDir{study},encodeDir,'glm4',subj_name{sn});
                outName=fullfile(outDir,sprintf('SNN_SC%d_%dcluster.mat',study,K));
            case 'leaveOneOut'
                outDir=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K));
                outName=fullfile(outDir,sprintf('SNN_leaveOut_%s.mat',subj_name{sn}));
                sn=returnSubjs(returnSubjs~=sn);
        end
        
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',sn,study);
        [F,G,Err1]=semiNonNegMatFac(X_C,K,'threshold',0.01);
        save(fullfile(outName),'F','G','volIndx','V');
    case 'EVAL:visualise_ICA_map'
        sn=varargin{1}; % [2] or 'group'
        study=varargin{2}; % 1 or 2 or [1,2]
        var=varargin{3}; % {.95,.90} etc
        type=varargin{4};
        metric=varargin{5}; % save out metric file: 'yes' or 'no'
        
        feat=1;
        
        % figure out if individual or group
        switch type,
            case 'group'
                outName=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100),'map.nii');
                load(fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100),'ICAs.mat'));
            case 'indiv'
                outName=fullfile(studyDir{study},encodeDir,'glm4',subj_name{sn},sprintf('map_SC%d_%dPOV.nii',study,var*100));
                load(fullfile(studyDir{study},encodeDir,'glm4',subj_name{sn},sprintf('ICAs_SC%d_%dPOV.mat',study,var*100)));
            case 'leaveOneOut'
                outName=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100),sprintf('map_leaveOut_%s.nii',subj_name{sn}));
                load(fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100),sprintf('ICAs_leaveOut_%s.mat',subj_name{sn})));
        end
        
        numFeat=size(S_PW,1);
        
        % get pos and neg comps
        for p=1:size(S_PW,2),
            signIdx=find(S_PW(:,p)>0);
            pwS=abs(S_PW(:,p));
            [~,i]=sort(pwS,'descend');
            % sort 'winner' into pos or neg group
            if find(signIdx==i(1)),
                groupFeat(p,:)=i(1);
            else
                groupFeat(p,:)=i(1)+numFeat;
            end
            % reassign some ICs
            if find(removeFeats{feat}==groupFeat(p,:)),
                if groupFeat(p,:)>numFeat,
                    groupFeat(p,:)=groupFeat(p,:)-numFeat;
                else
                    groupFeat(p,:)=groupFeat(p,:)+numFeat;
                end
            end
        end
        
        % rename data points
        featNum=unique(groupFeat);
        newFeat=[1:length(featNum)];
        for ff=1:length(featNum),
            groupFeat(groupFeat==featNum(ff))=newFeat(ff);
        end
        
        for f=1:length(featNum),
            if featNum(f)>numFeat,
                featNames(f)={sprintf('ICA%d:neg',featNum(f)-numFeat)};
            else
                featNames(f)={sprintf('ICA%d:pos',featNum(f))};
            end
        end
        
        % map features on group
        Yy=zeros(1,V.dim(1)*V.dim(2)*V.dim(3));
        Yy(1,volIndx)=groupFeat;
        Yy=reshape(Yy,[V.dim(1),V.dim(2),V.dim(3)]);
        Yy(Yy==0)=NaN;
        Vv{1}.dat=Yy;
        Vv{1}.dim=V.dim;
        Vv{1}.mat=V.mat;
        
        % save out vol of ICA feats
        exampleVol=fullfile(studyDir{study},suitDir,'glm4','s02','wdbeta_0001.nii'); % must be better way ??
        X=spm_vol(exampleVol);
        X.fname=outName;
        X.private.dat.fname=V.fname;
        spm_write_vol(X,Vv{1}.dat);
        
        switch metric,
            case 'yes'
                % map vol 2 surf
                M=caret_suit_map2surf(Vv,'space','SUIT','stats','mode','type','paint');
                M.paintnames=featNames;
                M.num_paintnames=length(featNum);
                M.column_name={sprintf('SC%d-%d%%',study,var*100)};
                
                caret_save(fullfile(studyDir{study},caretDir,'suit_flat','glm4',sprintf('funcMap_SC%d_sess%d_%s.paint',study,sess,parcelType)),M);
                
                %         % make areacolour
                cmap=load(fullfile(studyDir{1},encodeDir,'glm4','features.txt'));
                M.column_name=M.paintnames;
                M.column_color_mapping=repmat([-5 5],length(featNum),1);
                M.data=cmap(1:length(featNum),2:4);
                caret_save(fullfile(studyDir{1},caretDir,'suit_flat','glm4','ICA_features.areacolor'),M);
            case 'no'
                % do nothing
        end
    case 'EVAL:visualise_SNN_map'
        sn=varargin{1}; % [2] or 'group'
        study=varargin{2}; % 1 or 2 or [1,2]
        K=varargin{3}; % number of clusters
        type=varargin{4};
        
        % figure out if individual or group
        switch type,
            case 'group'
                outName=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K),'map.nii');
                load(fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K),'SNN.mat'));
            case 'indiv'
                outName=fullfile(studyDir{study},encodeDir,'glm4',subj_name{sn},sprintf('map_SC%d_%dcluster.nii',study,K));
                load(fullfile(studyDir{study},encodeDir,'glm4',subj_name{sn},sprintf('SNN_SC%d_%dcluster.mat',study,K)));
            case 'leaveOneOut'
                outName=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K),sprintf('map_leaveOut_%s.nii',subj_name{sn}));
                load(fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K),sprintf('SNN_leaveOut_%s.mat',subj_name{sn})));
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
    case 'EVAL:funcBound:GROUP'
        study=varargin{1}; 
        mapType=varargin{2}; % look in sc2/encoding/glm4 for map names (options are: 'lob10','lob26','SC1_95POV', 'SC1SC2_95POV' etc)
        data=varargin{3};
        
        % load in map
        mapName=fullfile(studyDir{study},encodeDir,'glm4',mapType,'map.nii');
        
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
        for sn=unique(T.SN)';
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
        end;
        save(fullfile(studyDir{study},'encoding','glm4',mapType,sprintf('spatialBound%s.mat',data)),'-struct','RR');
   
    case 'EVAL:funcBound:INDIV' % NEED TO UPDATE FOR SNN !!
        sn=varargin{1}; % [2]
        study=varargin{2}; % 1 or 2
        var=varargin{3}; % .95
        data=varargin{4}; % 'func1' or 'func2'
        
        outName=fullfile(studyDir{study},encodeDir,'glm4',subj_name{sn},sprintf('SC%d_%dPOV_spatialBound%s.mat',study,var*100,data));
        mapName=fullfile(studyDir{study},encodeDir,'glm4',subj_name{sn},sprintf('map_SC%d_%dPOV.nii',study,var*100));
        
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
    case 'EVAL:funcBound:LEAVEONEOUT' % NEED TO UPDATE FOR SNN !!
        sn=varargin{1}; % [2]
        study=varargin{2}; % 1 or 2
        var=varargin{3}; % .95
        data=varargin{4}; % 'func1' or 'func2
        
        outName=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100),sprintf('spatialBound%s_eval_%s.mat',data,subj_name{sn}));
        mapName=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dPOV',study,var*100),sprintf('map_leaveOut_%s.nii',subj_name{sn}));
        
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
    case 'ICA:PLOT:POV' % NEED TO UPDATE FOR SNN !!
        study=varargin{1}; % 1 or 2
        var=varargin{2}; % {95,90,85,80,75,70,65,60,55,50}
        data=varargin{3}; % 'func1' or 'func2'
        
        P=[];
        for t=1:length(var),
            T=load(fullfile(studyDir{study},'encoding','glm4',sprintf('groupEval_SC%d_%dPOV',study,var{t}*100),sprintf('spatialBound%s.mat',data)));
            
            W=getrow(T,T.crossval==1 & T.bwParcel==0); % crossVal + within
            B=getrow(T,T.crossval==1 & T.bwParcel==1); % crossVal + between
            
            S.wDist=unique(W.dist);
            S.bDist=unique(B.dist);
            for d=1:length(S.wDist),
                w=tapply(W,{'SN','corr'},{'bin','mean','subset',W.dist==S.wDist(d),'name','mBin'});
                b=tapply(B,{'SN','corr'},{'bin','mean','subset',B.dist==S.bDist(d),'name','mBin'});
                [S.t(d),S.p(d)]=ttest(w.corr,b.corr,2,'paired');
            end
            S.t=S.t'; S.p=S.p';
            S.POV=repmat(var{t},length(S.wDist),1);
            P=addstruct(P,S);
            clear S
        end
        
        figure('rend','painters','pos',[20 20 1000 800])
        lineplot(P.wDist,P.t,'split',P.POV,'leg','auto','style_shade')
        axis auto 
    case 'ICA:PLOT:GROUP' % NEED TO UPDATE FOR SNN !!
        study=varargin{1};
        var=varargin{2};
        data=varargin{3};
        type=varargin{4}; % 'group or 'leaveOneOut' or 'Indiv'
        
        switch type,
            case 'group'
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
        subplot(2,1,1);
        xyplot(T.dist,T.corr,T.bin,'split',T.bwParcel,'leg',{'within Parcel','between Parcels'},'subset',T.crossval==0 ,'style_thickline');
        set(gca,'YLim',[0 0.5],'XLim',[0 76]);
        xlabel('Spatial Distance (mm)');
        ylabel('Activity correlation');
        title(sprintf('SC%d(%dPOV)-%s-without crossval',study,var*100,data));
        subplot(2,1,2);
        xyplot(T.dist,T.corr,T.bin,'split',T.bwParcel,'leg',{'within Parcel','between Parcels'},'subset',T.crossval==1,'style_thickline');
        set(gca,'YLim',[0 0.5],'XLim',[0 76]);
        xlabel('Spatial Distance (mm)');
        ylabel('Activity correlation');
        title(sprintf('SC%d(%dPOV)-%s-with crossval',study,var*100,data));
        set(gcf,'PaperPosition',[2 4 10 12]);
        wysiwyg;
    case 'ICA:PLOT:INDIV' % NEED TO UPDATE FOR SNN !!
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
        
    case 'ENCODING:project_taskSpace'
        removeFeats=[1,2,6,12,8];
        
        % load raw Y (cereb)
        load(fullfile(encodeDir,'glm4','cereb_avrgDataStruct_sc1_sc2.mat'));
        
        % load ICs for cereb
        load(fullfile(regDir,'glm4','ICA_comps_cerebellum.mat'));
        
        % get pos + neg IC's
        S_PW=[S_PW;S_PW*-1];
        
        % make model
        Y=nanmean(Y,3); % get group
        [Y_C, ~]=remmean(Y); % center data
        
        % project back to task-space
        Uica=Y_C*S_PW'; % get U * S
        L=Uica'*Uica;
        S=diag(diag(sqrt(L)));
        X=Uica/S;
        Q=setdiff([1:size(S_PW,1)],removeFeats);
        
        % visualise ICc in task-space
        P=caret_load(fullfile(caretDir,'suit_flat','glm4','ICA_features.paint'));
        figure();imagesc_rectangle(X(:,Q),'YDir','reverse');
        caxis([0 1]);
        t=set(gca,'Ytick', 1:length(allCondNames),'YTickLabel',allCondNames,'FontSize',12,'FontWeight','bold',...
            'Xtick',1:length(Q),'XTickLabel',P.paintnames);
        t.Color='white';
        
        T.taskWeights=X(:,Q);
        T.taskNames=allCondNames;
        T.featNames=P.paintnames;
        
        save(fullfile(regDir,'glm4','ICA_taskLoadings.mat'),'T')
        
        % Plot left/right hands
        
        %         % scatterplots
        %         CAT.markertype='o';
        %         CAT.markersize=8;
        %
        %         % plot left/right hands
        %         scatterplot(X(:,14),X(:,10),'label',allCondNames,'CAT',CAT,'regression','linear','intercept',0,'printcorr','draworig')
        %         xlabel('right hand tasks');ylabel('left hand tasks');
    case 'ENCODING:features'
        D=dload(fullfile(rootDir,'featureTable_jd.txt'));        % Read feature table
        S=dload(fullfile(rootDir,'sc1_sc2_taskConds.txt'));   % List of task conditions
        % load(fullfile(sc2Dir,regDir,'glm4','ICA_taskLoadings.mat'));
        load(fullfile(rootDir,'ICA_taskLoadings.mat'));
        W=pivottablerow(S.condNumUni,T.taskWeights,'mean(x,1)');
        
        figure(1);
        subplot(2,2,1);
        imagesc_rectangle(W);
        
        set(gca,'YTick',[1:47],'YTickLabel',D.conditionName,'XTickLabel',T.featNames,'FontSize',6);
        
        % Make the feature matrix
        D.LeftHand    = D.LeftHand ./D.duration;
        D.RightHand   = D.RightHand ./D.duration;
        D.Saccade    = D.Saccade./D.duration;
        f=fieldnames(D);
        FeatureNames = f(5:end);
        F=[];
        for d=1:length(FeatureNames)
            F = [F D.(FeatureNames{d})];
        end;
        F= bsxfun(@rdivide,F,sum(F.^2,1));
        numCond = length(D.conditionName);
        numFeat = length(FeatureNames);
        numICAs = length(T.featNames);
        
        subplot(2,2,2);
        C=corr(F,W);
        imagesc_rectangle(C);
        set(gca,'YTick',[1:numFeat],'YTickLabel',FeatureNames,'XTickLabel',T.featNames,'FontSize',6);
        
        subplot(2,2,[3:4]);
        imagesc_rectangle(F);
        set(gca,'YTick',[1:numCond],'YTickLabel',D.conditionName,'XTick',[1:numFeat],'XTickLabel',FeatureNames,'FontSize',6);
        set(gca,'XTickLabelRotation',60);
        
        % Run a multiple regression analysis on ICAComp onto Features
        lambda = [0.001 0.001];
        X=bsxfun(@minus,F,mean(F,1));
        Y=bsxfun(@minus,W,mean(W,1));
        X=bsxfun(@rdivide,X,sqrt(mean(X.^2)));
        Y=bsxfun(@rdivide,Y,sqrt(mean(Y.^2)));
        XX=X'*X;
        XY=X'*Y;
        A = -eye(numFeat);
        b = zeros(numFeat,1);
        for p=1:numICAs
            U(:,p) = cplexqp(XX+lambda(2)*eye(numFeat),ones(numFeat,1)*lambda(1)-XY(:,p),A,b);
        end;
        
        
        % Present the list of the highest three correlation for each ICA
        for i=1:numICAs
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

