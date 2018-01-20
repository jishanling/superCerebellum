function varargout=sc1_sc2_functionalAtlas(what,varargin)

% Directories
baseDir          = '/Users/maedbhking/Documents/Cerebellum_Cognition';
% baseDir            = '/Volumes/MotorControl/data/super_cerebellum';
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
    case 'TASKSPACE:contrasts'
        
    case 'ATLAS:finalMap'
        % example: 'sc1_sc2_functionalAtlas('ATLAS:finalMap',[2],1,13,'leaveOneOut')
        K=varargin{1}; % number of clusters
        
        sn=returnSubjs;
        outName=fullfile(studyDir{2},encodeDir,'glm4','atlasFinal.nii');
        
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',sn,[1,2]);
        [F,G,Err1]=semiNonNegMatFac(X_C,K,'threshold',0.01);
        save(fullfile(studyDir{2},encodeDir,'glm4','atlasFinal_parcel.mat'),'F','G','volIndx','V');
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
        M.column_name={'atlasFinal'};
        
        caret_save(fullfile(studyDir{2},caretDir,'suit_flat','glm4','atlasFinal.paint'),M);
        
        % make areacolour
        cmap=load(fullfile(studyDir{2},caretDir,'suit_flat','glm4','features.txt'));
        M.column_name=M.paintnames;
        M.column_color_mapping=repmat([-5 5],length(M.column_name),1);
        M.data=cmap(1:length(M.column_name),2:4);
        caret_save(fullfile(studyDir{2},caretDir,'suit_flat','glm4','atlasFinal.areacolor'),M);
        
        figure()
        title('final atlas')
        M=caret_suit_map2surf(Vv,'space','SUIT','stats','mode');
        suit_plotflatmap(M.data,'type','label','cmap',colorcube(K))
        
    case 'EVAL:get_data'
        sn=varargin{1}; % Subje numbers to include
        study=varargin{2}; % 1 or 2 or [1,2]
        
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
            
            % if group - get mean
            UFull=nanmean(UFullAvrg,3);
  
            % remove mean of shared tasks
             UFullAvrg_C=bsxfun(@minus,UFull,mean(UFull(idx,:))); 
            
            % if func1+func2 - concatenate
            if length(study)>1,
                UFullAvrgAll=[UFullAvrgAll;UFullAvrg_C];
            else
                UFullAvrgAll=UFullAvrg_C;
            end
        end

        % center the data (remove overall mean)
        X_C=bsxfun(@minus,UFullAvrgAll,mean(UFullAvrgAll));
        varargout={X_C,volIndx,V};
    case 'EVAL:make_SNN_map'
        % example: 'sc1_sc2_functionalAtlas('EVAL:make_SNN_map',[2],1,13,'leaveOneOut')
        sn=varargin{1}; % subjNum or 'group'
        study=varargin{2}; % 1 or 2 or [1,2]
        K=varargin{3}; % number of clusters
        type=varargin{4}; % 'group','indiv','leaveOneOut'
        
        % figure out if individual or group or leave one out
        switch type,
            case 'group'
                sn=returnSubjs;
                outDir=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K));dircheck(outDir);
                outName=fullfile(outDir,'SNN.mat');
            case 'indiv'
                outDir=fullfile(studyDir{study},encodeDir,'glm4',subj_name{sn});
                outName=fullfile(outDir,sprintf('SNN_SC%d_%dcluster.mat',study,K));
            case 'leaveOneOut'
                outDir=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_SC%d_%dcluster',study,K)); dircheck(outDir);
                outName=fullfile(outDir,sprintf('SNN_leaveOut_%s.mat',subj_name{sn}));
                sn=returnSubjs(returnSubjs~=sn);
        end
        
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',sn,study);
        [F,G,Err1]=semiNonNegMatFac(X_C,K,'threshold',0.01);
        fprintf('SNN map (%d clusters) created for %s \n',K,type)
        save(fullfile(outName),'F','G','volIndx','V');
    case 'EVAL:visualise_SNN_map'
        sn=varargin{1}; % [2] or 'group'
        study=varargin{2}; % 1 or 2 or [1,2]
        K=varargin{3}; % number of clusters
        type=varargin{4}; % 'group','indiv',or 'leaveOneOut'
        
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
        
        figure()
        title(sprintf('%d clusters',K))
        M=caret_suit_map2surf(Vv,'space','SUIT','stats','mode');
        suit_plotflatmap(M.data,'type','label','cmap',colorcube(K))   
    case 'EVAL:bounds:GROUP'
        % code is written now so that func map from sc1 is always evaluated
        % on sc2 data (and vice versa)
        study=varargin{1}; % is map built on study [1] or [2] ?
        mapType=varargin{2}; % options are 'lob10','lob26','bucknerRest','SC<studyNum>_<num>cluster', or 'SC<studyNum>_POV<num>'
        data=varargin{3}; % evaluating data from study [1] or [2] ?
        eval=varargin{4}; % 'unique' or 'all'. Are we evaluating on all taskConds or just those unique to either sc1 or sc2 ?
        
        % load in map
        mapName=fullfile(studyDir{study},encodeDir,'glm4',sprintf('groupEval_%s',mapType),'map.nii');
        
        % load in func data to test (e.g. if map is sc1; func data should
        % be sc2)
        load(fullfile(studyDir{data},'encoding','glm4','cereb_avrgDataStruct.mat'));
        
        switch eval,
            case 'unique'
                % if funcMap - only evaluate unique tasks in sc1 or sc2
                D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
                D1=getrow(D,D.StudyNum==data);
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
        i1=idx;
        i2=idx+length(unique(T.cond)); % these indices are the same for every subj
        for sn=unique(T.SN)';
            %             i1 = find(T.SN==sn & T.sess==1);
            %             i2 = find(T.SN==sn & T.sess==2);
            B=(T.data(i1,voxIn)+T.data(i2,voxIn))/2; % why divide by 2 ?
            
            fprintf('%d cross\n',sn);
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
        save(fullfile(studyDir{study},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d.mat',data)),'-struct','RR');
    case 'EVAL:bounds:INDIV' % NEED TO UPDATE FOR SNN !!
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
    case 'EVAL:bounds:LEAVEONEOUT' % NEED TO UPDATE FOR SNN !!
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
        
    case 'EVAL:PLOT:allMaps'
        study=varargin{1};% is map built on study [1] or [2] ?
        mapType=varargin{2}; % {'lob10','bucknerRest','SC1_5cluster'}
        data=varargin{3}; % evaluating data from study [1] or [2] ?
        
        P=[];
        for m=1:length(mapType),
            T=load(fullfile(studyDir{study},'encoding','glm4',sprintf('groupEval_%s',mapType{m}),sprintf('spatialBoundfunc%d.mat',data)));
            
            A=getrow(T,T.crossval==1);
            A.type=repmat({mapType{m}},length(A.SN),1);
            A.m=repmat(m,length(A.SN),1);
            P=addstruct(P,A);
            clear A
        end
        
        % plot correlations (across clusters) for within
        figure('name','within')
        xyplot(P.dist,P.corr,P.bin,'split',P.type,'leg','auto','subset',P.bwParcel==0,'style_thickline');
        
        % plot correlations (across clusters) for between
        figure('name','between')
        xyplot(P.dist,P.corr,P.bin,'split',P.type,'leg','auto','subset',P.bwParcel==1,'style_thickline');
        
        % plot boxplot of different clusters
        W=getrow(P,P.bwParcel==0); % within
        B=getrow(P,P.bwParcel==1); % between
        W.diff=W.corr-B.corr;
        
        figure()
        myboxplot(W.m,W.diff) % within-between diff
        figure()
        lineplot(W.m,W.diff) % within-between diff
    case 'EVAL:PLOT:GROUP'
        study=varargin{1};% is map built on study [1] or [2] ?
        mapType=varargin{2}; % options are 'lob10','lob26','bucknerRest','SC<studyNum>_<num>cluster', or 'SC<studyNum>_POV<num>'
        data=varargin{3}; % evaluating data from study [1] or [2] ?
        type=varargin{4}; % 'group' or 'leaveOneOut'
        
        switch type,
            case 'group'
                T=load(fullfile(studyDir{study},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d.mat',data)));
            case 'leaveOneOut'
                T=[];
                for s=1:length(returnSubjs),
                    P=load(fullfile(studyDir{study},'encoding','glm4',sprintf('groupEval_%s',mapType),sprintf('spatialBoundfunc%d_eval_%s.mat',data,subj_name{returnSubjs(s)})));
                    T=addstruct(T,P);
                end
            otherwise
                fprintf('no such case')
        end
        
        figure()
        subplot(2,1,1);
        xyplot(T.dist,T.corr,T.bin,'split',T.bwParcel,'leg',{'within Parcel','between Parcels'},'subset',T.crossval==0 ,'style_thickline');
        set(gca,'YLim',[0 0.5],'XLim',[0 76]);
        xlabel('Spatial Distance (mm)');
        ylabel('Activity correlation');
        title(sprintf('%s-func%d-without crossval',mapType,data));
        subplot(2,1,2);
        xyplot(T.dist,T.corr,T.bin,'split',T.bwParcel,'leg',{'within Parcel','between Parcels'},'subset',T.crossval==1,'style_thickline');
        set(gca,'YLim',[0 0.5],'XLim',[0 76]);
        xlabel('Spatial Distance (mm)');
        ylabel('Activity correlation');
        title(sprintf('%s-func%d-with crossval',mapType,data));
        set(gcf,'PaperPosition',[2 4 10 12]);
        wysiwyg;
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
        
    case 'ENCODE:project_taskSpace'
        
        % load raw Y (cereb)
        [Y,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',returnSubjs,[1,2]);
       
        % load clusters for cereb
        load(fullfile(studyDir{2},encodeDir,'glm4','atlasFinal_parcel.mat'));
        
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
        D=dload(fullfile(rootDir,'featureTable_jd.txt'));        % Read feature table
        S=dload(fullfile(rootDir,'sc1_sc2_taskConds.txt'));      % List of task conditions
        load(fullfile(studyDir{2},encodeDir,'glm4','atlasFinal_condLoads.mat'));
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
        numClusters = length(T.featNames);
        
        subplot(2,2,2);
        C=corr(F,W);
        imagesc_rectangle(C);
        set(gca,'YTick',[1:numFeat],'YTickLabel',FeatureNames,'XTickLabel',T.featNames,'FontSize',6);
        
        subplot(2,2,[3:4]);
        imagesc_rectangle(F);
        set(gca,'YTick',[1:numCond],'YTickLabel',D.conditionName,'XTick',[1:numFeat],'XTickLabel',FeatureNames,'FontSize',6);
        set(gca,'XTickLabelRotation',60);
        
        % Run a multiple regression analysis on clustesr onto features
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
            U(:,p) = cplexqp(XX+lambda(2)*eye(numFeat),ones(numFeat,1)*lambda(1)-XY(:,p),A,b);
        end;

        % Present the list of the highest three correlation for each
        % cluster
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

