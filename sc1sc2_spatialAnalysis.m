function varargout=sc1sc2_spatialAnalysis(what,varargin)
% Does some specific spatial analyses on the cerebellum - to be added later
% to the final analysis routine
baseDir = '/Users/jdiedrichsen/Data/super_cerebellum_new';
baseDir = '/Volumes/MotorControl/Data/super_cerebellum_new';

studyDir{1}     =fullfile(baseDir,'sc1');
studyDir{2}     =fullfile(baseDir,'sc2');
suitDir         ='/suit';
caretDir        ='/surfaceCaret';
regDir          ='/RegionOfInterest/';
encodeDir       ='/encoding';
% rootDir = '/Users/jdiedrichsen/Projects/SuperCerebellum/';

returnSubjs=[2,3,4,6,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31];

switch(what)
    case 'avrgDataStruct_freq'
        exper=varargin{1};
        frequencyBands  = [0 0.5 1 1.5 2 inf];
        load(fullfile(studyDir{exper},'encoding','glm4','cereb_avrgDataStruct.mat'));
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
        save(fullfile(studyDir{exper},'encoding','glm4','cereb_avrgDataStruct_freq.mat'),'-struct','RR');
    case 'intersubjectCorrSpatial'
        exper = [1 2]; % %experiment
        glm   = 'glm4';
        vararginoptions(varargin,{'exper','glm'});
        C=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        A=[];
        for e=exper
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
        varargout={RR};
    case 'intersubjectCorrPlot'
        xlabels={'overall','0-0.5','0.5-1','1-1.5','1.5-2','>2'}
        T=varargin{1}; % Input from case 'intersubjectCorrSpatial'
        T=tapply(T,{'subj','freq'},{'withinSubj'},{'betweenSubj'},{'totSS'},'subset',ismember(T.subj,returnSubjs));
        T.freqK = T.freq>0;
        [ss,sn]=pivottable(T.subj,[],T.totSS,'mean','subset',T.freq==0);
        a(sn,1)=ss;
        T.relSS=T.totSS./a(T.subj);
        subplot(2,1,1);
        lineplot([T.freqK T.freq],[T.relSS],'style_thickline');
        set(gca,'XTickLabel',xlabels,'YLim',[0 0.35]);
        ylabel('Relative Power');
        xlabel('Cycles/cm');
        title('Relative amount of Variance per Frequency band');
        subplot(2,1,2);
        lineplot([T.freqK T.freq],[T.withinSubj T.betweenSubj],'style_thickline','leg',{'within','between'});
        ylabel('Correlationa');
        xlabel('Cycles/cm');
        drawline(0,'dir','horz');
        set(gca,'XTickLabel',xlabels);
        title('Within and between-subject correlation');
        set(gcf,'PaperPosition',[2 2 4.5 5.5]);
        wysiwyg;
    case 'freqPlot'
        exper= varargin{1};
        sn   = [26 26 27];
        se   = [1 2 1];
        con  = varargin{2};
        
        frequencyBands  = [0 0.5 1.5 2 inf];
        load(fullfile(studyDir{exper},'encoding','glm4','cereb_avrgDataStruct.mat'));
        RR=[];
        sn=unique(T.SN);
        for i=1:3
            S=getrow(T,T.SN == sn(i) & T.sess==se(i) & T.cond==con);
            X=zeros(V.dim);
            X(volIndx)=S.data;
            X(isnan(X))=0;
            % Y=mva_frequency3D(X,frequencyBands,'Voxelsize',[2 2 2],'plotSlice',15);
            figure(i);
            set(gcf,'PaperPosition',[2 2 4 12]);
            wysiwyg;
            Y=mva_frequency3D(X,frequencyBands,'Voxelsize',[2 2 2],'plotSlice',22,'numRows',5);
        end;
    case 'functionalBounds'
        suitDir=fileparts(which('spm_suit'));
        functionalParcel = fullfile(suitDir,'atlas','Cerebellum-SUIT.nii');
        load(fullfile(baseDir,'sc1','encoding','glm4','cereb_avrgDataStruct.mat'));
        % Now get the parcellation sampled into the same space
        [i,j,k]=ind2sub(V.dim,volIndx);
        [x,y,z]=spmj_affine_transform(i,j,k,V.mat);
        VA= spm_vol(functionalParcel);
        [i1,j1,k1]=spmj_affine_transform(x,y,z,inv(VA.mat));
        Parcel = spm_sample_vol(VA,i1,j1,k1,0);
        % Devide the voxel pairs into all the spatial bins that we want
        fprintf('parels\n');
        voxIn = Parcel>0;
        XYZ= [x;y;z];
        RR=[];
        [BIN,R]=mva_spatialCorrBin(XYZ(:,voxIn),'Parcel',Parcel(1,voxIn));
        clear XYZ i k l x y z i1 j1 k1 VA Parcel; % Free memory
        % Now calucalte the uncrossvalidated and crossvalidated
        % Estiamte of the correlation for each subject
        for sn=unique(T.SN)'
            i1 = find(T.SN==sn & T.sess==1);
            i2 = find(T.SN==sn & T.sess==2);
            D=(T.data(i1,voxIn)+T.data(i2,voxIn))/2;
            
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
        save(fullfile(baseDir,'sc1','encoding','glm4','cerebSpatialBoundLob.mat'),'-struct','RR');
        varargout={RR};
    case 'functionalBoundsPlot'
        T=varargin{1}; % REsults from functionalBounds
        subplot(2,1,1);
        xyplot(T.dist,T.corr,T.bin,'split',T.bwParcel,'leg',{'within Parcel','between Parcels'},'subset',T.crossval==0 ,'style_thickline');
        set(gca,'YLim',[0 0.5],'XLim',[0 76]);
        xlabel('Spatial Distance (mm)');
        ylabel('Activity correlation');
        title('without crossvalidation');
        subplot(2,1,2);
        xyplot(T.dist,T.corr,T.bin,'split',T.bwParcel,'leg',{'within Parcel','between Parcels'},'subset',T.crossval==1,'style_thickline');
        set(gca,'YLim',[0 0.5],'XLim',[0 76]);
        xlabel('Spatial Distance (mm)');
        ylabel('Activity correlation');
        title('with crossvalidation');
        set(gcf,'PaperPosition',[2 2 4 8]);
        wysiwyg;
    case 'SNN:simulate'
        P=100; % Number of voxels
        N=42; % Number of tasks
        Ktrue=10; % Number of clusters
        F=normrnd(0,1,N,Ktrue); % Cluster means
        G=rand(P,Ktrue); % Cluster weights: uniform distribution
        
        % Sort G
        [~,a]=max(G,[],2);
        [~,i]=sort(a);
        G=G(i,:);
        
        X=F*G'+normrnd(0,0.01,N,P);
        X=bsxfun(@minus,X,mean(X));
        
        % Compare the different algorithms
        K = 10;
        g = kmeans(X', K);
        Gkm = indicatorMatrix('identity',g);
        Fkm = X*pinv(G');
        
        % Semi-nonnegative MF:
        [Fsn,Gsn,Err1]=semiNonNegMatFac(X,K,'G0',Gkm,'threshold',0.01);
        
        % convex semi-nonnegative MF
        [Fcn,Gcn,Err2,W]=cnvSemiNonNegMatFac(X,K,'G0',Gkm,'threshold',0.01);
        
        subplot(2,4,1);
        imagesc(G);
        title('truth');
        subplot(2,4,2);
        imagesc(Gkm);
        title('kmeans');
        subplot(2,4,3);
        imagesc(Gsn);
        title('semi nonneg');
        subplot(2,4,4);
        imagesc(Gcn);
        title('convex sn');
        
        subplot(2,4,5);
        imagesc(F);
        subplot(2,4,6);
        imagesc(Fkm);
        subplot(2,4,7);
        imagesc(Fsn);
        subplot(2,4,8);
        imagesc(Fcn);
        Err1(end)
        Err2(end)
    case 'SNN:make_map'
        % example: 'sc1_sc2_functionalAtlas('ICA:make_map',[2],1,.95)
        sn=varargin{1}; % subject numbers
        study=varargin{2}; % 1 or 2 or [1,2]
        K=varargin{3};       % Number of clusters
        % figure out if individual or group or leave one out
        outDir=fullfile(studyDir{study},encodeDir,'glm4','SNN');
        outName=fullfile(outDir,sprintf('SNN_%d.mat',K));
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',sn,study,'build');
        [F,G,Info]=semiNonNegMatFac(X_C,K,'threshold',0.01);
        save(fullfile(outName),'F','G','volIndx','V');
    case 'SNN:check_convergence'
        % example: 'sc1_sc2_functionalAtlas('ICA:make_map',[2],1,.95)
        sn=varargin{1}; % subject numbers
        study=varargin{2}; % 1 or 2 or [1,2]
        K=varargin{3};       % Number of clusters
        D=[];
        numFits=10;
        % figure out if individual or group or leave one out
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',sn,study,'build');
        for i=1:numFits
            fprintf('%d\n',i);
            [F,G,Info]=semiNonNegMatFac(X_C,K,'threshold',0.01);
            [~,Cl(:,i)]=max(G,[],2);
            D=addstruct(D,Info);
        end;
        for i=1:numFits
            for j=1:numFits
                fprintf('%d,%d\n',i,j);
                RI(i,j)=RandIndex(Cl(:,i),Cl(:,j));
            end;
        end;
        subplot(3,1,[1:2])
        imagesc(RI)
        colorbar
        subplot(3,1,3)
        plot([1:10]',D.error);
        keyboard;
    case 'SNN:map_weights'
        % this function takes any labelled volume (already in SUIT space)
        % and plots to the surface
        inputMap=varargin{1}; % some options are 'Buckner_7Networks','SC1_9cluster','lob10', 'Cole_10Networks', 'SC2_90cluster' etc
        
        vararginoptions(varargin(2:end),{''});
        inputDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s',inputMap));
        load(fullfile(inputDir,'SNN.mat'));
        Vv{1}.dim=V.dim;
        Vv{1}.mat=V.mat;
        Vv{1}.dat=zeros(V.dim);
        G = bsxfun(@rdivide,bestG,sum(bestG,2));
        numComp = size(bestG,2);
        for i=1:numComp
            subplot(3,ceil(numComp/3),i);
            Vv{1}.dat(volIndx)=G(:,i);
            M=caret_suit_map2surf(Vv,'space','SUIT','stats','mode');
            suit_plotflatmap(M.data,'type','func','border',[],'cscale',[0 0.4])
        end;
    case 'SNN:bootstrap'
        type = 'oldF';      % {'oldF','newF'}: Estimate a new F every bootstrap interation?
        study=[1 2];        % Study 1 or 2 or [1,2]
        K=10;               % K=numClusters (i.e. 5);
        numBSIter = 100;     % Number of Bootstrap iterations
        algorithm = 'cnvf'; % Semi-nonengative mare
        smooth = [];  % Any postfix to the file name? (smoothing, etc)
        vararginoptions(varargin,{'type','study','K','numBSIter','algorithm','smooth'});
        
        % Set the String correctly
        studyStr = sprintf('SC%d',study);
        if length(study)>1
            studyStr='SC12'; % both studies combined
        end
        outDir=fullfile(studyDir{2},encodeDir,'glm4',sprintf('groupEval_%s_%s_%d',studyStr,algorithm,K));
        if ~isempty(smooth)
            outDir=strcat(outDir,sprintf('_s%2.1f',smooth));
        end;
        
        S = load(fullfile(outDir,sprintf('%s.mat',algorithm)));
        N =length(returnSubjs);
        
        
        for i=1:numBSIter
            % sample with replacement
            if (i==1)
                T.samp(i,:)=returnSubjs;
            else
                T.samp(i,:)=returnSubjs(unidrnd(N,1,N));
            end;
            [X,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',T.samp(i,:),study,'build','smooth',smooth);
            switch (type)
                case 'newF'
                    switch (algorithm)
                        case 'snn'
                            [F,G,Info]=semiNonNegMatFac(X,K,'G0',S.bestG,'threshold',0.001); % get a segmentation using
                        case 'cnvf'
                            [F,G,Info]=cnvSemiNonNegMatFac(X,K,'G0',S.bestG,'threshold',0.01,'maxIter',100); % get a segmentation using
                    end;
                    keyboard;
                case 'oldF'
                    [M,P]=size(X);
                    G=zeros(P,K);
                    for p=1:P
                        G(p,:) = lsqnonneg(double(S.bestF),X(:,p));
                    end;
                    R = X-double(S.bestF)*G';
                    T.Error(i,1)=sum(sum(R.*R));
                    [~,g]=max(G,[],2);
                    T.assign(i,:)=g';
                    fprintf('.');
            end;
        end;
        fprintf('\n');
        save(fullfile(outDir,'bootstrap_oldF.mat'),'-struct','T');
        varargout={T};
    case 'SNN:bootstrap_eval'
        mapsize = [ 2 2 5 4]; % For paper size
        T=load('bootstrap_oldF.mat');
        [~,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',2,1,'build'); % Get V and volindx
        N=size(T.Error,1);
        for n=2:N
            T.RandIndx(n,1)=RandIndex(T.assign(1,:)',T.assign(n,:)');
        end;
        % Generate map of assignments
        CON=bsxfun(@eq,T.assign(2:N,:),T.assign(1,:)); % Consistency of assignment
        con = mean(CON);
        
        
        % This is the histogram of Rand coefficients
        figure;
        histplot(T.RandIndx(2:end),'numcat',10);
        set(gcf,'PaperPosition',[2 2 5 3]);
        wysiwyg;
        
        figure;
        % If cmap is file, load and normalize
        cmap=dlmread('colourMap.txt');
        cmap = cmap(:,2:end);
        cmap = bsxfun(@rdivide,cmap,max(cmap));
        
        % Figure 2: Normal hard assignment
        LA=sc1sc2_spatialAnalysis('visualise_map',T.assign(1,:)',volIndx,V,'type','label','cmap',cmap);
        % title('Hard assignment');
        axis off;
        set(gcf,'PaperPosition',mapsize);
        wysiwyg;
        
        % Figure 3: Consistency
        figure;
        CO=sc1sc2_spatialAnalysis('visualise_map',con',volIndx,V,'type','func','cscale',[0.4 1]);
        % title('Consistency of assignment');
        axis off;
        set(gcf,'PaperPosition',mapsize);
        wysiwyg;
        
        % Figure 4: Certainty is expressed as a morph to gray (saturation)
        figure;
        i = ~isnan(LA.data);
        RGB = nan(size(LA.data,1),3);
        HSV = nan(size(LA.data,1),3);
        
        RGB(i,:)=cmap(LA.data(i),:); % color data
        HSV(i,:) = rgb2hsv(RGB(i,:)); % color data
        minCO = 0.5; % Setting the lower bound of consistency that will be totally gray
        scaleC = (CO.data(i,1)-minCO)./(1-minCO);  %
        scaleC(scaleC<0)=0;     % scaling factor between gray and color
        HSV(i,2) = HSV(i,2).*scaleC;                % Scale Saturation
        HSV(i,3) = (1-scaleC).*0.85+scaleC.*HSV(i,3); % Change towards a light gray
        mRGB(i,:) = hsv2rgb(HSV(i,:)); % color data
        suit_plotflatmap(mRGB,'type','rgb');
        % title('Weighted color map');
        axis off;
        set(gcf,'PaperPosition',mapsize);
        wysiwyg;
        
        
        % figure(5); % Different Alphas: We decided that this did not look
        % as good.
        % scaleA=nan(size(RGB,1),1);
        % scaleA(i,1)=scaleC;
        % suit_plotflatmap(RGB,'type','rgb','alpha',scaleA);
     
    case 'CLUSTER:GlobalRand'
        compare={'groupEval_SC12_cnvf_7','groupEval_SC12_cnvf_10','groupEval_SC12_cnvf_17',...
            'groupEval_Buckner_7Networks','groupEval_Cole_10Networks','groupEval_Buckner_17Networks'};
        numMaps = length(compare);
        load(fullfile(studyDir{2},encodeDir,'glm4','cereb_avrgDataStruct.mat'));  % Just to get V and volIndex
        clear T;
        for i=1:numMaps
            try
                T=load(fullfile(baseDir,'sc2','encoding','glm4',compare{i},'cnvf.mat'));
                [~,c(:,i)]=max(T.bestG,[],2);
            catch
                Vi=spm_vol(fullfile(baseDir,'sc2','encoding','glm4',compare{i},'map.nii'));
                [i1,j1,k1]=ind2sub(V.dim,volIndx');
                [x,y,z]=spmj_affine_transform(i1,j1,k1,V.mat);
                [i2,j2,k2]=spmj_affine_transform(x,y,z,inv(Vi.mat));
                c(:,i)=spm_sample_vol(Vi,i2,j2,k2,0);
            end;
        end;
        for i=1:numMaps-1
            for j=i+1:numMaps
                AR(i,j)=RandIndex(c(:,i),c(:,j));
                AR(j,i)=AR(i,j);
            end;
        end;
        varargout={AR};
        colormap(hot);
        imagesc(AR,[0 0.8]);
        colorbar;
        set(gca,'XTickLabel',{'C7','C10','C17','R7','R10','R17'},'YTickLabel',{'C7','C10','C17','R7','R10','R17'}); 
    case 'CLUSTER:LocalRand'
        type = 'searchlight'; % 'voxelwise','searchlight', or 'searchlightvoxel'
        compare={'groupEval_SC12_cnvf_7','groupEval_SC12_cnvf_10','groupEval_SC12_cnvf_17',...
            'groupEval_Buckner_7Networks','groupEval_Cole_10Networks','groupEval_Buckner_17Networks'...
            
            % ,'groupEval_SC12_9cluster',...
            % 'groupEval_SC12_10cluster','groupEval_SC12_11cluster','groupEval_SC12_12cluster',...
            % 'groupEval_SC1_7cluster','groupEval_SC1_8cluster','groupEval_SC1_9cluster',...
            % 'groupEval_SC1_10cluster','groupEval_SC1_11cluster','groupEval_SC1_12cluster',...
            % 'groupEval_SC2_7cluster','groupEval_SC2_8cluster','groupEval_SC2_9cluster',...
            % 'groupEval_SC2_10cluster','groupEval_SC2_11cluster','groupEval_SC2_12cluster',...
            };
        split = [1 1 1 2 2 2];
        radius = 10;
        vararginoptions(varargin,{'radius','compare','type','split'});
        numMaps = length(compare);
        load(fullfile(studyDir{2},encodeDir,'glm4','cereb_avrgDataStruct.mat'));  % Just to get V and volIndex
        clear T;
        
        for i=1:numMaps
            try
                T=load(fullfile(baseDir,'sc2','encoding','glm4',compare{i},'cnvf.mat'));
                [~,c(:,i)]=max(T.bestG,[],2);
            catch
                Vi=spm_vol(fullfile(baseDir,'sc2','encoding','glm4',compare{i},'map.nii'));
                [i1,j1,k1]=ind2sub(V.dim,volIndx');
                [x,y,z]=spmj_affine_transform(i1,j1,k1,V.mat);
                [i2,j2,k2]=spmj_affine_transform(x,y,z,inv(Vi.mat));
                c(:,i)=spm_sample_vol(Vi,i2,j2,k2,0);
            end;
        end;
        ar=[];
        [x,y,z]=ind2sub(T.V.dim,T.volIndx);
        [x,y,z]=spmj_affine_transform(x,y,z,T.V.mat);
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
            sc1sc2_spatialAnalysis('visualise_map',mean(ar(:,pair==i),2),T.volIndx,T.V,'type','func');
        end;
        Output.type= type;
        Output.compare= compare;
        Output.split= split;
        Output.ar= ar;
        Output.pair= pair;
        varargout = {Output};
    case 'CLUSTER:FigureRand' % Figure of the Rand index
        %Out=sc1sc2_spatialAnalysis('CLUSTER:LocalRand');
        % save(fullfile(studyDir{2},encodeDir,'glm4','LocalRand_TaskRest.mat'),'-struct','Out');
        set(gcf,'PaperPosition',[2 2 8 6]);
        wysiwyg;
        load(fullfile(studyDir{2},encodeDir,'glm4','LocalRand_TaskRest.mat'));
        load(fullfile(studyDir{2},encodeDir,'glm4','cereb_avrgDataStruct.mat'));  % Just to get V and volIndex
        clear T;
        subplot(2,2,1);
        sc1sc2_spatialAnalysis('CLUSTER:GlobalRand');
        
        titles = {'Task-Task','Rest-Rest','Task-Rest'};
        for i=1:3
            subplot(2,2,i+1);
            sc1sc2_spatialAnalysis('visualise_map',mean(ar(:,pair==i),2),volIndx,V,'type','func','cscale',[0 0.8]);
            axis equal;
            title(titles{i});
        end;
        
    case 'visualise_map' % Plot any data on the cerebellar flatmap
        data = varargin{1};     % Data to plot
        volIndx = varargin{2};  % Indices into the volume (mask)
        V = varargin{3};        % Cerebellar suit volume
        cmap = [];
        cscale = [];
        type = 'label';     % func / label
        vararginoptions(varargin(4:end),{'cmap','type','cscale'});
        
        % map features on group
        V.dat=zeros([V.dim(1) V.dim(2) V.dim(3)]);
        V.dat(volIndx)=data;
        
        switch (type)
            case 'label'
                stats = 'mode';
                if (isempty(cmap))
                    cmap = colorcube(max(data));
                end;
            case {'func','rgb'}
                stats = 'nanmean';
                if (isempty(cmap))
                    cmap = hot;
                end;
        end;
        
        % map the data and display on flatmap
        data=suit_map2surf(V,'space','SUIT','stats',stats);
        suit_plotflatmap(data,'type',type,'cmap',cmap,'cscale',cscale);
        
        P=caret_struct('paint','data',data);
        varargout={P};
    case 'EVALBOUND:getBoundaries'
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
        numParcel = 28; % max(Parcel);
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
    case 'EVALBOUND:evalBoundaries'
        mapType = varargin{1};
        study   = varargin{2}; % evaluating data from study [1] or [2] ?
        
        spatialBins = [0:3:20];
        condType    = 'all';
        
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
        
    case 'EVALBOUND:fullEval'   % Evaluate a parcellation on both studies and save
        mapType = varargin{1};
        EvalDir = fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_%s',mapType));
        
        T1=sc1sc2_spatialAnalysis('EVALBOUND:evalBoundaries',mapType,1);
        T2=sc1sc2_spatialAnalysis('EVALBOUND:evalBoundaries',mapType,2);
        T1.study = ones(length(T1.SN),1)*1;
        T2.study = ones(length(T1.SN),1)*2;
        T=addstruct(T1,T2);
        save(fullfile(EvalDir,'BoundariesFunc3_all.mat'),'-struct','T');
        varargout={T};
    case 'EVALBOUND:visualize' % Generates metric file and border file for the pacellation
        mapType = varargin{1};
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
    case 'EVALBOUND:unCrossval:GROUP'
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
    case 'EVALBOUND:Figure'
        set(gcf,'position',[10 10 1600 800]);
        
        subplot(2,3,1);
        sc1_sc2_functionalAtlas('STRENGTH:visualise_bound','lob10');
        set(gca,'XLim',[-100 100],'YLim',[-100 100],'Color',[0 0 0],'Visible','on','Xticklabel',[],'Yticklabel',[],'Box','off');
        axis equal;
        subplot(2,3,2);
        sc1_sc2_functionalAtlas('STRENGTH:visualise_bound','Buckner_17Networks');
        set(gca,'XLim',[-100 100],'YLim',[-100 100],'Color',[0 0 0],'Visible','on','Xticklabel',[],'Yticklabel',[],'Box','off');
        axis equal;
        subplot(2,3,3);
        sc1_sc2_functionalAtlas('STRENGTH:visualise_bound','Buckner_7Networks');
        set(gca,'XLim',[-100 100],'YLim',[-100 100],'Color',[0 0 0],'Visible','on','Xticklabel',[],'Yticklabel',[],'Box','off');
        axis equal;
        subplot(2,3,4);
        sc1_sc2_functionalAtlas('STRENGTH:visualise_bound','SC12_10cluster');
        set(gca,'XLim',[-100 100],'YLim',[-100 100],'Color',[0 0 0],'Visible','on','Xticklabel',[],'Yticklabel',[],'Box','off');
        axis equal;
end;
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

% Within/Between Subj/Sess Corr
function [DiffSubjSameSess,SameSubjDiffSess,DiffSubjDiffSess]=bwse_corr(C)
N=size(C,1);
n=N/2;
sess = [ones(1,n) ones(1,n)*2];
subj = [1:n 1:n];
SameSess = bsxfun(@eq,sess',sess);
SameSubj = bsxfun(@eq,subj',subj);
DiffSubjSameSess=mean(C(~SameSubj & SameSess));
SameSubjDiffSess=mean(C(SameSubj & ~SameSess));
DiffSubjDiffSess=mean(C(~SameSubj & ~SameSess));
%                 ROI