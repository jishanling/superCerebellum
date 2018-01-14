function varargout=sc1sc2_spatialAnalysis(what,varargin)
% Does some specific spatial analyses on the cerebellum - to be added later
% to the final analysis routine
baseDir = '/Users/jdiedrichsen/Data/super_cerebellum_new';
studyDir{1}     =fullfile(baseDir,'sc1');
studyDir{2}     =fullfile(baseDir,'sc2');
suitDir         ='/suit';
caretDir        ='/surfaceCaret';
regDir          ='/RegionOfInterest/';
encodeDir       ='/encoding';
% rootDir = '/Users/jdiedrichsen/Projects/SuperCerebellum/';

returnSubjs=[2,3,4,6,7,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31];

switch(what)
    case 'avrgDataStruct_freq'
        frequencyBands  = [0 0.5 1 1.5 2 inf];
        load(fullfile(baseDir,'sc2','encoding','glm4','cereb_avrgDataStruct.mat'));
        RR=[];
        sn=unique(T.SN);
        for s = 1:length(sn)
            fprintf('subject %d\n',sn(s));
            for se=1:2
                S=getrow(T,T.SN == sn(s) & T.sess==se);
                for c=1:29
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
        save(fullfile(baseDir,'sc1','encoding','glm4','cereb_avrgDataStruct_freq.mat'),'-struct','RR');
    case 'intersubjectCorrSpatial'
        T=load(fullfile(baseDir,'sc1','encoding','glm4','cereb_avrgDataStruct_freq.mat'));
        D=load(fullfile(baseDir,'sc1','encoding','glm4','cereb_avrgDataStruct.mat'));
        D.T.freq = ones(length(D.T.SN),1)*0;
        D.T.freqLow = -1*ones(length(D.T.SN),1);
        T=addstruct(T,D.T);
        D=[];
        sn=unique(T.SN);
        numSubj = length(sn);
        RR=[];
        for f=unique(T.freq)'
            for s = 1:numSubj
                for se=1:2
                    D(:,:,s+(se-1)*length(sn))=T.data(T.SN==sn(s) & T.sess==se & T.freq==f,:);
                end;
            end;
            C=intersubj_corr(D);
            R.sess = [ones(1,numSubj) ones(1,numSubj)*2]';
            R.subj = [1:numSubj 1:numSubj]';
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
        xlabels={'0-0.5','0.5-1','1-1.5','1.5-2','>2'}
        T=varargin{1}; % Input from case 'intersubjectCorrSpatial'
        T=tapply(T,{'subj','freq'},{'withinSubj'},{'betweenSubj'},{'totSS'});
        [ss,sn]=pivottable(T.subj,[],T.totSS,'mean','subset',T.freq==0);
        a(sn,1)=ss;
        T.relSS=T.totSS./a(T.subj);
        subplot(2,1,1);
        lineplot([T.freq],[T.relSS],'subset',T.freq>0,'style_thickline');
        set(gca,'XTickLabel',xlabels);
        ylabel('Relative Power');
        xlabel('Cycles/cm');
        title('Relative amount of Variance per Frequency band');
        subplot(2,1,2);
        lineplot([T.freq],[T.withinSubj T.betweenSubj],'subset',T.freq>0,'style_thickline','leg',{'within','between'});
        a=mean(T.withinSubj(T.freq==0));
        b=mean(T.betweenSubj(T.freq==0));
        drawline([0 a b],'dir','horz');
        ylabel('Correlationa');
        xlabel('Cycles/cm');
        set(gca,'XTickLabel',xlabels);
        title('Within and between-subject correlation');
        set(gcf,'PaperPosition',[2 2 4 5]);
        wysiwyg;
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
        sn=varargin{1}; % subjNum or 'group'
        study=varargin{2}; % 1 or 2 or [1,2]
        type=varargin{3}; % 'group'
        K=varargin{4};       % Number of clusters
        % figure out if individual or group or leave one out
        outDir=fullfile(studyDir{study},encodeDir,'glm4','SNN');
        switch type,
            case 'group'
                sn=returnSubjs;
                outName=fullfile(outDir,sprintf('SNN_%d.mat',K));
            case 'indiv'
                outName=fullfile(outDir,sprintf('SNN_%d_indiv_%s.mat',K,subj_name{sn}));
            case 'leaveOneOut'
                outName=fullfile(outDir,sprintf('SNN_%d_loo_%s.mat',K,subj_name{sn}));
                sn=returnSubjs(returnSubjs~=sn);
        end
        
        [X_C,volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',sn,study);
        [F,G,Err1]=semiNonNegMatFac(X_C,K,'threshold',0.01);
        save(fullfile(outName),'F','G','volIndx','V');
    case 'SNN:visualise_map'
        sn=varargin{1}; % subjNum or 'group'
        study=varargin{2}; % 1 or 2 or [1,2]
        type=varargin{3};
        K=varargin{4};       % Number of clusters
        
        % figure out if individual or group
        outDir=fullfile(studyDir{study},encodeDir,'glm4','SNN');
        switch type,
            case 'group'
                sn=returnSubjs;
                outName=fullfile(outDir,sprintf('SNN_%d.mat',K));
            case 'indiv'
                outName=fullfile(outDir,sprintf('SNN_%d_indiv_%s.mat',K,subj_name{sn}));
            case 'leaveOneOut'
                outName=fullfile(outDir,sprintf('SNN_%d_loo_%s.mat',K,subj_name{sn}));
        end
        
        load(outName);
        [~,groupFeat]=max(G,[],2);
        
        % map features on group
        Yy=zeros([V.dim(1) V.dim(2) V.dim(3)]);
        Yy(volIndx)=groupFeat;  % JD: No reshaping necessary
        Yy(Yy==0)=NaN;
        Vv{1}.dat=Yy;
        Vv{1}.dim=V.dim;
        Vv{1}.mat=V.mat;
        
        % save out vol of ICA feats
        M=caret_suit_map2surf(Vv,'space','SUIT','stats','mode');
        suit_plotflatmap(M.data,'type','label','cmap',colorcube(K))
        
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