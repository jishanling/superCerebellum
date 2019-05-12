function varargout=sc1_sc2_GroupIndivid(what,varargin)

% Directories
% baseDir          = '/Users/maedbhking/Remote/Documents2/Cerebellum_Cognition';
baseDir            = '/Volumes/MotorControl/data/super_cerebellum_new';
% baseDir          = '/Users/jdiedrichsen/Data/super_cerebellum_new';

atlasDir='/Users/maedbhking/Documents/Atlas_templates/';

studyDir{1}     =fullfile(baseDir,'sc1');
studyDir{2}     =fullfile(baseDir,'sc2');
studyStr        = {'SC1','SC2','SC12'};
behavDir        ='/data';
suitDir         ='/suit';
caretDir        ='/surfaceCaret';
regDir          ='/RegionOfInterest/';
encodeDir       ='/encoding';

subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11',...
    's12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24',...
    's25','s26','s27','s28','s29','s30','s31'};

returnSubjs=[2,3,4,6,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31];

hem={'lh','rh'};
hemName={'LeftHem','RightHem'};

switch what    
    case 'BuildAndEval'    % Uses different amounts of data and different weightings of group and individual data to build and evaluate individual parcellations 
        mapName = 'cnvf_10';  % Use the features that have been learned over both experiments
        subjs=returnSubjs;
        groupW = 0; 
        vararginoptions(varargin,{'mapName','subjs','groupW'}); 
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        % load group features (these don't need to be crossvalidated)
        for exp=1:2
            GP=load(fullfile(studyDir{2},'encoding','glm4',sprintf('groupEval_SC%d_%s',exp,mapName),'cnvf.mat'));
            groupF{exp}=bsxfun(@minus,GP.bestF,mean(GP.bestF,1)); % recenter if necessary
            numConds(exp)   = size(groupF{exp},1);
            numCluster(exp) = size(groupF{exp},2);
        
            % Get all the data from all the individuals 
            %load(fullfile(studyDir{exp},'encoding','glm4','cereb_avrgDataStruct.mat'));
            % TT(exp)=
            [X{exp},volIndx,V] = sc1_sc2_functionalAtlas('EVAL:get_data',returnSubjs,exp,'eval');
            GroupX{exp}=mean(X{exp},3); 
        end; 
        

        % Prepare the evaluation by binning the voxel pairs 
        [i,j,k]=ind2sub(V.dim,volIndx);
        [x,y,z]=spmj_affine_transform(i,j,k,V.mat);
        XYZ= [x;y;z];
        % fprintf('calculating distances\n');
        % Dist=single(surfing_eucldist(XYZ,XYZ));
        RR=[];
        clear i j k x y z; % Free memory

        RR=[]; 
        % Loop over Subjects and a)
        for s=1:length(subjs)
            for exp=1:2 
                fprintf('%d %d\n',subjs(s),exp);
                tic; 
                
                % Select data for parcellation 
                ClData = (1-groupW) * X{exp}(:,:,s) + groupW * GroupX{exp};
                numVox=size(ClData,2);
                individG=nan(numCluster(exp),numVox);
                i1=[]; 
                i2=[]; 
                % Now express the X_C= groupF * indivdG
                for i=1:numVox
                    if (mod(i,2000)==0)
                        fprintf('.');
                    end;
                    if sum(abs(ClData(:,i)))>0
                        individG(:,i)=lsqnonneg(double(groupF{exp}),ClData(:,i));
                    end; 
                end
                fprintf(' %d done\n',subjs(s));
                [~,ClusterI]=max(individG,[],1);
                ClusterI(sum(individG,1)==0 | isnan(sum(individG,1)))=nan;

                % Now evaluate this on the data of unique tasks 
                voxIn = find(~isnan(ClusterI)); 
                fprintf('calculating bins\n');
                [BIN,R]=mva_spatialCorrBin(XYZ(:,voxIn),'Parcel',ClusterI(voxIn));  % ,'Dist',Dist(voxIn,voxIn)
                eexp = 3-exp; % Use the other experiment to evaluate 
                load(fullfile(studyDir{eexp},'encoding','glm4','cereb_avrgDataStruct.mat'));
                D1=getrow(D,D.StudyNum~=exp);
                idx=D1.condNum(D1.overlap==0); 
                for c=1:length(idx),
                    i1(c) = find(T.SN==subjs(s) & T.sess==1 & T.cond==idx(c));
                    i2(c) = find(T.SN==subjs(s) & T.sess==2 & T.cond==idx(c));
                end
                Data=(T.data(i1,voxIn)+T.data(i2,voxIn))/2;
            
                fprintf('cross\n');
                R.SN = ones(length(R.N),1)*s;
                R.exp = ones(length(R.N),1)*exp; 
                R.corr = mva_spatialCorr(T.data([i1;i2],voxIn),BIN,...
                    'CrossvalPart',T.sess([i1;i2],1),'excludeNegVoxels',1);
                R.crossval = ones(length(R.corr),1);
                RR = addstruct(RR,R);
%                 fprintf('correl\n');
%                 R.corr=mva_spatialCorr(Data,BIN);
%                 R.crossval = zeros(length(R.corr),1);
%                 RR = addstruct(RR,R);
                toc; 
            end; 
        end        
        varargout={RR}; 
    case 'plotDiffCurves' % group versus indiv for multi-task map
        EvalFiles={'GW0_Sess2_AllC','GW1_Sess2_AllC','GW.3_Sess2_AllC'};
        
        % Aesthetics
        CAT.markertype='none';
        CAT.errorwidth=.5;
        CAT.linecolor={'r','k','b'};
        CAT.errorcolor={'r','k','b'};
        CAT.linewidth=2;
        CAT.linestyle='-';
        D=[]; 
        % Make summary structure 
        for m=1:length(EvalFiles),
            T=load(['Eval_' EvalFiles{m} '.mat']); 
            DD=tapply(T,{'SN','distmin','distmax'},{'corr','mean','subset',T.bwParcel==0,'name','betweenCorr'},...
                                      {'corr','mean','subset',T.bwParcel==1,'name','withinCorr'},...
                'subset',T.crossval==1 & T.distmin<36); 
            DD.type = ones(size(DD.SN,1),1)*m;
            D=addstruct(D,DD); 
        end 
        D.DCBC=D.betweenCorr-D.withinCorr;
       % Labelling
       D.distM = [D.distmin+D.distmax]/2; 
       lineplot(D.distM,D.DCBC,'split',D.type,'leg',EvalFiles,'CAT',CAT); 
        set(gca,'YLim',[0 0.22],'FontSize',12,'XLim',[0 38],'xtick',[0:5:35],'XTickLabel',{'0','','','','','','','35'});
        xlabel('Spatial Distances (mm)');
        ylabel('DCBC');
        set(gcf,'Units','centimeters','PaperPosition',[2,2,6,6]);
        wysiwyg; 
    case 'allCalculation' 
        T=sc1_sc2_GroupIndivid('BuildAndEval','groupW',0);
        save Eval_GW0_Sess2_AllC.mat -struct T
        T=sc1_sc2_GroupIndivid('BuildAndEval','groupW',1);
        save Eval_GW1_Sess2_AllC.mat -struct T
        T=sc1_sc2_GroupIndivid('BuildAndEval','groupW',0.3);
        save Eval_GW.3_Sess2_AllC.mat -struct T         
        T=sc1_sc2_GroupIndivid('BuildAndEval','groupW',0.5);
        save Eval_GW.5_Sess2_AllC.mat -struct T         
        T=sc1_sc2_GroupIndivid('BuildAndEval','groupW',0.1);
        save Eval_GW.1_Sess2_AllC.mat -struct T         
end