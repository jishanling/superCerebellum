function varargout=sc1sc2_repSpace(what,varargin)
% Function to experiment with different representational structures 
% from the sc1sc2 data set 
% Directories
baseDir          = '/Users/maedbhking/Documents/Cerebellum_Cognition';
baseDir            = '/Volumes/MotorControl/data/super_cerebellum_new';
% baseDir          = '/Users/jdiedrichsen/Data/super_cerebellum_new';

studyDir{1}     =fullfile(baseDir,'sc1');
studyDir{2}     =fullfile(baseDir,'sc2');
suitDir         ='suit';
caretDir        ='surfaceCaret';
regDir          ='RegionOfInterest';
encodeDir       ='encoding';

subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11',...
    's12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24',...
    's25','s26','s27','s28','s29','s30','s31'};

returnSubjs=[2,3,4,6,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31];

switch what
    case 'TASKSPACE:selectSubset' % Select subset of task that best span the representational space 
        filename = fullfile(studyDir{2},regDir,'glm4','G_hat_sc1_sc2_cerebellum.mat'); 
        load(filename); 
        G=(mean(G_hat,3)); 
        numTasks=10; 
        
        for i=1:100 
            T.indx(i,:) = sample_wor([1:size(G,1)],numTasks,1); 
            G_red = G(T.indx(i,:),T.indx(i,:)); 
            G_red = bsxfun(@minus,G_red,mean(G_red)); 
            G_red = bsxfun(@minus,G_red,mean(G_red)');
            [V,L]=eig(G_red); 
            l = diag(L); 
            T.trL(i,:) = sum(1./l); 
            T.sumL0(i,:) = sum(l>0.00001); 
            T.trL1(i,:) = sum(1./l(l>0.00001)); 
        end;  
            
        keyboard; 
        
    case 'TASKSPACE:overlap'
        
        
end
