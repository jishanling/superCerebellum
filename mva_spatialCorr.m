function [corr,T,varargout]=mva_spatialCorr(Data,BIN,varargin)
% function T=mva_spatialCorr(Data,varargin)
% INPUT:
%   Data: KxP matrix of voxel / vertex data to compute spatial correlation
%           over
%          Note the function computes corr = x'*y/sqrt(x'*x * y'*y); 
%          So Data needs to be mean-subtracted if you want to get a
%          correlation, other
%   BINS: A PxP matrix dividing the any voxel / vertex pair into a category. 
%  VARARGIN:  
%       between voxels / vertices 
%  OUTPUT: 
%       Corr: numBins-vector of spatial correlations 
%           T: Data frame with the fields 
%           T.corr= mean distance between voxels for this bin 
%           T.meanCOV = mean covariance 
%           T.meanVAR = mean variance 

% First compute the covariances between tuning function
% Maedbh King edited this function to compute voxel-wise correlations
% (as an addition to the bin-wise corr)
CrossvalPart =[]; 
excludeNegVoxels = 1; 
saveVox = 'no'; 
numBins =[]; 
verbose =0; 
vararginoptions(varargin,{'CrossvalPart','excludeNegVoxels','numBins','verbose','saveVox'}); 


if (isempty(numBins))
    numBins = max(BIN(:)); 
end; 

if (verbose) 
    fprintf('cov'); 
end; 

if (~isempty(CrossvalPart)) 
    numPart = max(CrossvalPart); 
    if numPart>2 
        error('currently only split-half partitions - numPart = 2'); 
    end; 
    i1 = CrossvalPart == 1; 
    i2 = CrossvalPart == 2; 
    K  = sum(i1); 
    K2 = sum(i2); 
    if K~=K2 
        error('partitions need to be exactly the same size'); 
    end; 
    VAR = sum(Data(i1,:).*Data(i2,:))/K;
    
    % Now select only voxels with reliable tuning 
    if (excludeNegVoxels) 
        subindx = VAR>0; 
        VAR = VAR(subindx); 
        Data = Data(:,subindx); 
        BIN  = BIN(subindx,subindx); 
    end; 
else
    K = size(Data,1);
    VAR = sum(Data.^2)/K;
end; 

VAR = (VAR'*VAR);
for i=1:numBins
    T.meanVAR(i,1)=nanmean(VAR(BIN==i));
end; 
%clear VAR;

% And finally Covariance 
if (verbose) 
    fprintf('cov'); 
end; 
if (~isempty(CrossvalPart)) 
    COV = Data(i1,:)'*Data(i2,:)/K; 
else
    COV = Data'*Data/K; 
end;
for i=1:numBins
    T.meanCOV(i,1)=nanmean(COV(BIN==i)); 
end; 
T.corr=T.meanCOV./sqrt(T.meanVAR); 
corr = T.corr; 

% save out voxel-wise corr (binned)
if strcmp('saveVox','yes'),
    % voxel-wise correlations
    corrVox = COV./sqrt(VAR);
    % binned voxel-wise correlations (probably a faster way of doing this ..)
    for v=1:size(corrVox,1),
        for i=1:numBins,
            binCorrVox(v,i)=nanmean(corrVox(v,BIN(v,:)==i));
            fprintf('bin%d: voxel %d \n',i,v)
        end
    end
end
varargout={'binCorrVox'}; 
