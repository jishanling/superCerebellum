function [A_PW,S_PW,W_PW,winner]=pca_ica(X_C,varargin)
% function [A_PW,S_PW,W_PW,winner]=pca_ica(X_C,varargin)
% INPUT:
%   Data: RxP matrix of voxel / vertex data to do pca + ica
%   Matrix should be centered before being submitted
%  VARARGIN:  
%       threshold: percent of variance to keep. Default is .9
%  OUTPUT: 
%       A_PW: R by P mixing matrix. R - number of regressors 
%       S_PW: K by P separated signal matrix. K - number of indep. comps
%       W_PW: K by R separating matrix. 
%       winner: 1 x P. Best (i.e. max) PC per voxel. 
%  Maedbh King. March 2018

threshold=.9; % default threshold

vararginoptions(varargin,{'threshold'}); 

% calculate the pca for data
[E,D] = pcamat(X_C,1,size(X_C,1),'off','off');

% get PCs explaining <var> of POV
[EV,I]=sort(diag(D),'descend');
for i=1:length(EV),
    POV(i)=EV(i)/trace(D);
end
EVindx=find(cumsum(POV)<=threshold);

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

numFeat=size(S_PW,1);

% do winner-take-all
for p=1:size(S_PW,2),
    signIdx=find(S_PW(:,p)>0);
    pwS=abs(S_PW(:,p));
    [~,i]=sort(pwS,'descend');
    % sort 'winner' into pos or neg group
    if find(signIdx==i(1)),
        winner(p,:)=i(1);
    else
        winner(p,:)=i(1)+numFeat;
    end
end

% rename data points
featNum=unique(winner);
newFeat=[1:length(featNum)];
for ff=1:length(featNum),
    winner(winner==featNum(ff))=newFeat(ff);
end
