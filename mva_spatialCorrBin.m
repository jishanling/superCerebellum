function [BIN,T]=mva_spatialCorrBin(XYZ,varargin)
% function [BIN,T]=mva_spatialCorrBin(XYZ,Data,varargin)
% INPUT:
%   XYZ: Coordinates of voxels, a 3xP matrix (in mm)
%  VARARGIN:  
%   Parcel: A 1xP parcellation vector - by default compartment 0 will be
%         ignored - then computes the correlation for between and within
%         parcel voxel pairs 
%   spatialBins: K+1 boundaries between the different bins of distances 
%       between voxels / vertices 
%  OUTPUT: 
%       BIN: large PxP matrix of bin indices 
%       T: Data frame with the fields 
%           T.dist = mean distance between voxels for this bin 
%           T.N    = number of voxel pairs 
Parcel = []; 
spatialBins = [0:5:60 inf]; 
vararginoptions(varargin,{'Parcel','spatialBins'}); 
P = size(XYZ,2); 
D=surfing_eucldist(XYZ,XYZ);        % Eucledian distance 
BIN = zeros(P,P,'uint8');  % Category label 

% First drop all voxels pairs into a different category 
numBins = length(spatialBins)-1; 
for i=1:numBins
    ind = D>spatialBins(i) & D<=spatialBins(i+1); 
    BIN(ind)=i;
end; 

% split by within and between 
if (~isempty(Parcel))
    A=bsxfun(@ne,Parcel',Parcel) & BIN>0; 
    BIN=BIN+uint8(A)*uint8(numBins);
    T.bwParcel = [zeros(numBins,1);ones(numBins,1)]; 
    T.bin     = kron([1;1],[1:numBins]'); 
    T.distmin = kron([1;1],spatialBins(1:numBins)'); 
    T.distmax = kron([1;1],spatialBins(2:numBins+1)'); 
    numBins=numBins*2; 
else
    T.bin     = [1:numBins]'; 
    T.distmin = spatialBins(1:numBins); 
    T.distmax = spatialBins(2:numBins+1); 
end;

% Determine number and mean difference 
for i=1:numBins
    ind = BIN==i; 
    T.N(i,1) = sum(ind(:)); 
    T.dist(i,1) = sum(D(ind))/T.N(i); 
end; 
