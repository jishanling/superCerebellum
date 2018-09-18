function [R2_vox,R_vox,Ypred,varargout]=encode_crossval(Y,X,partition,method,varargin)
%
% INPUT:
%    Y: NxP matrix Data
%    X: NxQ matrix for random effects
%    partition: Nx1 vector of partitions
%    method:
%    'linRegress'
%    'nonnegExpSlow'
%    'nonnegExp'
%    'nonnegExp_L2'
%    'nonnegExp_L1'
%    'winnerTakeAll'
%    'ridgeFixed'
% OUTPUT:
%    R2         : correlation value between Y-actual and Y-pred (overall)
%    R          : portion of correctly specified variance (overall)
%    R2_vox     : portion of correctly specified variance (voxels)
%    R_vox      : correlation values between Y-actual and Y-pred (voxels)
%
% Maedbh King (26/08/2016)

[N,P] = size(Y);
[N,Q] = size(X);
lambda=0;
<<<<<<< HEAD
y=[];
=======
C=[];
>>>>>>> e63178d691078ec65a55edaa0221b60335df8b4a
features=[1:Q];

vararginoptions(varargin,{'lambda'});
Ypred=zeros(size(Y));

% Loop over partitions
part    = unique(partition)';  % get unique partitions
% part = part(part~=0)';
for i=1:length(part)
    trainIdx = partition~=part(i);
    testIdx  = partition==part(i);
    
    trainY = Y(trainIdx,:);
    trainX = X(trainIdx,:);
    testX  = X(testIdx,:);
    
    % Estimate the weights
    switch method
        case 'linRegress'
            u = (trainX'*trainX)\(trainX'*trainY);
        case 'lsqnonneg'               % matlab internal LSQ non-neg
            [N,P] = size(Y);
            for p=1:P
                u(:,p) = lsqnonneg(trainX,trainY(:,p));
            end;
        case 'cplexqp'                 %  Non-neg least-squares over quadratic programming
            [N,P]= size(Y);
            [N,Q]= size(X);
            XX=trainX'*trainX;
            XY=trainX'*trainY;
            A = -eye(Q);
            b = zeros(Q,1);
            for p=1:P
                u(:,p) = cplexqp(XX,ones(Q,1)*lambda-XY(:,p),A,b);
            end;
        case 'cplexqp_L2'              %  Non-neg least-squares over quadratic programming
            [N,P]= size(Y);
            [N,Q]= size(X);
            XX=trainX'*trainX;
            XY=trainX'*trainY;
            A = -eye(Q);
            b = zeros(Q,1);
            for p=1:P
                u(:,p) = cplexqp(XX+lambda*eye(Q),-XY(:,p),A,b);
            end;
        case 'l1'
            [N,P] = size(Y);
            for p=1:P
                [u(:,p)]=l1_ls(trainX,trainY(:,p),lambda,[],1);
            end;
        case 'lasso'
            [N,P] = size(Y);
            for p=1:P
                u(:,p) = lasso(trainX,trainY(:,p),'Lambda',0.015);
            end;
        case 'l1_nonneg'
            [N,P] = size(Y);
            for p=1:P
                [u(:,p)]=l1_ls_nonneg(trainX,trainY(:,p),lambda,[],1);
            end;
        case 'nonNegExpSlow'           % Slow version of the nonnegative regression
            u = (trainX'*trainX)\(trainX'*trainY);
            u(u<0) = 1e-5;
            theta0  = log(u);
            [theta,fX,iter] = minimize(theta0,@sc1_nonnegExpSlow,1000,trainY,trainX);
            u=exp(theta);
        case 'nonNegExp'               % Fast version - always use this one
            u = (trainX'*trainX)\(trainX'*trainY);
            u(u<0) = 1e-5;
            theta0  = log(u);
            XY = trainX'*trainY;   % Precompute for speed
            XX = trainX'*trainX;   % Precompute for speed
            
            % checkderiv(@sc1_nonnegExp,theta0,0.0001,XY,XX);
            [theta,fX,iter] = minimize(theta0,@sc1_nonnegExp,1000,XY,XX);
            u=exp(theta); % added my mking u=exp(Q);
        case 'nonNegExp_L2'
            u = (trainX'*trainX)\(trainX'*trainY);
            u(u<0) = 1e-5;
            theta0  = log(u);
            XY = trainX'*trainY;   % Precompute for speed
            XX = trainX'*trainX;   % Precompute for speed
            
            G=u*u';
            sigma2=G*(1-.1)/.1; % .1 is pattern consistency value for cerebellum (.1=G/G+sigma2)
            lambda=mean(mean(pinv(G)*sigma2)); % lambda=(G^-1)*sigma2
            % checkderiv(@sc1_nonnegExp_L2,theta0,0.0001,XY,XX,0.1);
            [theta,fX,iter] = minimize(theta0,@sc1_nonnegExp_L2,100,XY,XX,lambda); % mking added lambda
            u=exp(theta);
        case 'nonNegExp_L1'
            u = (trainX'*trainX)\(trainX'*trainY);
            u(u<0) = 1e-5;
            theta0  = log(u);
            XY = trainX'*trainY;   % Precompute for speed
            XX = trainX'*trainX;   % Precompute for speed
            G=u*u';
            sigma2=G*(1-.1)/.1; % .1 is pattern consistency value for cerebellum (.1=G/G+sigma2)
            lambda=mean(mean(pinv(G)*sigma2)); % lambda=(G^-1)*sigma2
            % checkderiv(@sc1_nonnegExp_L1,theta0,0.0001,XY,XX,0.1);
            [theta,fX,iter] = minimize(theta0,@sc1_nonnegExp_L1,100,XY,XX,lambda); %mking added lambda
            u=exp(theta);
        case 'winnerTakeAll'
            % get correlation for each network
            for q=1:Q,
                C(q,:)=corr(trainX(:,q),trainY);
            end
            % get model feature weights
            u=(trainX'*trainX)\trainX'*trainY;
            % limit model feature weights to "winner" network
            for p=1:P,
<<<<<<< HEAD
                [x,y]=max(C,[],1);
                u=zeros(size(u,1),size(u,2)); 
                u(y)=x;
=======
                I=find(C(:,p)==max(C(:,p)));
                other=find(features~=I);
                u(other,p)=1e-5;
>>>>>>> e63178d691078ec65a55edaa0221b60335df8b4a
            end
        case 'winnerTakeAll_nonNeg'
            trainX(trainX<0)=1e-5;
            % get correlation for each network
            for q=1:Q,
                C(q,:)=corr(trainX(:,q),trainY);
            end
            % get model feature weights
            u=(trainX'*trainX)\trainX'*trainY;
            % limit model feature weights to "winner" network
            for p=1:P,
                I=find(C(:,p)==max(C(:,p)));
                other=find(features~=I);
                u(other,p)=1e-5;
            end
        case 'ridgeFixed'
            %           u = G*trainX'*((trainX*G*trainX'+eye(sum(trainIdx))*sigma2)\trainY);
            u = (trainX'*trainX + eye(Q)*lambda)\(trainX'*trainY);
        otherwise
            error ('unknown Method');
    end
    
    % Make prediction
    Ypred(testIdx,:)=testX*u;
end;


% Evaluate prediction by calculating R2 and R
% -------------------------------------------------
SST = nansum(Y.*Y);
res=Y-Ypred;
SSR = nansum(res.^2);
R2_vox(1,:) = 1-SSR./SST;
R2          = 1-nansum(SSR)/nansum(SST);

% R (per voxel)
SYP = nansum(Y.*Ypred,1);
SPP = nansum(Ypred.*Ypred);

R_vox(1,:) = SYP./sqrt(SST.*SPP);
R          = nansum(SYP)./sqrt(nansum(SST).*nansum(SPP));

<<<<<<< HEAD
varargout={R2,R,y};
=======
varargout={R2,R,C};
>>>>>>> e63178d691078ec65a55edaa0221b60335df8b4a

% Derivative functions for models
% Basic non-negative regression without a prior
% This is the explicit, slow version
function [f,d]=sc1_nonnegExpSlow(theta,Y,X)
u=exp(theta);
res = Y - X*u;
f  = sum(sum(res.*res));        % Sum of square errors
d  = (-2*X'*Y + 2*X'*X *u).*u;  % Derivative of f in respect to theta

% This is the corresponding fast version of the optimisation - here the
% products of XY and XX are precomputed outside the function. We drop the
% term trace(Y*Y') from the squared error, as it does not depend on the
% parameters.
function [f,d]=sc1_nonnegExp(theta,XY,XX)
u=exp(theta);
f  = -2*sum(sum(XY.*u))+sum(sum(XX.*(u*u')));  % Sum of square errors
d  = 2*(-XY + XX *u).*u;

% Now add a L2-norm penality on exp(theta)
function [f,d]=sc1_nonnegExp_L2(theta,XY,XX,lambda)
u=exp(theta);
u2 = u.*u;
f  = -2*sum(sum(XY.*u))+sum(sum(XX.*(u*u')))+lambda*sum(sum(u2));
d  = 2*(-XY + XX *u).*u+lambda*2*u2;

% Add a L1-norm penality on exp(theta)
function [f,d]=sc1_nonnegExp_L1(theta,XY,XX,lambda)
u=exp(theta);
f  = -2*sum(sum(XY.*u))+sum(sum(XX.*(u*u')))+lambda*sum(sum(u));  % Sum of square errors
d  = 2*(-XY + XX *u).*u+lambda*u;

