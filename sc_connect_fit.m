function [u,R2,R,R2_vox,R_vox,varargout]=sc_connect_fit(Y,X,method,varargin)
% function [u,R2,R,R2_vox,R_vox,varargout]=sc1_encode_fit(Y,X,method,varargin)
%
% INPUT:
%    Y: NxP matrix Data
%    X: NxQ matrix for random effects
%    method: 'linRegress','winnerTakeAll','ridgeFixed'
% VARARGIN:
%    'threshold': Thresholds the u-coefficient at a particular value(s)
%                 before evaluating the prediction (all < threshold -> 0)
% OUTPUT:
%    R2      : correlation value between Y-actual and Y-pred (overall)
%    R       : portion of correctly specified variance (overall)
%    R2_vox  : portion of correctly specified variance (voxels)
%    R_vox   : correlation values between Y-actual and Y-pred (voxels)
%    u       : regression coefficients
% Maedbh King (26/08/2016)
i=0;
[N,P] = size(Y);
[N,Q] = size(X);
lambda=0;
u=zeros(Q,P);
features=[1:Q];

vararginoptions(varargin,{'lambda'});

% Estimate the weights
switch method
    case 'linRegress'              %  Normal linear regression
        u = (X'*X)\(X'*Y);
    case 'nonnegExpSlow'           %  nonnegative regression on log-transform
        u = (X'*X)\(X'*Y);
        u(u<0) = 1e-5;
        theta0  = log(u);
        [theta,fX,iter] = minimize(theta0,@sc1_nonnegExpSlow,100,Y,X);
        u=exp(theta);
    case 'nonNegExp'               %  nonnegative regression on log-transform
        u = (X'*X)\(X'*Y);
        u(u<0) = 1e-5;
        theta0  = log(u);
        XY = X'*Y;   % Precompute for speed
        XX = X'*X;   % Precompute for speed
        
        % checkderiv(@sc1_nonnegExp,theta0,0.0001,XY,XX);
        [theta,fX,iter] = minimize(theta0,@sc1_nonnegExp,1000,XY,XX);
        fX=fX+sum(sum(Y.^2)); % Add data offset;
        u=exp(theta);
    case 'lsqnonneg'               %  matlab internal non-neg least-squares
        [N,P] = size(Y);
        for p=1:P
            u(:,p) = lsqnonneg(X,Y(:,p));
        end;
    case 'cplexnonneg'             %  CPLEX non-neg least-squares
        [N,P] = size(Y);
        if (isnan(lambda) || lambda == 0)
            for p=1:P
                u(:,p) = cplexlsqnonneglin(X,Y(:,p));
            end;
        else
            Aineq=ones(1,Q);
            for p=1:P
                u(:,p) = cplexlsqnonneglin(X,Y(:,p),Aineq,1/lambda);
            end;
        end;
    case {'cplexqp','cplexqpL1'}   %  Non-neg least-squares over quadratic programming - L1
        [N,P]= size(Y);
        [N,Q]= size(X);
        XX=X'*X;
        XY=X'*Y;
        A = -eye(Q);
        b = zeros(Q,1);
        for p=1:P
            u(:,p) = cplexqp(XX,ones(Q,1)*lambda(1)-XY(:,p),A,b);
        end;
    case 'cplexqp_L2'              %  Non-neg least-squares over quadratic programming - L2
        [N,P]= size(Y);
        [N,Q]= size(X);
        XX=X'*X;
        XY=X'*Y;
        A = -eye(Q);
        b = zeros(Q,1);
        for p=1:P
            u(:,p) = cplexqp(XX+lambda(2)*eye(Q),-XY(:,p),A,b);
        end;
    case 'cplexqpL1L2'             %  Non-neg least-squares over quadratic programming - Elastic net
        [N,P]= size(Y);
        [N,Q]= size(X);
        XX=X'*X;
        XY=X'*Y;
        A = -eye(Q);
        b = zeros(Q,1);
        u=nan(Q,P);  % Make non-calculated to nan to keep track of missing voxels
        for p=find(~isnan(sum(Y)))
            u(:,p) = cplexqp(XX+lambda(2)*eye(Q),ones(Q,1)*lambda(1)-XY(:,p),A,b);
        end;
    case 'quadraticProg'           %  Non-neg least-squares over quadratic programming (cplex)
        [N,P] = size(Y);
        [N,Q]= size(X);
        OPT=optimoptions(@quadprog);
        OPT.Display='off';
        XX=X'*X;
        XY=X'*Y;
        A = -eye(Q);
        b = zeros(Q,1);
        for p=1:P
            u(:,p) = quadprog(XX+lambda(2)*eye(Q),ones(Q,1)*lambda(1)-XY(:,p),A,b,[],[],[],[],[],OPT);
        end;
    case 'elasticNet'              %  Matlab's elastic net - determines optimal lambdas 
        [N,P] = size(Y);
        for p=1:P
            u(:,p) = lasso(X,Y(:,p),'Lambda',lambda(1),'Alpha',.5);
        end;
    case 'l1'                      %  Matlab's l1 - determines optimal lambdas
       [N,P] = size(Y);
        for p=1:P
%             u(:,p) = lasso(X,Y(:,p),'Alpha',1); % 1- l1
            [B,S]=lasso(X,Y(:,p)); 
        end;
    case 'l2'                      %  Matlab's l2 - determines optimal lambas
        [N,P] = size(Y);
        for p=1:P
            u(:,p) = lasso(X,Y(:,p),'Alpha',.1); % l2
        end;    
    case 'l1_nonneg'               % l1 ls nonneg -
        [N,P] = size(Y);
        for p=1:P
            [u(:,p)]=l1_ls_nonneg(X,Y(:,p),lambda(1),[],1);
        end;
    case 'winnerTakeAll'
        % get correlation for each network
        yy=sum(Y.*Y,1);
        xx=sum(X.*X,1);
        C=(X'*Y)./sqrt(bsxfun(@times,yy,xx'));
        % get model feature weights for winning network only
        u=zeros(Q,P);
        % limit model feature weights to "winner" network
        for p=1:P,
            [~,I]=max(abs(C(:,p)));
            u(I,p)=X(:,I)'*X(:,I)\X(:,I)'*Y(:,p);
        end
    case 'winnerTakeAll_nonNeg'
        % get correlation for each network
        yy=sum(Y.*Y,1);
        xx=sum(X.*X,1);
        C=(X'*Y)./sqrt(bsxfun(@times,yy,xx'));
        % get model feature weights for winning network only
        u=zeros(Q,P);
        % limit model feature weights to "winner" network
        for p=1:P,
            [~,I]=max(C(:,p));
            u(I,p)=X(:,I)'*X(:,I)\X(:,I)'*Y(:,p);
        end
        u(u<0)=0;
    case 'ridgeFixed'               % L2 regression
        %             u = G*trainX'*((trainX*G*trainX'+eye(sum(trainIdx))*sigma2)\trainY);
        u = (X'*X + eye(Q)*lambda(2))\(X'*Y);
    otherwise
        error ('unknown Method');
end

% Evaluate prediction by calculating R2 and R
SST = nansum(Y.*Y);

Ypred=X*u;
res =Y-Ypred;
SSR = nansum(res.^2);
R2_vox(1,:) = 1-SSR./SST;
R2          = 1-nansum(SSR)/nansum(SST);

% R (per voxel)
SYP = nansum(Y.*Ypred,1);
SPP = nansum(Ypred.*Ypred);

R_vox(1,:) = SYP./sqrt(SST.*SPP);
R          = nansum(SYP)./sqrt(nansum(SST).*nansum(SPP));


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

% Add a L1-norm penality on exp(theta)
function [f,d]=L1_norm(theta,XY,XX,lambda)
f  = -2*sum(sum(XY.*theta))+sum(sum(XX.*(theta*theta')))+lambda*sum(sum(theta));  % Sum of square errors
d  = 2*(-XY + XX *theta).*theta+lambda*theta;

