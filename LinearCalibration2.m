function [parsEstimate,SigPars,x0,ux0] = LinearCalibration2(x,y,Vx,Vy,Vxy,y0,uy0,options)
%LinearCalibration2 Estimated parameters of the linear calibration function
% and their covariance matrix. Calculation based on the algorithms from ISO 
% TS28037:2010.  
%
% SYNTAX:
% [par,V,x0,ux0] = LinearCalibration2(x,y,Vx,Vy,Vxy,y0,uy0,options)
%
%
% INPUT:
%  x        - vector of x calibration data x = [x1,...,xm]
%  y        - vector of y calibration data y = [y1,...,ym]
%  Vx       - covariance matrix of x measurements, if Vx is a vector it is
%             assumed that instead of the covariance matrix Vx is given the
%             vector ux, i.e. the vector of uncertainties of x
%  Vx       - covariance matrix of y measurements, if Vy is a vector it is
%             assumed that instead of the covariance matrix Vy is given the
%             vector uy, i.e. the vector of uncertainties of y
%  Vxy      - covariance matrix cov(x,y), if Vxy is empty, then it is
%             assumed that Vxy = 0
%  y0       - vector of new indications by the calibrated device
%  uy0      - vector of new indications uncertainties uy0
%  options  - structure with the following parameters:
%             options.isPlot = false % logical indicator for plotting 
%             options.isReversePrediction = true % logical indicator for
%             using the reverse prediction (oposite id the direct
%             prediction method). 
%
% OUTPUT:
%  par      - vector of the estimated parameters of the linear calibration
%             function 
%  V        - covariance matrix of the estimated parameters 
%  x0       - best predictor of the stimulus values based on observed
%             indications y0
%  ux0      - uncertainties of the predicted stimulus values based on observed
%             indications y0 (with given uncertainties uy0)
%
% EXAMPLE 1
% x   = [1.2 1.9 2.9 4.0 4.7 5.9]';
% ux  = [0.2 0.2 0.2 0.2 0.2 0.2]';
% y   = [3.4 4.4 7.2 8.5 10.8 13.5]'; 
% uy  = [0.2 0.2 0.2 0.4 0.4 0.4]';
% y0  = [9 11];
% uy0 = [0.4 0.4];
% [a,Va,x0a,ux0a] = LinearCalibration2(x,y,ux,uy,[],y0,uy0);
%
% EXAMPLE 2
% y0  = [9 11];
% uy0 = [0.4 0.4];
% options.isReversePrediction = false; 
% [b,Vb,x0b,ux0b] = LinearCalibration2(y,x,uy,ux,[],y0,uy0,options);
% 
% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Aug-2021 09:16:39
%% CHECK/SET THE INPUT PARAMETERS
narginchk(2, 8);
if nargin < 8, options = []; end
if nargin < 7, uy0 = []; end
if nargin < 6, y0  = []; end
if nargin < 5, Vxy = []; end
if nargin < 4, Vy  = []; end
if nargin < 3, Vx  = []; end

if ~isfield(options,'isPlot')
    options.isPlot = true;
end

if ~isfield(options,'isReversePrediction')
    options.isReversePrediction = true;
end

if ~isfield(options,'order')
    options.order = 1;
end

if ~isfield(options,'maxiter')
    options.maxiter = 100;
end

if ~isfield(options,'tolerance')
    options.tolerance = 1e-8;
end


if isempty(Vx)
    Vx = ones(size(x));
end


if isempty(Vy)
    Vy = ones(size(y));
end

if isempty(Vxy)
    Vxy = zeros(size(x));
end

if isvector(Vx)
    Vx = diag(Vx.^2);
end

if isvector(Vy)
    Vy = diag(Vy.^2);
end

if isvector(Vxy)
    Vxy = diag(Vxy);
end

U  = [Vx, Vxy; Vxy' Vy];
% ux = sqrt(diag(Vx));
% uy = sqrt(diag(Vy));

%% ALGORITHM
x = x(:);
y = y(:);
m = length(x);
xMin = min(x);
xMax = max(x);
yMin = min(y);
yMax = max(y);

order   = options.order;

SigX    = Vx;
SigY    = Vy;

mu0 = x;
nu0 = y;
pars0 = flip(polyfit(mu0,nu0,order))';

% Initialization
maxiter      = options.maxiter;
tolerance    = options.tolerance;
crit         = 1;
loops        = 1;
muEstimate   = mu0;
nuEstimate   = nu0;
parsEstimate = pars0;

while (crit > tolerance) && (loops <= maxiter)
    mu0         = muEstimate;
    nu0         = nuEstimate;
    pars0       = parsEstimate;
    muOld       = muEstimate;
    nuOld       = nuEstimate;
    parsOld     = parsEstimate;
    [B0,D0]     = update(mu0,nu0,pars0,m,order);
    W0          = D0*SigX*D0 + SigY;
    W0inv       = W0\eye(m);
    Q0          = [W0 B0; B0' zeros(order+1)]\eye(m+order+1);
    Q011        = Q0(1:m,1:m);
    Y           = D0*(x-mu0) - y;
    muEstimate  = x - SigX*D0*Q011*Y;
    nuEstimate  = y + SigY*Q011*Y;
    parsEstimate = -(B0'*W0inv*B0)\(B0'*W0inv*Y);
    crit = norm([muEstimate;nuEstimate;parsEstimate] - ...
        [muOld;nuOld;parsOld])^2/(2*m+order+1);
    loops = loops+1;
end

[B,D,A,c] = update(muEstimate,nuEstimate,parsEstimate,m,order);
W         = D*SigX*D + SigY;
Q         = [W B; B' zeros(order+1)]\eye(m+order+1);
Winv      = W\eye(m);
KY        = (B'*Winv*B)\(B'*Winv);
LX        = -KY * D;
id        = (m+1):(m+order+1);
SigPars   = -Q(id,id);



%% Prediction

a = parsEstimate(1);
b = parsEstimate(2);
u2a = SigPars(1,1);
u2b = SigPars(2,2);
uab = SigPars(1,2);

isReversePrediction = options.isReversePrediction;
if ~isempty(y0)
    if isReversePrediction
        x0  = (y0 - a)/b;
        ux0 = sqrt((u2a + u2b*((y0-a)./b).^2 + 2*uab*(y0-a)./b + uy0.^2)./b^2);
    else
        x0  = a + b*y0;
        ux0 = sqrt(u2a + u2b*y0.^2 + 2*uab*y0 + uy0.^2*b^2);
    end
else
    x0 = [];
    ux0 = [];
end


%% Plot
if options.isPlot
    xx = linspace(floor(min(x)),ceil(max(x)));
    %     set(0, 'DefaultLineLineWidth', 2)
    %     set(0, 'DefaultAxesFontSize', 12)
    figure
    hold on
    grid on
    %errorbar(x,y,uy,uy,ux,ux,'o')
    plot(x,y,'o')
    plot(xx,a+b*xx,'-')
    plot(xx,(a+b*xx)+2*sqrt(u2a+u2b*xx.^2+2*uab*xx),'r--')
    plot(xx,(a+b*xx)-2*sqrt(u2a+u2b*xx.^2+2*uab*xx),'r--')
    if ~isempty(x0)
        if isReversePrediction
            errorbar(x0,y0,2*uy0,2*uy0,2*ux0,2*ux0,'o')
            %errorbar(x0,y0,uy0,uy0,ux0,ux0,'o')
            %errorbar(x0,y0,[],[],2*ux0,2*ux0,'o')
        else
            errorbar(y0,x0,2*ux0,2*ux0,2*uy0,2*uy0,'o')
            %errorbar(y0,x0,ux0,ux0,uy0,uy0,'o')
            %errorbar(y0,x0,2*ux0,2*ux0,[],[],'o')
        end
    end
    axis('square')
    hold off
end

end
%% Function update
function [B,D,A,C] = update(mu0,nu0,a0,n,order,q,Aijk)
% update calculates the update of the matrices A, B and the vector c,
% for a linear calibration model of given order, locally at the
% specified mu0, nu0, and a0. If the coefficient (n,order+1,q)-dimensional
% array Aijk is given, then the updated matrices A, B and the vector c are
% calculated for the generalized linear calibration model.
%
% SYNTAX
%  [B,D,A,C] = update(mu0,nu0,a0,n,order,q,Aijk)
%
% EXAMPLE
%  mu0   = [-19.6843   -9.8142    0.0989   10.0149   19.9634   29.8131]';
%  nu0   = [92.1489   96.0965  100.0499  103.9924  107.9354  111.8316]';
%  a0    = [100.0108    0.3982   -0.0001]';
%  n     = 6;
%  order = 2;
%  q     = [];
%  Aijk  = [];
% [B,D,A,C] = update(mu0,nu0,a0,n,order,q,Aijk)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Apr-2018 13:19:57

%%
if nargin < 7
    Aijk = [];
end

if nargin < 6
    q = [];
end

polyPars = flip(a0);
if isempty(Aijk)
    polyParsD1 = polyder(polyPars);
    D = diag(polyval(polyParsD1,mu0));
    A = [D -eye(n)];
    B = ones(n,order+1);
    for j = 1:order
        B(:,j+1) = mu0.*B(:,j);
    end
    C = polyval(polyPars,mu0) - nu0;
else
    B = zeros(n,order+1);
    D = zeros(n,order+1);
    for i = 1:n
        for j = 1:(order+1)
            for k = 1:q
                B(i,j) = B(i,j) + Aijk(i,j,k)*mu0(i)^(k-1);
                if k > 1
                    D(i,j) = D(i,j) + Aijk(i,j,k)*(k-1)*mu0(i)^(k-2);
                end
            end
        end
    end
    if isempty(a0)
        D = [];
        A = [];
        C = [];
    else
        D = diag(D*a0);
        A = [D -eye(n)];
        C = B*a0 - nu0;
    end
end
end