function [ab,Uab,x0,ux0] = LinearCalibration(x,y,Ux,Uy,Uxy,y0,uy0,options)
%LinearCalibration Estimates the parameters of the linear calibration
%  function and their covariance matrix. Calculation is based on using the
%  algorithms from ISO TS28037:2010. 
%
%  The original algorithms have been implemented at the National Physical
%  Laboratory (NPL), the UK's National Metrology Institute, as the NPL's
%  software to support ISO/TS 28037:2010(E), see [1], and are freely
%  available at https://www.npl.co.uk/resources/software/iso-ts-28037-2010e.    
%
% SYNTAX:
% [ab,Uab,x0,ux0] = LinearCalibration(x,y,Ux,Uy,Uxy,y0,uy0,options)
%
% INPUTS:
%  x        - vector of best estimates x - calibration data x = [x1,...,xm]
%  y        - vector of best estimates y - calibration data y = [y1,...,ym]
%  Ux       - covariance matrix of x measurements, if Ux is a vector it is
%             assumed that instead of the covariance matrix Ux is given the
%             vector ux, i.e. the vector of uncertainties of x
%  Ux       - covariance matrix of y measurements, if Uy is a vector it is
%             assumed that instead of the covariance matrix Uy is given the
%             vector uy, i.e. the vector of uncertainties of y
%  Uxy      - covariance matrix cov(x,y), if Uxy is empty, then it is
%             assumed that Uxy = 0
%  y0       - vector of new indications by the calibrated device
%  uy0      - vector of new indications uncertainties uy0
%  options  - structure with the following parameters:
%             options.isPlot = false % logical indicator for plotting 
%             options.tol = 1e-8     % tolerance for the stopping rule
%             options.isReversePrediction = true % logical indicator for
%             using the reverse prediction (oposite is the direct
%             prediction method). 
%
% OUTPUTS:
%  ab       - vector of the estimated parameters of the linear calibration
%             function 
%  Uab      - covariance matrix of the estimated parameters 
%  x0       - best predictor of the stimulus values based on observed
%             indications y0
%  ux0      - uncertainties of the predicted stimulus values based on
%             observed indications y0 (with given uncertainties uy0)
%
% EXAMPLE 1
%  x   = [1.2 1.9 2.9 4.0 4.7 5.9]';
%  ux  = [0.2 0.2 0.2 0.2 0.2 0.2]';
%  y   = [3.4 4.4 7.2 8.5 10.8 13.5]'; 
%  uy  = [0.2 0.2 0.2 0.4 0.4 0.4]';
%  y0  = [9 11];
%  uy0 = [0.4 0.4];
%  [ab,Uab,x0a,ux0a] = LinearCalibration(x,y,ux,uy,[],y0,uy0);
%
% EXAMPLE 2
%  y0  = [9 11];
%  uy0 = [0.4 0.4];
%  options.isReversePrediction = false; 
%  [b,Vb,x0b,ux0b] = LinearCalibration(y,x,uy,ux,[],y0,uy0,options);
% 
% REFERENCES
% [1] Smith I., Harris, P. NPL's software to support ISO/TS 28037:2010(E).
%     Available: https://www.npl.co.uk/resources/software/iso-ts-28037-2010e. 
% [2] Palencar, J.; Palencar, R.; Witkovsky, V.; Wimmer, G. Linear
%     Comparative Calibration and Uncertainty Analysis of the Measurement
%     Result Obtained With the Calibrated Device. Version August 4, 2021.
%     Submitted to Energies. https://www.mdpi.com/journal/energies. 
% 
% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 11-Aug-2021 14:23:41

%% CHECK/SET THE INPUT PARAMETERS
narginchk(2, 8);
if nargin < 8, options = []; end
if nargin < 7, uy0 = []; end
if nargin < 6, y0  = []; end
if nargin < 5, Uxy = []; end
if nargin < 4, Uy  = []; end
if nargin < 3, Ux  = []; end

if ~isfield(options,'isPlot')
    options.isPlot = true;
end

if ~isfield(options,'tol')
    options.tol = 1e-8;
end

if ~isfield(options,'isReversePrediction')
    options.isReversePrediction = true;
end

if isempty(Ux)
    Ux = ones(size(x));
end


if isempty(Uy)
    Uy = ones(size(y));
end

if isempty(Uxy)
    Uxy = zeros(size(x));
end

if isvector(Ux)
    Ux = diag(Ux.^2);
end

if isvector(Uy)
    Uy = diag(Uy.^2);
end

if isvector(Uxy)
    Uxy = diag(Uxy);
end

Uab  = [Ux, Uxy; Uxy' Uy];

%% ALGORITHM
x = x(:);
y = y(:);

m = length(x);
[ai, bi] = algm_wls_steps_1_to_8(x, y, sqrt(diag(Uy)));
tt{1} = [x; ai; bi]; 
tol   = options.tol;
dt{1} = []; 
f{1}  = []; 
J{1}  = []; 
ft{1} = []; 
Jt{1} = []; 
g{1}  = []; 
H{1}  = []; 
M{1}  = []; 
q{1}  = []; 

ind   = 1;
[tt, dt, f, J, ft, Jt, g, H, M, q] ... 
  = algm_ggmr_cholesky_steps_2_to_9(x,y,Uab,tt,dt,f,J,ft,Jt,g,H,M,q,ind); 

while any(abs(dt{ind}) > tol)
    ind = ind + 1;
    [tt, dt, f, J, ft, Jt, g, H, M, q] ...
        = algm_ggmr_cholesky_steps_2_to_9(x,y,Uab,tt,dt,f,J,ft,Jt,g,H,M,q,ind);
end

a = tt{ind+1}(m+1);
b = tt{ind+1}(m+2);

M22 = M{ind}(m+1:m+2, m+1:m+2); 
m11 = M22(1, 1); 
m21 = M22(2, 1); 
m22 = M22(2, 2); 
u2a = (m22^2 + m21^2)/(m11^2*m22^2); 
u2b = 1/(m22^2); 
uab = -m21/(m11*m22^2); 

ab = [a;b];
Uab   = [u2a,uab;uab,u2b];

%% Prediction
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
    figure
    hold on
    grid on
    plot(x,y,'o')
    plot(xx,a+b*xx,'-')
    plot(xx,(a+b*xx)+2*sqrt(u2a+u2b*xx.^2+2*uab*xx),'r--')
    plot(xx,(a+b*xx)-2*sqrt(u2a+u2b*xx.^2+2*uab*xx),'r--')
    if ~isempty(x0)
        if isReversePrediction
            errorbar(x0,y0,2*uy0,2*uy0,2*ux0,2*ux0,'o')
        else
            errorbar(y0,x0,2*ux0,2*ux0,2*uy0,2*uy0,'o')
        end
    end
    axis('square')
    hold off
end

%% Function to solve the weighted least squares (WLS) problem (Clause 6)
function [a, b, u2a, u2b, uab, w, g0, h0, g, h, G2, r, R] = ...
    algm_wls_steps_1_to_8(x, y, uy)
%  The original algorithms have been implemented at the National Physical
%  Laboratory (NPL), the UK's National Metrology Institute, as the NPL's
%  software to support ISO/TS 28037:2010(E), see [1], and are freely
%  available at https://www.npl.co.uk/resources/software/iso-ts-28037-2010e.

% Step 1. 
m = length(x);
w = ones(m, 1)./uy;
F2 = sum(w.*w);

% Step 2. 
g0 = (sum(w.*w.*x))/F2; 
h0 = (sum(w.*w.*y))/F2; 

% Step 3. 
g = w.*(x - g0); 
h = w.*(y - h0); 

% Step 4. 
G2 = sum(g.*g); 

% Step 5. 
b = (sum(g.*h))/G2; 
a = h0 - b*g0; 

% Step 6. 
u2a = 1/F2 + (g0^2)/G2; 
u2b = 1/G2; 
uab = -g0/G2; 

% Step 7. 
r = w.*(y - a - b*x); 

% Step 8. 
R = sum(r.*r); 

% End of algm_wls_steps_1_to_8.m
end
%% Function to undertake steps 2 to 9 of the GGMR algorithm (Clause 10)
function [tt, dt, f, J, ft, Jt, g, H, M, q] = ... 
  algm_ggmr_cholesky_steps_2_to_9(x, y, U, tt, dt, f, J, ft, Jt, g, H, M, q, ind)
%  The original algorithms have been implemented at the National Physical
%  Laboratory (NPL), the UK's National Metrology Institute, as the NPL's
%  software to support ISO/TS 28037:2010(E), see [1], and are freely
%  available at https://www.npl.co.uk/resources/software/iso-ts-28037-2010e.

% Step 2. 
m = length(x); 
f{ind} = [x - tt{ind}(1:m); 
          y - (tt{ind}(m+1) + tt{ind}(m+2)*tt{ind}(1:m))]; 
J{ind} = [-eye(m), zeros(m, 2); 
          -tt{ind}(m+2)*eye(m), -ones(m, 1), -tt{ind}(1:m)]; 

% Step 3. 
L = chol(U, 'lower'); 

% Step 4. 
ft{ind} = L\f{ind}; 
Jt{ind} = L\J{ind}; 

% Step 5. 
g{ind} = Jt{ind}'*ft{ind}; 
H{ind} = Jt{ind}'*Jt{ind}; 

% Step 6. 
M{ind} = chol(H{ind}, 'lower'); 

% Step 7. 
q{ind} = -M{ind}\g{ind}; 

% Step 8. 
dt{ind} = M{ind}'\q{ind}; 

% Step 9. 
tt{ind+1} = tt{ind} + dt{ind}; 
% End of algm_ggmr_cholesky_steps_2_to_9.m
end
end