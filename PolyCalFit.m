%% Determination and use of straight-line calibration function
%  
%  For more details see
%  Palenčár J., Palenčár R., Chytil M., Witkovský V., Wimmer G. Jr., and
%  Wimmer G.: On Linear Comparative Calibration and Uncertainty Analysis of
%  the Measurement Result ObtainedWith the Calibrated Instrument, Working
%  Paper 2022.   
%
%  EVALUATION EXAMPLE
%  Linear calibration of the pressure transducer using a pressure standard
%  according to the ISO Technical Specification ISO/TS 28037:2010 and 
%  Characteristic Functions Approach
%  Date: 20-Aug-2021 09:37:23

clear
close all

%% MEASUREMENT DATA 
%  The best estimates observed in the calibration experiment of the TSZ
%  pressure transducer realized at the Faculty of Mechanical Engineering of
%  the Slovak University of Technology in Bratislava.

DATA = [...
1   0.0000  4.0030;
2   10.191  6.7160;
3   20.102  9.3710;
4   30.170  12.053;
5   42.230  15.266;
6   50.050  17.351;
7   60.070  20.036;
8   50.080  17.369;
9   40.115  14.718;
10  30.089  12.039;
11  20.095  9.3760;
12  10.070  6.6970;
13  0.0000  4.0080];

%% x => Calibrated Device; y => Standard Device
n   = 13;
x   = DATA(:,3);            % mA
uxA = sqrt(0.00001444);     % mA
uxB = 0.0014;               % mA
Ux  = uxA^2*eye(n) + uxB^2*ones(n);


y   = DATA(:,2);            % kPa
uyA = sqrt(0.000036);       % kPa
uyB = sqrt(0.000036);       % kPa
Uy  = uyA^2*eye(n) + uyB^2*ones(n);

Uxy = [];



%% Linear Calibration according to ISO/TS 28037:2010
NewIndication = 7.4970;
uxE0 = 0.0038;
uxD0 = 0.0150;
uNewIndication = sqrt(uxE0^2 + uxD0^2);
%NewIndication  = mean(x);
%uNewIndication = sqrt(uxA^2 + uxB^2);
options.isReversePrediction  = false;    % Set DIRECT prediction!

[ab,Uab,NewStimulus,uNewStimulus] = ...
    LinearCalibration(x,y,Ux,Uy,Uxy,NewIndication,uNewIndication,options);


%% Control calculations

x0  = NewIndication;
%x0 =   7.4970

ux0 = uNewIndication;
%ux0 = 0.0155

a   = ab(1);
%a = -15.0167

b   = ab(2);
%b =    3.7481 

ua  = sqrt(Uab(1,1));
%ua =    0.0133

ub  = sqrt(Uab(2,2));
%ub =   8.4567e-04

uab = Uab(1,2);
%uab =  -8.1970e-06

nu0  = a + b * x0;
%nu0 =   13.0829

unu0 = sqrt(ua^2 + x0^2*ub^2 + 2*x0*uab + b^2*ux0^2);
%unu0 =  0.0588

alpha = 0.05;
CInu0 = [NewStimulus + norminv(alpha/2)*uNewStimulus , NewStimulus + norminv(1-alpha/2)*uNewStimulus];
%CInu0 =   12.9676   13.1981

%% CF New indication / set the xFit values and the weights 

NewIndication = 7.4970;
uxE0 = 0.0038;
uxD0 = 0.0150;
cfx0 = @(t)cf_RectangularSymmetric(sqrt(3)* uxD0*t) .* cf_Normal(uxE0*t) .* exp(1i*t*NewIndication);
nN   = 51;
%prob = linspace(1e-10,1-1e-10,nN);
[prob,wNew] = ChebPoints(nN,[0,1]);
resultx0 = cf2DistGP(cfx0,[],prob);
xNew = resultx0.qf;

%% PolyCal FIT

clear options
options.order = 1;
options.uX    = uxA*ones(size(x));
options.uX0   = uxB;
options.uY    = uyA*ones(size(x));
options.uY0   = uyB;
options.cfX   = [];
options.cfX0  = @(t) cf_RectangularSymmetric(sqrt(3)*uxB*t);
options.cfY   = [];
options.cfY0  = [];

res = PolyCal(x,y,xNew,[],[],options)


%% Plot a

cfa = @(t)res.cfPars{1}(t) .* exp(1i*t*res.pars(1));
resulta = cf2DistGP(cfa)

figure
plot(resulta.x,resulta.pdf,'g','LineWidth',2)

%% Plot b

cfb = @(t)res.cfPars{2}(t) .* exp(1i*t*res.pars(2))
resultb = cf2DistGP(cfb)

figure
plot(resultb.x,resultb.pdf,'g','LineWidth',2)

%% Plot yNew

cfyNew = @(t) wNew(1) *  res.cfFit{1}(t);
for i = 2:nN
    cfyNew = @(t) (cfyNew(t)  + wNew(i) * res.cfFit{i}(t));
end
%cfyNew = @(t)cfyNew(t) / nN;


prob = [0.025,0.975];
resultyNew = cf2DistGP(cfyNew,[],prob)

figure
plot(resultyNew.x,resultyNew.pdf,'g','LineWidth',2)

%% Plot CF

cfyNew0 = @(t) cfyNew(t) .* exp(- 1i*t*resultyNew.xMean)

t = linspace(-200,200,1001);
figure
plot(t,real(cfyNew0(t)),t,imag(cfyNew0(t)))

%% PDF from cfyNew0

resultyNew0 = cf2DistGP(cfyNew0,[],prob)
