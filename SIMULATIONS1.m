%% Determination and use of straight-line calibration function
%  AUTHORS: 
%  Jakub Palenčár, Rudolf Palenčár, Viktor Witkovský and Gejza Wimmer 
%  TITLE: 
%  Linear Comparative Calibration and Uncertainty Analysis of the
%  Measurement Result Obtained With the Calibrated Device 
%  SUBMITTED TO: ENERGIES (MDPI)
%
%  SIMULATIONS
%  Linear calibration of the pressure transducer using a pressure standard
%  according to the ISO Technical Specification ISO/TS 28037:2010 and Monte
%  Carlo method
%  Date: 11-Aug-2021 14:23:41

clear
close all
rng(0)

%%  Set the inputs

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


%% Initial estimate

% NewIndication  = mean(x);
% uNewIndication = sqrt(uxA^2 + uxB^2);
NewIndication = 7.4970;
uxE0 = 0.0038;
uxD0 = 0.015;
uNewIndication = sqrt(uxE0^2 + uxD0^2);

%options.isPlot = false;
options.isReversePrediction  = false;    % Set DIRECT prediction!
[par0,Uab,NewStimulus,uNewStimulus] = ...
    LinearCalibration(x,y,Ux,Uy,Uxy,NewIndication,uNewIndication,options);

ParA0  = par0(1);
%ParA0 =  -15.0167
uParA0 = sqrt(Uab(1,1));
%uParA0 =    0.0133

ParB0  = par0(2);
%ParB0 =   3.7481
uParB0 = sqrt(Uab(2,2));
%uParB0 =  8.4567e-04

%NewStimulus =   13.0829
%uNewStimulus =  0.0588

uParAB0 = Uab(1,2);

alpha = 0.05;
CIa  = [ParA0 + norminv(alpha/2)*uParA0 , ParA0 + norminv(1-alpha/2)*uParA0];
%CIa =  -15.0427  -14.9907

CIb  = [ParB0 + norminv(alpha/2)*uParB0 , ParB0 + norminv(1-alpha/2)*uParB0];
%CIb =    3.7464    3.7498

CInu = [NewStimulus + norminv(alpha/2)*uNewStimulus , NewStimulus + norminv(1-alpha/2)*uNewStimulus];
%CInu =   12.9676   13.1981

%% Simulate the parameters a and b and the stimulus nu0

options.isPlot = false;
xNew  = 7.4970;
uxNew = uNewIndication;
N     = 10000;
par   = zeros(3,N);
for i = 1:N
    xN = x + uxA*randn(n,1) + uxB*sqrt(3)*(rand*2-1);
    yN = y + uyA*randn(n,1) + uyB*randn;
    x0 = xNew + uxE0*randn + uxD0*sqrt(3)*(rand*2-1);
    par(1:2,i) = LinearCalibration(xN,yN,Ux,Uy,[],[],[],options);
    par(3,i) = par(1,i) + par(2,i)*x0;
end

aSort  = sort(par(1,:));
bSort  = sort(par(2,:));
nuSort = sort(par(3,:));
alpha  = 0.05;

%% Histogram / State-of-knowledge distribution of the estimated parameter a

ParA  = mean(aSort);
%ParA =  -15.0167
uParA = std(aSort);
%uParA =  0.0132

CIaMC  = [ParA + norminv(alpha/2)*uParA , ParA + norminv(1-alpha/2)*uParA];
%CIaMC =  -15.0426  -14.9909
CIparA = [aSort(floor(N*(alpha/2))) aSort(ceil(N*(1-alpha/2)))];
%CIparA =  -15.0424  -14.9909

%xA    = linspace(floor(100*min(aSort))/100,ceil(100*max(aSort))/100);
xA    = linspace(-15.08,-14.94)
yA1   = normpdf(xA,ParA,uParA);
yA0   = normpdf(xA,ParA0,uParA0);

figure('DefaultAxesFontSize',12)
histogram(aSort,30,'Normalization','pdf')
hold on
p1 = plot(xA,yA1);
p1.LineWidth = 4;
p1.LineStyle = '-';
p1.Color = [0 0 1];
p2 = plot(xA,yA0);
p2.LineWidth = 2;
p2.LineStyle = '--';
p2.Color = [1 0 0];
xlabel('a')
ylabel('pdf')
title('State-of-knowledge distribution about the coefficient a')
axis([ -15.08     -14.95     0  35])
xticks([-15.08 -15.06 -15.04 -15.02 -15.00 -14.98 -14.96 -14.94])
xticklabels({'-15.08' '-15.06' '-15.04' '-15.02' '-15.00' '-14.98' '-14.96' '-14.94'})


print -dpdf -painters Figure01

%% Histogram / State-of-knowledge distribution of the estimated parameter b

ParB  = mean(bSort);
%ParB =  3.7481
uParB = std(bSort);
%uParB = 8.4896e-04

CIbMC  = [ParB + norminv(alpha/2)*uParB , ParB + norminv(1-alpha/2)*uParB];
%CIbMC =  3.7464    3.7498
CIparB   = [bSort(floor(N*(alpha/2))) bSort(ceil(N*(1-alpha/2)))];
%CIparB =    3.7464    3.7497

%xB    = linspace(floor(1000*min(bSort))/1000,ceil(1000*max(bSort))/1000);
xB    = linspace(3.744,3.752);
yB1   = normpdf(xB,ParB,uParB);
yB0   = normpdf(xB,ParB0,uParB0);

figure('DefaultAxesFontSize',12)
histogram(par(2,:),30,'Normalization','pdf')
hold on
p1 = plot(xB,yB1);
p1.LineWidth = 4;
p1.LineStyle = '-';
p1.Color = [0 0 1];
p2 = plot(xB,yB0);
p2.LineWidth = 2;
p2.LineStyle = '--';
p2.Color = [1 0 0];
xlabel('b')
ylabel('pdf')
title('State-of-knowledge distribution about the coefficient b')
hold off
axis([  3.744    3.752     0  500])
xticks([3.744,3.746,3.748,3.750,3.752])
xticklabels({'3.744','3.746','3.748','3.750','3.752'})

print -dpdf -painters Figure02

%% Histogram / State-of-knowledge distribution of the new stimulus nu_New

nuNew  = mean(nuSort);
%nuNew =   13.0837
unuNew = std(nuSort);
%unuNew =    0.0585

CInuMC  = [nuNew + norminv(alpha/2)*unuNew , nuNew + norminv(1-alpha/2)*unuNew];
%CInuMC =  12.9691   13.1984
CInuNew   = [nuSort(floor(N*(alpha/2))) nuSort(ceil(N*(1-alpha/2)))];

%CInuNew_MC2 =  12.9812  13.1855
%CInuNew_CF  =  12.9810  13.1847
%CInuNew_ISO =  12.9676  13.1981


%xC    = linspace(floor(100*min(nuSort))/100,ceil(100*max(nuSort))/100);
xC    = linspace(12.85,13.3);
yC1   = normpdf(xC,nuNew,unuNew);
yC0   = normpdf(xC,NewStimulus,uNewStimulus);

figure('DefaultAxesFontSize',12)
histogram(nuSort,30,'Normalization','pdf')
hold on
p1 = plot(xC,yC1);
p1.LineWidth = 4;
p1.LineStyle = '-';
p1.Color = [0 0 1];
p2 = plot(xC,yC0);
p2.LineWidth = 2;
p2.LineStyle = '--';
p2.Color = [1 0 0];
xlabel('\nu_0')
ylabel('pdf')
title('State-of-knowledge distribution about the stimulus \nu_0')
hold off
axis([ 12.85   13.3        0   7])
xticks(12.85 :0.05:  13.3 );
xticklabels({'12.85' '12.90' '12.95' '13.00' '13.05' '13.10' '13.15' '13.20' '13.25' '13.30'})


print -dpdf -painters Figure03
% PDF Bouding CUT 93 93 36 38
