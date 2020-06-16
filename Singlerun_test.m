% Test the blp on a single randomly generated noisy dipolar signal
% requires:
%   DEERLab v0.8
% 
% Author: Jason W. Sidabras (jason.sidabras@gmail.com)
% Initial writing: 04/05/2020 JWS
% New BLP algorithm: 16/06/2020 JWS
% GPLv3 License
rmean = 2.; %nm
width = 0.25; %nm
tmin = 0; %us
tmax = 1; %us
% 2 ns steps
N = 500;
% number of points to predict. In typical use, you would set the 
% to what you think zero should be. In this test script, this automatically 
% removes these points from the simulated data so we already know 
% where zero will be. 
t_cut = 10;
% Extra points to predict before zero. Typically, setting this to 10 is
% good idea to get enough points to find the point that the signal comes
% back down.
extra = 10; 

% To generate an approximate SNR of 50, set this to 0.025
% the range is typically an SRN of 46-55
noisenum = 0.025;
% parameters end
L = t_cut+extra;


%Time-domain axis
t = linspace(tmin,tmax,N); %us
stp = (t(2)-t(1));

rmin = 1; %nm
rmax = 4; %nm
N = 500;
%Time-domain axis
r = linspace(rmin,rmax,N); %nm

% Generate a distance distribution
P = rd_onegaussian(r,[rmean width]);

% generate noiseless 
tfree = linspace(tmin-extra*0.002,tmax,N+extra); %us
Sfree = dipolarsignal(tfree,r,P);

% generate noisy signal
y = dipolarsignal(t,r,P,'noiselevel',noisenum );
y = y';
y = y(t_cut+1:end);

% use the BLP algorithm to find the missing points!
[g, gblp] = blp_epr(y,L,25);
figure(1)
clf
plot(g)
hold on
plot(gblp,'r')
title('Data and predicted values')


SNR = (max(g)-min(g))/std(g(400:end));
fprintf('SNR value: %f1000.0 \n', SNR)
