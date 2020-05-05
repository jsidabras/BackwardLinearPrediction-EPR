% Test the blp on real 3Pulse DEER data collected at the MPI-CEC
% requires:
%   DEERLab v0.8
% 
% Author: Jason W. Sidabras (jason.sidabras@gmail.com)
% Initial writing: 04/05/2020 JWS
% GPLv3 License.

clear, clf
[traw,Vraw] = deerload('./data/X320030405_JD2060.DTA');

n = 15;
q = 10*n;
M = 20;

Vraw(1:6) = [];
Vraw(200:end) = [];
traw(1:6) = [];
traw(200:end) = [];
% Optimization & Correction of phase
V = correctphase(Vraw);

% best guess zero time needed for background 
t = correctzerotime(V, traw);

% Optimization & Correction of Y-axis scale
V = V/max(real(V));
% find background and subtract from data
[B,lambda] = fitbackground(V,t,@td_strexp);
Vsub = V - (1 - lambda)*(B);

% use the blp algorithm to find a back projection
backpred = blp_epr(Vsub, n, q, M);
tnew = traw(1)-4*M:4:traw(1)-4;

% concatenate the projection solution with the experimental data
Vfull = [backpred Vsub];
Vfull = Vfull/max(Vfull);
tfull = [tnew traw'];

% Optimization & Correction of the found zero-time
% ns -> us
t = correctzerotime(Vfull,tfull)/1000;
r = time2dist(t);

KB = dipolarkernel(t,r);


% Pfit = fitregmodel(Vfull,KB,r,'tikh','aicc');
maxGauss = 2;
[Pfit,param,Nopt,metrics,Peval]  = fitmultigauss(Vfull,KB,r,maxGauss,'aic');
Vfit = KB*Pfit;

%Plot results
subplot(211)
plot(t,Vfull,'k.',t,Vfit)
xlabel('t [\mus]')
ylabel('V(t)')
legend('data','fit')
axis tight, grid on, box on
set(gca,'FontSize',14)

subplot(212)
plot(r,Pfit)
xlabel('r [nm]')
ylabel('P(r)')
legend('distance')
axis tight, grid on, box on
set(gca,'FontSize',14)

% maxF = max(Pfit);  % Find max value over all elements.
% indexOfFirstMax = find(Pfit == maxF, 1, 'first');  % Get first element that is the max.
% % Get the x and y values at that index.
% maxX = r(indexOfFirstMax)
param