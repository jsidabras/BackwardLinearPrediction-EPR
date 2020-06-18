% Test the blp on real 3Pulse DEER data collected at the MPI-CEC
% requires:
%   DEERLab v0.8
% 
% Author: Jason W. Sidabras (jason.sidabras@gmail.com)
% Initial writing: 04/05/2020 JWS
% New BLP algorithm: 16/06/2020 JWS
% GPLv3 License.

clear, clf
% Just plot the data to make sure you are looking at the right one.
[traw,Vraw] = deerload('./data/X320060201_JD2073_helix_3pRectDEER_p14_o20-20_overnight.DTA');

plot(traw,Vraw)

%% 3Pulse DEER data without BLP
% Test the blp on real 3Pulse DEER data collected at the MPI-CEC
% requires:
%   DEERLab v0.8
% 
% Author: Jason W. Sidabras (jason.sidabras@gmail.com)
% Initial writing: 04/05/2020 JWS
% New BLP algorithm: 16/06/2020 JWS
% GPLv3 License.
% %%%% Uncomment the 2 outputmatrix commands to save data %%%%%%%

clear, clf
[traw,Vraw] = deerload('./data/X320060201_JD2073_helix_3pRectDEER_p14_o20-20_overnight.DTA');


remove = 15;

Vraw(1:remove) = [];
traw(1:remove) = [];

% Optimization & Correction of phase
V = correctphase(Vraw);

% best guess zero time needed for background 
% ns -> us
t = correctzerotime(V, traw)/1000;
r = time2dist(t);

% Correction of Y-axis scale due to dead-time start
V = V/max(real(V));
% find background and subtract from data
[B,lambda] = fitbackground(V,t,@td_strexp);
Vsub = V - (1 - lambda)*(B);
Vsub = Vsub/max(real(Vsub));

KB = dipolarkernel(t,r);


% regression model fit
Pfit = fitregmodel(Vsub,KB,r,'tikh','aicc');
% gaussian model fit
% maxGauss = 2;
% [Pfit,param,Nopt,metrics,Peval]  = fitmultigauss(Vfull,KB,r,maxGauss...
%     ,'aic', 'Lower',[1.5 0.1]);
Vfit = KB*Pfit;

%Plot results
subplot(211)
plot(t,Vsub,'k.',t,Vfit)
xlabel('t [\mus]')
ylabel('V(t)')
legend('data','fit')
axis tight, grid on, box on
set(gca,'FontSize',14)

% output the signal data
%writematrix([t Vsub' Vfit], '4pMS3psignal_out.csv')

subplot(212)
plot(r,Pfit)
xlabel('r [nm]')
ylabel('P(r)')
legend('distance')
axis tight, grid on, box on
set(gca,'FontSize',14)

% output the solved distance data
% writematrix([r' Pfit], '4pMS3solved_out.csv')

%% Use BLP to fix 3Pulse DEER data 
% Test the blp on real 3Pulse DEER data collected at the MPI-CEC
% requires:
%   DEERLab v0.8
% 
% Author: Jason W. Sidabras (jason.sidabras@gmail.com)
% Initial writing: 04/05/2020 JWS
% New BLP algorithm: 16/06/2020 JWS
% GPLv3 License.
% %%%% Uncomment the 2 outputmatrix commands to save data %%%%%%%

clear, clf
close all
[traw,Vraw] = deerload('./data/X320060201_JD2073_helix_3pRectDEER_p14_o20-20_overnight.DTA');
% parameters
% number of points to cut before analysis
remove = 15;
% number of points to search for.
L = remove+20;
% parameters end
stp = traw(5)-traw(4); % time step


Vraw(1:remove) = [];
traw(1:remove) = [];

% use the blp algorithm to find a back projection
[Vfull] = blp_epr(Vraw',L,100);
Vfull = Vfull';
% Data needs to be transposed to concatenate correctly. 
tnew = traw(1)-stp*L:stp:traw(1)-stp;

% concatenate the projection solution with the experimental data
Vfull = Vfull/max(Vfull);
tfull = [tnew traw'];

% Optimization & Correction of phase
V = correctphase(Vfull);

% best guess zero time needed for background 
% ns -> us
t = correctzerotime(V, tfull)/1000;

% Optimization & Correction of Y-axis scale
V = V/max(real(V));
% find background and subtract from data
[B,lambda] = fitbackground(V,t,@td_strexp);
Vsub = V - (1 - lambda)*(B);
Vsub = Vsub/max(Vsub);

% Optimization & Correction of the found zero-time
% ns -> us
r = time2dist(t);

KB = dipolarkernel(t,r);

% regression model fit
Pfit = fitregmodel(Vsub,KB,r,'tikh','aicc');
% gaussian model fit
% maxGauss = 2;
% [Pfit,param,Nopt,metrics,Peval]  = fitmultigauss(Vfull,KB,r,maxGauss...
%     ,'aic', 'Lower',[1.5 0.1]);
Vfit = KB*Pfit;

%Plot results
subplot(211)
plot(t,Vsub,'k.',t,Vfit)
xlabel('t [\mus]')
ylabel('V(t)')
legend('data','fit')
axis tight, grid on, box on
set(gca,'FontSize',14)

% output the signal data
% writematrix([t' Vfull' Vfit], 'signal_out_n20.csv')

subplot(212)
plot(r,Pfit)
xlabel('r [nm]')
ylabel('P(r)')
legend('distance')
axis tight, grid on, box on
set(gca,'FontSize',14)

% output the solved distance data
% writematrix([r' Pfit], 'solved_out_n20.csv')
