% Test the blp on real 3Pulse DEER data collected at the MPI-CEC
% requires:
%   DEERLab v0.8
% 
% Author: Jason W. Sidabras (jason.sidabras@gmail.com)
% Initial writing: 04/05/2020 JWS
% GPLv3 License.

clear, clf

% parameters
n_std = 0.025; % noise standard deviation 0.025 ~ SNR 50
rmean = 2.; %nm
width = 0.15; %nm
t_cut = [4 6 8 10 12 15 20 30]; % series of points to cut
nstat = 500; % number of independant trials
n = 150;
q = 25;
M = t_cut+20;
% parameters end

tmin = 0; %us
tmax = 1; %us
% 2 ns steps
N = 500;
%Time-domain axis
t = linspace(tmin,tmax,N); %us
stp = (t(2)-t(1));

rmin = 1; %nm
rmax = 4; %nm
N = 500;
%Time-domain axis
r = linspace(rmin,rmax,N); %nm

%Generate a distance distribution
P = rd_onegaussian(r,[rmean width]);
Sfree = dipolarsignal(t,r,P);

rmse_ans = [];
for icut = 1:length(t_cut)
    list_rmean = zeros(1,nstat);
    list_sigma = zeros(1,nstat);
    for i = 1:nstat
        %Simulate dipolar evolution function
        Snoise = dipolarsignal(t,r,P,'noiselevel',n_std);
        Snoise = Snoise';
        Snoisecut = Snoise(t_cut(icut)+1:end);
        timecut = t(t_cut(icut)+1:end);

        
        % use the blp algorithm to find a back projection
        backpred = blp_epr(Snoisecut,n,q,M(icut));
        tnew = timecut(1)-stp*M(icut):stp:timecut(1)-stp;

        % concatenate the projection solution with the experimental data
        Vfull = [backpred Snoisecut];
        Vfull = Vfull/max(Vfull);
        tfull = [tnew timecut];
        % best guess zero time needed for background 
        % ns -> us
        tfix = correctzerotime(Vfull, tfull);

        KB = dipolarkernel(tfix,r);

        % regression model fit
        % Pfit = fitregmodel(Vfull,KB,r,'tikh','aicc');
        % gaussian model fit
        maxGauss = 2;
        [Pfit,param,Nopt,metrics,Peval]  = fitmultigauss(Vfull,KB,r,maxGauss,'aic','Lower',[1.5 0.1]);
        Vfit = KB*Pfit;

        list_rmean(i) = param(1);
        list_sigma(i) = param(2);
        
    end
%Plot results
subplot(211)
plot(tfull,Vfull,'k.',tfull,Vfit)
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
    fprintf('Points cut: %i \n', t_cut(icut))
    fprintf('r_mean = %.4f +/- %.4f\n', mean(list_rmean), std(list_rmean))
    fprintf('sigma_mean = %.4f +/- %.4f\n', mean(list_sigma), std(list_sigma))
end

%Plot
% figure('position',[0 0 500 200])
% plot(t,Snoise,'k',t,Sfree,'r','Linewidth',1.5)
% set(gca,'fontsize',14)
% axis tight, grid on
% xlabel('t [\mus]'),ylabel('S(t)')