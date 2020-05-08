% Test the blp on real 3Pulse DEER data collected at the MPI-CEC
% requires:
%   DEERLab v0.8
% 
% Author: Jason W. Sidabras (jason.sidabras@gmail.com)
% Initial writing: 04/05/2020 JWS
% GPLv3 License
clear, clf
clear all

% parameters
% snrval = [10 20 50 100];
% n_std = 0.025*[5 2.5 1 0.5]; % noise standard deviation 0.025 ~ SNR 50
snrval = [50 100];
n_std = 0.025*[2 2]; % noise standard deviation 0.025 ~ SNR 50

rmean = 2.; %nm
width = 0.25; %nm
t_cut = [4 6 8 10 12 15 20 30]; % series of points to cut
t_cut = [6 6];

nstat = 2; % number of independant trials
n = 150;
q = 15;
extra = 10; % extra points to predict before zero
% parameters end
M = t_cut+extra;
hold off
tmin = 0; %us
tmax = 1; %us
% 2 ns steps
N = 250;
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

%generate noiseless 
tfree = linspace(tmin-extra*0.002,tmax,N+extra); %us
Sfree = dipolarsignal(tfree,r,P);

for q = 1:length(snrval)
    rmean_out = zeros(1,2*length(t_cut));
    sigma_out = zeros(1,2*length(t_cut));
    rmse_out = zeros(1,2*length(t_cut));
    for icut = 1:length(t_cut)
        list_rmean = zeros(1,nstat);
        list_sigma = zeros(1,nstat);
        list_rmse = zeros(1,nstat);
        for i = 1:nstat
            %Simulate dipolar evolution function
            Snoise = dipolarsignal(t,r,P,'noiselevel',n_std(q));
            Snoise = Snoise';
            Snoisecut = Snoise(t_cut(icut)+1:end);
            timecut = t(t_cut(icut)+1:end);

            % use the blp algorithm to find a back projection
            backpred = blp_epr(Snoisecut,n,q,M(icut));
            tnew = timecut(1)-stp*M(icut):stp:timecut(1)-stp;

            % concatenate the projection solution with the experimental data
            Vfull = [backpred Snoisecut];
%             Vfull = Vfull/max(Vfull);
            tfull = [tnew timecut];
            % best guess zero time needed for background 
            % ns -> us
            tfix = correctzerotime(Vfull, tfull);

            KB = dipolarkernel(tfix,r);

            % gaussian model fit
            maxGauss = 2;
            [Pfit,param,Nopt,metrics,Peval]  = fitmultigauss(Vfull,KB,r,maxGauss,'aic','Lower',[1.5 0.05]);
            Vfit = KB*Pfit;

            plot(tfull,Vfull,'k',tfree,Sfree,'r')
            hold on
            list_rmean(i) = param(1);
            list_sigma(i) = param(2);
            backlen = length(backpred);
            list_rmse(i) = sqrt((1/backlen)*sum((backpred-Sfree(1:backlen)').^2));

        end
        rmse_out(2*icut-1) = mean(list_rmse);
        rmse_out(2*icut) = std(list_rmse);
        rmean_out(2*icut-1) = mean(list_rmean);
        rmean_out(2*icut) = std(list_rmean);
        sigma_out(2*icut-1) = mean(list_sigma);
        sigma_out(2*icut) = std(list_sigma);
        fprintf('Points cut: %i \n', t_cut(icut))
        fprintf('rmse = %.4f +/- %.4f\n', mean(list_rmse), std(list_rmse))
        fprintf('r_mean = %.4f +/- %.4f\n', mean(list_rmean), std(list_rmean))
        fprintf('sigma_mean = %.4f +/- %.4f\n', mean(list_sigma), std(list_sigma))
    end
% 
%     formatStr = ['%d ' repmat('%f ', 1, length(rmean_out)) '\n'];
%     fid = fopen('./output/rmean.txt', 'a+');
%     fprintf(fid, formatStr, snrval(q), rmean_out);
%     fclose(fid);
% 
%     fid = fopen('./output/sigma.txt', 'a+');
%     fprintf(fid, formatStr, snrval(q), sigma_out);
%     fclose(fid);
%     
% 
%     fid = fopen('./output/timedomain_rmse.txt', 'a+');
%     fprintf(fid, formatStr, snrval(q), rmse_out);
%     fclose(fid);
end