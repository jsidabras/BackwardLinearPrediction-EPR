% Test the blp by generating (a lot) of data to run as statistics. How good
% is this algorithm? 
% requires:
%   DEERLab v0.8
% 
% Author: Jason W. Sidabras (jason.sidabras@gmail.com)
% Initial writing: 04/05/2020 JWS
% New BLP algorithm: 16/06/2020 JWS
% GPLv3 License
clear, clf
clear all


% parameters
q = 25;
snrval = [10 20 50 100 500];
n_std = 0.025*[5 2.5 1 0.5 0.1]; % noise standard deviation 0.025 ~ SNR 50

rmean = 2.0; %nm
width = 0.15; %nm
t_cut = [4 6 8 10 12 15 20 30]; % series of points to cut

nstat = 500; % number of independant trials
extra = 10; % extra points to predict before zero
% parameters end
M = t_cut+extra;
hold off
tmin = 0; %us
tmax = 1; %us
% 2 ns steps
N = 500;
NtimeSteps = N;
%Time-domain axis
t = linspace(tmin,tmax,N); %us
stp = (t(2)-t(1));

rmin = 1; %nm
rmax = 5; %nm
N = 500;
%Time-domain axis
r = linspace(rmin,rmax,N); %nm

%Generate a distance distribution
P = rd_onegaussian(r,[rmean width]);

%generate noiseless 
tfree = linspace(tmin-extra*0.002,tmax,N+extra); %us
Sfree = dipolarsignal(tfree,r,P);

for qi = 1:length(snrval)
    fprintf('SNR value: %i \n', snrval(qi))
    rmean_out = zeros(1,2*length(t_cut));
    sigma_out = zeros(1,2*length(t_cut));
    rmse_out = zeros(1,2*length(t_cut));
    for icut = 1:length(t_cut)
        list_rmean = zeros(1,nstat);
        list_sigma = zeros(1,nstat);
        list_rmse = zeros(1,nstat);
        for i = 1:nstat
            %Simulate dipolar evolution function
            Snoise = dipolarsignal(tfree,r,P,'noiselevel',n_std(qi));
            Snoise = Snoise';
            Snoisecut = Snoise(M(icut)+1:end);
            timecut = tfree(M(icut)+1:end);
            Sfree = Snoise';

            % use the blp algorithm to find a back projection
            tnew = timecut(1)-stp*M(icut):stp:timecut(1)-stp;

            L = M(icut);

            % concatenate the projection solution with the experimental data
            [Vfull, backpred] = blp_epr(Snoisecut,L,q);
%             Vfull = Vfull/max(Vfull);
            tfull = [tnew timecut];
            % best guess zero time needed for background 
            % ns -> us
            tfix = correctzerotime(Vfull, tfull);

            KB = dipolarkernel(tfix,r);

            % gaussian model fit
            maxGauss = 1;
            [Pfit,param,Nopt,metrics,Peval]  = fitmultigauss(Vfull,KB,r,...
                maxGauss,'aic','Lower',[0.5 0.05],'Upper',[5. 0.75]);
            Vfit = KB*Pfit;

            list_rmean(i) = param(1);
            list_sigma(i) = param(2);
            backlen = length(backpred);
            list_rmse(i) = sqrt(mean((backpred-Sfree(1:backlen)).^2));
            hold on
plot(backpred)
plot(Sfree(1:backlen))
        end
        rmse_out(2*icut-1) = mean(list_rmse);
        rmse_out(2*icut) = std(list_rmse);
        rmean_out(2*icut-1) = mean(list_rmean);
        rmean_out(2*icut) = std(list_rmean);
        sigma_out(2*icut-1) = mean(list_sigma);
        sigma_out(2*icut) = std(list_sigma);
        fprintf('Points cut: %i \n', t_cut(icut))
%         fprintf('rmse = %.4f +/- %.4f\n', mean(list_rmse), std(list_rmse))
%         fprintf('r_mean = %.4f +/- %.4f\n', mean(list_rmean), std(list_rmean))
%         fprintf('sigma_mean = %.4f +/- %.4f\n', mean(list_sigma), std(list_sigma))
    end

    formatStr = ['%d ' repmat('%f ', 1, length(rmean_out)) '\n'];
    fmeanfn = [ './output/DEERrmean-' num2str(rmean) 'nm-' num2str(width) 'sigma-' num2str(NtimeSteps) 'pts' num2str(tmax) 'us2.csv'];
    fid = fopen(fmeanfn, 'a+');
    fprintf(fid, formatStr, snrval(qi), rmean_out);
    fclose(fid);

    fsigmafn = [ './output/DEERsigma-' num2str(rmean) 'nm-' num2str(width) 'sigma-' num2str(NtimeSteps) 'pts' num2str(tmax) 'us2.csv'];
    fid = fopen(fsigmafn, 'a+');
    fprintf(fid, formatStr, snrval(qi), sigma_out);
    fclose(fid);
    
    tdrmsefn = [ './output/timedomain-RMSE-' num2str(rmean) 'nm-' num2str(width) 'sigma-' num2str(NtimeSteps) 'pts' num2str(tmax) 'us2.csv'];
    fid = fopen(tdrmsefn, 'a+');
    fprintf(fid, formatStr, snrval(qi), rmse_out);
    fclose(fid);

end