% clog_sim.m generates the sheets in the worksheet files 'cp_n.xlsx' and
% 'cf_n.xlsx'. These sheets can later be used in 'cf_n_plotter'...
% and 'cp_n_plotter' that generate figures 3 and 7. This file performs stochastic simulation of clogging for constant total flow and constant
%pressure difference. Collects the clogging time and fraction of clogged channels at the time of clogging and exports them in an Excel Spreadsheet  


close all; clear all;
syms x z(x)




markers = {'o','d','^'};
marker_i = 0;
for k1 = 256% [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500 ] % 
    for k2 =  [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100] %plot for different columns and rows ,512,2048
        
        
        Ntrials = 300; %number of trials
        N = 10000; %number of timesteps
        T = 2000; %final endpoint for time
        dt = T/N; %timestep size
        tvec = 0:dt:T; %time vector


        m = k1; %number of rows
        n = k2 %number of columns
        Q0 = 200; % initial total flow (for a constant total flow device) 
        P = 200; % initial pressure (for a constant pressure device)
        N0 = m*n;    % total number of channels
        e1 = (1/m);  % epsilon
% cloggedN = zeros(N,1);
        Gmat = ones(m,n); %matrix of conductances, scaled to 1 here
        Qmat = ones(m,n); % placeholder for the matrix of channel flows
%  masksumvec = ones(N,1);
        fvec  = zeros(N,1);  % fraction of clogged channels at each time step
        lamb = zeros(N,1); %lambda at each time step
        mask = ones(m,n); %matrix of channels, 0 for clogged channels, 1 for open channels
        Gsum = sum(mask.*Gmat,1); % column conductances
        G = 1/sum(Gsum.^-1); % net conductance
        c_mode = 0; % device conditions, 0 for consant flow, 1 for constant p
        Tclog = zeros(Ntrials,1); %
        fclog = zeros(Ntrials,1);
        if c_mode == 1
           % 'init flow is'
            Q0 =P*G; %cp
        else
         %   'init P is'
            P = Q0/G; %cf

        end
        for trial = 1:Ntrials
            mask = ones(m,n); %matrix of clogged channels
% if (mod(trial-1,20)==0)
%     trial
% end
 

            for step=1:N
   
               Gsum = sum(mask.*Gmat,1); %recalculate the conductance for each column
               G = 1/(sum(1./Gsum)); %recalculate the conductance for the device
               fvec(step) = (N0 - nnz(mask))/N0; %calculate the fraction of clogged channels at this time step


                if find(all(mask==0))>0 %  see if there is a column that is completely clogged

                    Tclog(trial) = dt*(step-1); % store time to complete clog
                    fclog(trial) = fvec(step);  % store the fraction of clogged channels at the time of clogging

                    break
                end
                

                if c_mode == 1 %constant pressure
                    for column_i = 1:1:n
                        Qmat(1:m,column_i) = G*P./repmat(nnz(mask(1:m,column_i)),m,1);
                    end
                else
                    for column_i = 1:1:n %constant flow
                        Qmat(1:m,column_i) = Q0./repmat(nnz(mask(1:m,column_i)),m,1);
                    end
                end
                q0 = Q0/m;


                fvec(step) = (N0 - nnz(mask))/N0; %store the number of clogged channels at this time step

               randmat = rand(m,n); %create a random matrix

%                  sum_lamb = sum((Qmat/q0).*mask)./sum(mask);
            %         sum_lamb = sum(m.*(Qmat/q0).^2.*mask)./sum(mask);
%                     sum_lamb = sum((Qmat/q0).*mask.*(1.-Qmat./(2*Q0)))./sum(mask);

%                  probs = dt.*((Qmat/q0).*mask);
%                probs = dt.*(Qmat./Q0).*mask;
%                probs = dt.*m.*(Qmat./Q0).^2.*mask;
               probs = dt.*(Qmat/(Q0)).*(1.-Qmat./(2*Q0)).*mask; %create the matrix of probability of clogging for each channel
%                 probs = dt.*(Qmat/(q0)).*(1.-Qmat./(2*Q0)).*mask;
            %    lamb(step) = mean(sum_lamb);

               mask(randmat<probs) = 0; %clog the channels according to their probability of clogging
            end
        end
    marker_i = marker_i +1;
    AA=strcat('A',num2str(marker_i),':J',num2str(marker_i));
     table2 = table([mean(Tclog),mean(fclog),std(Tclog),std(fclog),T,N,m,n,Ntrials,Q0]);
        %uncomment the following line (according to the choice of lambda)
        %to store the clogging time data
     
%      writetable(table2, 'cf_n.xlsx','Sheet','q1-Q0','Range',AA,'WriteVariableNames', false)
%      writetable(table2, 'cf_n.xlsx','Sheet','q2Q0','Range',AA,'WriteVariableNames', false)
%      writetable(table2, 'cf_n.xlsx','Sheet','qQ0','Range',AA,'WriteVariableNames', false)
     

    end
end
















