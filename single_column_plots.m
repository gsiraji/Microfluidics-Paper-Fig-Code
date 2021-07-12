% single_column_plots.m generates figures 2 and 6- stochastic simulation of clogging for single-column devices 
close all; clear all;
syms z(x) t

% plot settings
markers = {'o','d','^'};
markercolors = {'r','b','m'};
marker_i = 0;
linecolors = {	'#A2142F'; '#0072BD'; '#7E2F8E';'k'};
linestyles = {'-','-.', '--', 'k:' };
line_i = 1;
linecolor = string(linecolors(line_i));
line_style = string(linestyles(line_i));

% set of rows and columns
row_set =[1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000];
column_set = 1;

% initialize vector of exact results for independent channels in fig  2
Ttheory = zeros(max(row_set)-1,2);

% iterate to simulate using 3 different clogging rate (lambda) functions
for lambdatype = 1:1:3
    
    % initialize the table of simulation data
    table1 = table;
   
    
    for k1 = row_set 
        for k2 = column_set
            
            c_mode = 1; % device conditions, 0 for consant flow, 1 for constant p
            
            Ntrials = 300; %300 number of trials
            N = 20000; %number of timesteps
            T = 10000; %final endpoint for time
            dt = T/N; %timestep size
            tvec = 0:dt:T; %time vector


            m = k1 %number of rows
            n = k2; %number of columns
            Q0 = 200; % initial total flow (for a constant total flow device) 
            P = 200; % initial pressure (for a constant pressure device)
            N0 = m*n;    % total number of channels
            e1 = (1/m);  % epsilon
            
            Gmat = ones(m,n); %matrix of resistances, scaled to 1 here
            Qmat = ones(m,n); % placeholder for the matrix of channel flows
            fvec  = zeros(N,1);  % fraction of clogged channels at each time step
            lamb = zeros(N,1); %lambda at each time step
            mask = ones(m,n); %matrix of channels, 0 for clogged channels, 1 for open channels
            Gsum = sum(mask.*Gmat,1); % column conductances
            G = 1/sum(Gsum.^-1); % net conductance
            
            if c_mode == 1
               % 'init flow is'
                Q0 =P*G; %cp
            else
             %   'init P is'
                P = Q0/G; %cf

            end
            Tclog = zeros(Ntrials,1); 
            fclog = zeros(Ntrials,1);
    
            for trial = 1:Ntrials
                mask = ones(m,n); %matrix of open channels, 1 open, 0 clogged


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
                     %create the matrix of probability of clogging for each
                     %channel according to the lambda function
                   if lambdatype == 1

                       probs = dt.*(Qmat./Q0).*(1.-Qmat./(2*Q0)).*mask;
                    else 
                        if lambdatype == 2

                            probs = dt.*(Qmat./Q0).*mask;
                        else

                            probs = dt*m.*(Qmat./Q0).^2.*mask;
                        end
                    end




                   mask(randmat<probs) = 0; %clog the channels according to their probability of clogging
                end
            end
            marker_i = marker_i +1;

            table2 = table([mean(Tclog),mean(fclog),std(Tclog),std(fclog),T,N,m,n,Ntrials,Q0]);
            table1 = [table1;table2];


        end
    end
    T =  table2array(table1(:,1));
    markertype = string(markers(lambdatype));
    markercolor = string(markercolors(lambdatype));
    markersz = 7;

    plot(row_set,T(:,1), markertype, 'MarkerSize', markersz, 'MarkerFaceColor', markercolor,'MarkerEdgeColor', markercolor )
    hold on;
    e = errorbar(row_set,T(:,1), T(:,3),'vertical','.','HandleVisibility','off', 'MarkerEdgeColor','none');
    hold on;
    
    for mi = 1:1:max(row_set)
        if c_mode == 1
            Q0 = P/n*mi;
            q = Q0/mi;
            if lambdatype == 1
            lamb_theo = 1/(mi)*(1-1/(2*mi));
            lambtheory = q/Q0*(1- q/(2*Q0));
            else

                if lambdatype == 2
                    lamb_theo = 1/mi;
                    %% Time to failure predicted by reliability engineering      
                    lambtheory = q/Q0;% prob of failure, constant
                     Rp = 1 - (1-exp(-lambtheory*t))^mi ;
                     R = Rp^n;
                     R = matlabFunction(R) ;
                    Ttheory(mi,2) = quadgk(R,0,inf);
                    Ttheory(mi,1) = mi; 
                else
                    lamb_theo = 1/mi;
                    lambtheory = mi.*(q/Q0).^2;% prob of failure, constant
                end
            end
            opts = odeset('stats', 'on');
            eqns = diff(z,x) == lamb_theo*(1-z);
            conds = z(0) == 0;
            S= dsolve(eqns, conds);
            Tmeanfield(mi,1) = mi;
            Tmeanfield(mi,2) = double(solve(S==((mi-1)/mi),x));
        else
            if lambdatype == 1
                lamb_theo = 1/((1-z)*mi)*(1-1/((1-z)*2*mi));
            else
                if lambdatype == 2
                    lamb_theo = 1/((1-z)*mi);
                else
                    lamb_theo = 1/(((1-z).^2)*mi);
                end
            end
        opts = odeset('stats', 'on');
        eqns = diff(z,x) == lamb_theo*(1-z);
        conds = z(0) == 0;
        S= dsolve(eqns, conds);
        Tmeanfield(mi,2) = double(solve(S==1,x));
        Tmeanfield(mi,1) = mi;
        end
    end
    
    plot(Tmeanfield(:,1),real((Tmeanfield(:,2))),line_style, 'lineWidth',3, 'color', linecolor)
    line_i = line_i +1;
    linecolor = string(linecolors(line_i));
    line_style = string(linestyles(line_i));
    hold on;




end

label1 = 'stochastic, \hspace{0.0005cm} $ \lambda = \alpha \tilde{q}(1-\tilde{q}/2)$';
label11 = 'mean-field, $ \lambda =  \alpha \tilde{q}(1-\tilde{q}/2)$';
label2 = 'stochastic, \hspace{0.0005cm} $ \lambda =  \alpha \tilde{q}$';
label22 = 'mean-field, $ \lambda =  \alpha \tilde{q}$';
label3 = 'stochastic, \hspace{0.0005cm} $ \lambda =  \alpha m \tilde{q}^2$';
label33 = 'mean-field, $ \lambda =  \alpha m \tilde{q}^2$';
label4 = 'exact,   \hspace{0.64cm}    $ \lambda = $ const. ';
if c_mode == 1
    plot(Ttheory(:,1),Ttheory(:,2), line_style, 'lineWidth', 3)
    legend(label1,label11,label2,label22,label3,label33, label4,'Interpreter','latex','Color','none','EdgeColor','none' );
else 
    legend(label1,label11,label2,label22,label3,label33,'Interpreter','latex','Color','none','EdgeColor','none' );
end
legend(  'Location', 'NorthWest' );
xlabel('$m$','Interpreter','latex');
ylabel('$ T_\textbf{clog} $ (a.u.) ','Interpreter','latex');
set(gca, 'fontsize', 16);
 set(gca, 'fontName','Times New Roman');
 set(gca, 'YScale', 'log')
 set(gca, 'XScale', 'log')
%  ylim([1 1600])