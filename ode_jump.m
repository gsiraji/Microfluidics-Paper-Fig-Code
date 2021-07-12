% ode_jump.m generates figures 4 and 8 by running a stochastic simulation of clogging for
% 3 devices with different number of columns and solving the ODE system numerically
clear all
syms z(x)  %set up the symbolic function for the mean-field equation

n = 10; m = 256; %set the number of columns n and rows m
yjump = 1/n; % set the magnitude of the stochastic jump
cp_mode = 1; % choose constant pressure (1) or constant total flow (0)


% solve the system of equations numerically with the stochastic jump
% imposed
if cp_mode == 1
    tjump = 500; % set the time at which the jump occurs
    tfin = 1500; % set the final time for the numerical solver
    [t1, y1] = ode45(@cp,[0 tjump],zeros(n,1)) ; %solve the ode up to the jump
    tsz = length(t1);
    [t2, y2] = ode45(@cp,[tjump tfin],[y1(tsz,1)+yjump, y1(tsz,2:n)]) ; %solve the ode after the jump
   
    % solve the symbolic solution to the meanfield equation
    lambda_theory = 1/m; % set lambda
    eqn = diff(z,x) == lambda_theory*(1-z); % set the differential equation
    conds = z(0) == 0; % pick the initial conditions
    Sol = matlabFunction(dsolve(eqn , conds)); %solve the ode symbolically
else 
    tjump = 100; % set the time at which the jump occurs
    tfin = 230; % set the final time for the numerical solver
    [t1, y1] = ode23(@cf,[0 tjump],zeros(n,1)) ; %solve the ode until the jump
    tsz = length(t1); %store the size of the time vector
    [t2, y2] = ode23(@cf,[tjump tfin],[y1(tsz,1)+yjump, y1(tsz,2:n)]) ; %solve the ode after the jump
    
    % solve the symbolic solution to the meanfield equation
    lambda_theory = 1/(m.*(1-z)); % set the clogging rate function  
    eqn = diff(z,x) == lambda_theory*(1-z); % define the differential equatio
    conds = z(0) == 0; % set the initial conditions
    Sol = matlabFunction(dsolve(eqn , conds)); % solve the ODE
end

markersizes = [6,5.5,5.5]; %set the marker sizes for the plot
markercolors = [[0.26 0.265 0.275];[0.175 0.76 0.84];[1 0.6 0.5]]; %set the marker colors for the plot
markers = {'o','d','^'}; % select the markers for the plot
marker_i = 0; % set the marker index

tilez = tiledlayout(1,2);   % set the tiles for the plot
nexttile 

for k2 =  [1,10,256] % create the plot for different number of columns 
        
        c_mode = cp_mode; % device conditions, 0 for consant flow, 1 for constant p
        Ntrials = 1; %number of trials
        if c_mode == 1 
        N = 2800; %number of timesteps
        T = 1600; %final endpoint for time
        else
            N = 1000; %number of timesteps
            T = 300; %final endpoint for time
        end
        dt = T/N; %timestep size
        tvec(:,1) = 0:dt:T; %time vector

        n = k2 % set the number of columns
        Q0 = 100; % set the initial total flow (for a constant total flow device) 
        P = 200; % set the initial pressure (for a constant pressure device)
        N0 = m*n;    % calculate the  total number of channels


        Gmat = ones(m,n); %initiate the matrix of conductances, scaled to 1 here
        Qmat = ones(m,n); %initiate the  placeholder for the matrix of channel flows

        fvec  = zeros(N,1);  %initiate the vector of fraction of clogged channels at each time step
        lamb = zeros(N,1); %initiate the vector of lambda at each time step
        mask = ones(m,n); %initiate the matrix of channels, 0 for clogged channels, 1 for open channels
        Gsum = sum(mask.*Gmat,1); % calculate the column conductances
        G = 1/sum(Gsum.^-1); % calculate net conductance
        
        Tclog = zeros(Ntrials,1); % initiate the vector of time of clogging
        fclog = zeros(Ntrials,1); %initiate the vector of fraction of clogged channels at the time of clogging
        if c_mode == 1
           % 'init flow is'
            Q0 =P*G; % calculate the flow for a constant pressure device
        else
         %   'init P is'
            P = Q0/G; % calculate the pressure for a constant flow device

        end
        for trial = 1:Ntrials
            mask = ones(m,n); %reset the matrix of clogged channels at each trial
 
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
                        Qmat(1:m,column_i) = G*P./repmat(nnz(mask(1:m,column_i)),m,1); %calculate and set the channel flow
                    end
                else
                    for column_i = 1:1:n %constant flow
                        Qmat(1:m,column_i) = Q0./repmat(nnz(mask(1:m,column_i)),m,1);%calculate and set the channel flow
                    end
                end
                q0 = Q0/m;


               fvec(step) = (N0 - nnz(mask))/N0; %store the number of clogged channels at this time step

                
               randmat = rand(m,n); %create a random matrix
               probs = dt.*(Qmat./Q0).*mask; %create the matrix of clogging probabilities for the channels
               mask(randmat<probs) = 0; %clog the channels according to their probability of clogging
               
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    marker_i = marker_i +1;

     plot([tvec(1:50:step); tvec(step,1)], [fvec(1:50:step); fvec(step,1)],string(markers(marker_i)),'MarkerFaceColor',[markercolors(marker_i,:)],'MarkerEdgeColor','none','DisplayName', ['$ n = $ ' num2str(n)], 'MarkerSize', markersizes(marker_i))
     hold on;
    legend( 'Interpreter','latex','Color','none','EdgeColor',[0.5 0.5 0.5] );
    legend(  'Location', 'southEast' , 'FontSize',12 );
    xlabel('$t$ (a.u.)','Interpreter','latex');
    ylabel('$ f $','Interpreter','latex');
    set(gca, 'fontsize', 12);
    set(gca, 'fontName','Times New Roman');
    if c_mode == 1
        xlim([0 1500])
        else
    xlim([0 258])
    end
    ylim([0 1])
    axis 'square'
end



plot(tvec, Sol(tvec),'k','LineWidth',2,'displayName', 'mean-field')


%% f plots from the ODE

nexttile




plot([t1; t2], [y1(:,2);y2(:,2)],'-','color',[0.8500 0.3250 0.0980],'displayName', 'other columns', 'lineWidth', 5) %,'MarkerFaceColor',[0.175 0.76 0.84],'MarkerEdgeColor','none','MarkerSize',12)
hold on;
plot([t1; t2], [mean(y1, 2);mean(y2, 2)],'--y','displayName', 'all columns', 'lineWidth', 4)%,'MarkerFaceColor',[1 0.6 0.5],'MarkerEdgeColor','none','MarkerSize',9)
hold on;
plot(tvec, Sol(tvec),'k','LineWidth',2,'displayName', 'mean-field')
hold on;
plot([t1(1:5:length(t1)); t2(1:5:length(t2))], [y1((1:5:length(t1)),1);y2((1:5:length(t2)),1)],':s','color',[0.3 0.4 0.4],'lineWidth',1.5,'displayName', 'column 1','MarkerFaceColor',[0.3 0.4 0.4],'MarkerEdgeColor',[0.3 0.4 0.4],'MarkerSize',6)


%  
legend( 'Interpreter','latex','Color','none','EdgeColor',[0.5 0.5 0.5] );
leg = legend(  'Location', 'southEast' , 'FontSize',12 );
leg.ItemTokenSize = [50,50,50,50];
xlabel('$t$ (a.u.)','Interpreter','latex');
ylabel('$ f $','Interpreter','latex');
if c_mode == 1
    xlim([0 1500])
else
    xlim([0 258])
end
ylim([0 1])
set(gca,'FontName','Times New Roman','FontSize',12);
axis 'square'





% remove the extra white space
tilez.TileSpacing = 'none';
tilez.Padding = 'compact';

%% differential equation functions
function dy = cp(t,y)
m = 256;
alpha = 1;
r = 1;
P = 200;
n = 10;
Q0 = m.*P/(r.*n);
q = ones(length(y),1);
dy = ones(length(y),1);
Ry = r./(m - m.*y);
R = sum(Ry, 'all');
for i=1:1:length(y)
    q(i,1) = (P/R)./(m - m.*y(i));
    dy(i,1) = alpha.*q(i,1)./Q0.*(1 - y(i,1));
end
end

function dy = cf(t,y)
m = 256;
alpha = 1;
Q0 = 100;
q = ones(length(y),1);
dy = ones(length(y),1);
for i=1:1:length(y)
    q(i,1) = Q0./(m - m.*y(i));
    dy(i,1) = alpha.*q(i,1)./Q0.*(1 - y(i,1));
end
end