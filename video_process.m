% video_process.m generates figure 6. To create the exact figure, ...
% the following 2 files must be available: ... 
% 'vid_tracking_data.xlsx' manually tracked f data and 'du.avi' the
% experiment video. The link to get the table of manually tracked data is
% https://brandeis.box.com/s/3fbqlsslyryjjeey56u8kf3lwxiefvj0 with password
% brightness. Your operating system may need codec to read the .avi file.


close all; clear all;

prompt1 = 'Did you download the data table? ["y"/"n"]';
answer1 = input(prompt1);
if answer1 == 'n'
    error('Please download the data table from the links in this code before running this file')
end






% read the brightness/f data manually collected
plot_data = readtable('vid_tracking_data.xlsx','Sheet',1);
vtime = table2array(plot_data(:,1)); % store the time data
% download the data from url2 =
% 'https://brandeis.box.com/s/3fbqlsslyryjjeey56u8kf3lwxiefvj0' before
% running this line (pw: brightness)
tracked = table2array(plot_data(:,2)); % store the manually tracked f. 




% download the video from the online source. If the link is broken, check https://www.pnas.org/content/112/5/1422/tab-figures-data, video 2 for the latest link
url = "https://static-movie-usa.glencoesoftware.com/source/10.1073/884/ba65404748f24af2a4b39f93f759f7bd00e9c81a/pnas.1424111112.sm02.avi";
filename = 'du.avi';
outfilename = websave(filename, url);



%read video and extract the normalized brightness slope
vid = VideoReader('du.avi');
% get the length of the video in seconds
info = mmfileinfo('du.avi');
vid_duration = info.Duration;
% clogging occurs from time 20.7666s to time 37.5
cloggingstarts = 20.7666;
cloggingends = 37.5;
% max fraction of clogged channel is 7/10
max_f = 0.7; 
% device characteristics
alpha = 4e15;
mu = 4e-3;
L = 15e-6;
h = 4e-6;
w = 5e-6;
% channel resistance
r = 12*mu*L/(w*h^3*(1-0.63*h/w));
% device pressure 
P = 222;
% P = 1;
%% read brightness of frames

% choose the x-y coordinates of the line profile
x = [685 685];
y = [1 264];


% iterate over the frames to collect line profile data

% initiate the index
step = 0;
% initiate the background brightness
bg = 0;
while hasFrame(vid)
    % update the index
    step = step +1;
	% read video frame
	vidFrame = readFrame(vid);
    % set the line profile
    linepr = improfile(vidFrame, x,y);
    % calculate the total brightness
    bv = sum(linepr);
    % store the relative line profile birghtness of the frame  
    bvec(step,1) = bv(:,:,1) - bg;
    % calculate the background brightness
    if step == 1
        bg = bvec(step,1);
        bvec(step,1) = 0;
    end

    
end

% 623-1125

%% calculate the slope in figure 6.B

% calculate the conversion rate for time to frame number
t2frame = step/vid_duration;
% calculate the initial frame number
iframe = round(cloggingstarts*t2frame);
% calculate the final frame number
fframe = round(cloggingends*t2frame);
% normalize the brightness 
normbr = -max_f.*bvec(iframe:fframe)/(max(abs(bvec(iframe:fframe))));
% create the time vector
timeframe = linspace(0,60,step);
normt = linspace(0,(37.5-timeframe(1,623)),503)';
vtime = vtime - timeframe(1,623).*ones(length(vtime),1);

%% Fit: 'cb fit 1'.
% cb = createFit(tracked, vtime); % cb=0.0020; %from tracked data
cb = createFit(normbr, normt); % cb=0.001701; %from brightness data


%% do Ntrials trials of the stochastic simulation to create fig 5.c
Ntrials = 1000; % set the number of trials
N = 50000; %number of timesteps

T = 500; %set the total simulation time
dt = T/N; % calculate the time-step length
tvec = 0:dt:T;  % create the sim. time vector

% the device properties
m=10; % set the number of rows
n=10; %set the number of columns

% the total number of channels for column 1
N0 = m;   

fvec_many  = ones(N+1,Ntrials); % initilize the f vector for all trials
Tplot = zeros(Ntrials,1); % initilize the vector of clogging times


for trial = 1:Ntrials
    fvec  = ones(N+1,1); % initilize the vector of fraction of clogged channels for each trial
    mask = ones(m,n); % matrix of open channels, 1 open, 0 blocked
    Qmat = ones(m,n);
    fvec(1) = (N0 - nnz(mask(:,1)))/N0; % calculate the initial f
%     lamb = zeros(N,1); %
    Gmat = ones(m,n);
    Gsum = sum(mask.*Gmat,1); % column conductances
	G = 1/sum(Gsum.^-1); % net conductance
    % at each time step, iterate until a column is clogged
    for step=1:N
        Gsum = sum(mask.*Gmat,1); %recalculate the conductance for each column
        G = 1/(sum(1./Gsum)); %recalculate the conductance for the device
        % stop the trial when one column is clogged
        if find(all(mask(:,1)==0))>0  
            Tplot(trial) = dt*(step-1); % store the clogging time
            break
        end
  
        for column_i = 1:1:n
            Qmat(:,column_i) = G*P./repmat(nnz(mask(:,column_i)),m,1);
        end
    
    
   
 
   % calculate the clogging rate lambda at this time step
%         lamb(step) = c1/(1-dt*step*c1);
        lamb = zeros(m,n);
        lamb(:,1) = (alpha*cb/r).*Qmat(:,1);
%         lamb(:,1) = a.*Qmat(:,1);
%         if lamb(step) < 0 . 
%             lamb(step) = c1;
%         end
        
        probs = dt.*lamb.*mask; % calculate the clogging probablity
        % clog a channel if the clogging probability is greater than the
        % random number for a channel
        randmat = rand(m,1);
        mask1 = mask(:,1);
        mask1(randmat<probs) = 0;
        mask(:,1) = mask1;
        % store the fraction of clogged channels for column 1
        fvec(step+1) = (N0 - nnz(mask(:,1)))/N0;
    end
    fvec_many(:, trial) = fvec; % store the vector of fraction of clogged channels for this trial

end

mean_fvec_many = [];
std_fvec_many = [];

for j=1:1:N
    mean_fvec_many = [mean_fvec_many,mean(fvec_many(j,:))]; % calculate the mean of fvec

    std_fvec_many = [std_fvec_many,std(fvec_many(j,:))]; % calculate the standard deviation over the trials
end
% solve the mean-field equation
syms z(t)
lambda_theory = (alpha*cb*P/r)./((n-1).*(1-z)+1); % set the clogging rate function .
eqn = diff(z,t) == lambda_theory*(1-z); % define the differential equatio
conds = z(0) == 0; % set the initial conditions
Sol = matlabFunction(dsolve(eqn , conds)); % solve the ODE
%a*(1-lambertw(-9*exp(-9-0.0731*x))/9) -- 1 - lambertw(0, -9*exp(log(-exp(9)) - (5269856008713839*t)/72057594037927936))/9

%% set up 
curve1 = mean_fvec_many + std_fvec_many;
curve2 = mean_fvec_many - std_fvec_many;
x2 = [tvec(1:nnz(mean_fvec_many)), fliplr(tvec(1:nnz(mean_fvec_many)))];
inBetween = [curve1(1:nnz(mean_fvec_many)), fliplr(curve2(1:nnz(mean_fvec_many)))];
% flip the video thumbnail to show cells flowing from left to right
vidFrame = vidFrame(:,end:-1:1,:);



%% plot
tiles = tiledlayout(2,2,'TileSpacing','none','Padding', 'compact');
% tile for fig 6A
nexttile(1, [1 2])
imshow(vidFrame);
% draw a line to demonstrate the line profile
annotation('line',[0.153 0.153],[0.966 0.56],'Color',[1 1 0.06],'LineWidth',3);
% tile for fig 6B 
nexttile
plot(vtime, tracked, 'o','MarkerFaceColor','k','MarkerEdgeColor','k', 'displayName','manually tracked', 'MarkerSize', 6)
hold on; plot(normt, normbr, '^','MarkerFaceColor','none','MarkerEdgeColor','r','displayName','calculated from brightness', 'MarkerSize', 4 )
xlim([0 25])
ylim([0 1])
set(gca, 'fontsize', 12, 'fontname', 'Times New Roman');
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$f$', 'interpreter', 'latex')
legend('location', 'northwest','Units','Pixels', 'color', 'none', 'edgeColor','none','fontSize',12)
% tile for fig 6C
nexttile
fill(x2, inBetween, [0.92,0.88,0.94], 'displayName', 'standard deviation', 'EdgeColor', [0.92,0.88,0.94]);
hold on;
plot(tvec(1:nnz(mean_fvec_many)), mean_fvec_many(1:nnz(mean_fvec_many)),':','LineWidth',3, 'displayName', 'stochastic mean', 'color', [94,60,153]./256)
hold on;
plot(tvec(1:step+400), fvec(1:step+400),'-.','LineWidth',2,'DisplayName', 'sample trial', 'color', [230,97,1]./256)
hold on;
plot(tvec, Sol(tvec),'k','LineWidth',2,'displayName', 'mean-field')
xlim([0 60])
ylim([0 1])
set(gca, 'fontsize', 12, 'fontname', 'Times New Roman');
xlabel('$t$ (s)', 'interpreter', 'latex')
legend('location', 'southeast', 'color', 'none','Units','Pixels', 'edgeColor','none','fontSize',12)%'location', 'northwest',

% Create the label for fig 6A
annotation('textbox',...
    [0.0258756218905473 0.974137931034483 0.0404593698175788 0.0431034482758617],...
    'String',{'A'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create the label for fig 6B
annotation('textbox',...
    [0.0209004975124378 0.5 0.0553847429519071 0.0646551724137931],...
    'String',{'B'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create the label for fig 6C
annotation('textbox',...
    [0.506804311774462 0.5 0.0553847429519072 0.0646551724137926],...
    'String','C',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');


function cb = createFit(normbr, normt)
%CREATEFIT
%  Create a fit.
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%  Auto-generated by MATLAB 
%% Fit: 

% normbrn = (r/(P*alpha)).*(9.*normbr - log(1 - normbr));
[xData, yData] = prepareCurveData(  normbr, normt);

% Set up fittype and options.
ft = fittype( '(0.0051/cb)*(9*x-log(1-x))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.504913705185244;

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, opts );

% Extract fit coefficients
cb = coeffvalues(fitresult); 


end
