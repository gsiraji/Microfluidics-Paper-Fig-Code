% vid_process.m generates figure 5. To create the exact figure, ...
% the following 2 files must be available: ... 
% 'vid_data.xlsx' manually tracked f data and 'du.avi' the experiment video,


close all; clear all;

% read the brightness/f data manually collected
plot_data = readtable('vid_data.xlsx','Sheet',1);
vtime = table2array(plot_data(:,1)); % store the time data
tracked = table2array(plot_data(:,2)); % store the manually tracked f 

% download the video from the online source. If the link is broken, check https://www.pnas.org/content/112/5/1422/tab-figures-data for the latest link
url ="https://static-movie-usa.glencoesoftware.com/source/10.1073/884/ba65404748f24af2a4b39f93f759f7bd00e9c81a/pnas.1424111112.sm02.avi";
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

%% calculate the slope in figure 5.b

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
% create a fit to estimate the slop of the normalized brightness 
mdl = fitlm(timeframe(iframe:fframe),normbr,'linear');
% extract the slope coefficient
c1 = table2array(mdl.Coefficients(2,1)); 


%% do Ntrials trials of the stochastic simulation to create fig 5.c


Ntrials = 1000; % set the number of trials
N = 5000; %number of timesteps

T = 200; %set the total simulation time
dt = T/N; % calculate the time-step length
tvec = 0:dt:T;  % create the sim. time vector

% the device properties
m=10; % set the number of rows
n=1; %set the number of columns

% the total number of channels
N0 = m*n;   

fvec_many  = ones(N,Ntrials); % initilize the f vector for all trials
Tplot = zeros(Ntrials,1); % initilize the vector of clogging times


for trial = 1:Ntrials
    fvec  = ones(N,1); % initilize the vector of fraction of clogged channels for each trial
    mask = ones(m,n); % matrix of open channels, 1 open, 0 blocked
 
    fvec(1) = (N0 - nnz(mask))/N0; % calculate the initial f
    lamb = zeros(N,1); %
    
    % at each time step, iterate until a column is clogged
    for step=1:N
        % stop the trial when one column is clogged
        if find(all(mask==0))>0  
            Tplot(trial) = dt*(step-1); % store the clogging time
            break
        end
  
   
    
    
   
 
   % calculate the clogging rate lambda at this time step
        lamb(step) = c1/(1-dt*step*c1);
        if lamb(step) < 0
            lamb(step) = c1;
        end
        
        probs = dt*lamb(step).*mask; % calculate the clogging probablity
        % clog a channel if the clogging probability is greater than the
        % random number for a channel
        randmat = rand(m,n);
        mask(randmat<probs) = 0;
        % store the fraction of clogged channels
        fvec(step+1) = (N0 - nnz(mask))/N0;
    end
    fvec_many(:, trial) = fvec; % store the vector of fraction of clogged channels for this trial

end

mean_fvec_many = [];
std_fvec_many = [];

for j=1:1:N
    mean_fvec_many = [mean_fvec_many,mean(fvec_many(j,:))]; % calculate the mean of fvec

    std_fvec_many = [std_fvec_many,std(fvec_many(j,:))]; % calculate the standard deviation over the trials
end

%% set up 
curve1 = mean_fvec_many + std_fvec_many;
curve2 = mean_fvec_many - std_fvec_many;
x2 = [tvec(1:nnz(mean_fvec_many)), fliplr(tvec(1:nnz(mean_fvec_many)))];
inBetween = [curve1(1:nnz(mean_fvec_many)), fliplr(curve2(1:nnz(mean_fvec_many)))];
% flip the video thumbnail to show cells flowing from left to right
vidFrame = vidFrame(:,end:-1:1,:);

normt = linspace(0,(37.5-timeframe(1,623)),503);
vtime = vtime - timeframe(1,623).*ones(length(vtime),1);

%% plot
tiles = tiledlayout(2,2,'TileSpacing','none','Padding', 'compact');
% tile for fig 5A
nexttile(1, [1 2])
imshow(vidFrame);
% draw a line to demonstrate the line profile
annotation('line',[0.153 0.153],[0.966 0.56],'Color',[1 1 0.06],'LineWidth',3);
% tile for fig 5B 
nexttile
plot(vtime, tracked, 'o','MarkerFaceColor','k','MarkerEdgeColor','k', 'displayName','manually tracked', 'MarkerSize', 6)
hold on; plot(normt, normbr, '^','MarkerFaceColor','none','MarkerEdgeColor','r','displayName','calculated from brightness', 'MarkerSize', 4 )
xlim([0 25])
ylim([0 1])
set(gca, 'fontsize', 12, 'fontname', 'Times New Roman');
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$f$', 'interpreter', 'latex')
legend('location', 'northwest','Units','Pixels', 'color', 'none', 'edgeColor','none','fontSize',12)
% tile for fig 5C
nexttile
fill(x2, inBetween, [0.925 0.91 0.91], 'displayName', 'standard deviation', 'EdgeColor', [0.925 0.9 0.9]);
hold on;
plot(tvec(1:nnz(mean_fvec_many)), mean_fvec_many(1:nnz(mean_fvec_many)),'LineWidth',3, 'displayName', 'mean')
hold on;
plot(tvec(1:step+400), fvec(1:step+400),'-.','LineWidth',2,'DisplayName', 'sample trial', 'color', [0.7,0.7,0.7])
xlim([0 35])
ylim([0 1])
set(gca, 'fontsize', 12, 'fontname', 'Times New Roman');
xlabel('$t$ (s)', 'interpreter', 'latex')
legend('location', 'northwest', 'color', 'none','Units','Pixels', 'edgeColor','none','fontSize',12)%'location', 'northwest',

% Create the label for fig 5A
annotation('textbox',...
    [0.0258756218905473 0.974137931034483 0.0404593698175788 0.0431034482758617],...
    'String',{'A'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create the label for fig 5B
annotation('textbox',...
    [0.0209004975124378 0.5 0.0553847429519071 0.0646551724137931],...
    'String',{'B'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create the label for fig 5C
annotation('textbox',...
    [0.506804311774462 0.5 0.0553847429519072 0.0646551724137926],...
    'String','C',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

