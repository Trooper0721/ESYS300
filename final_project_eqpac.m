% Clear all variables and close all figures

clearvars, close all
%%
% Define path

%User = 'Tremblay'
User = 'Student'
%User = 'Planat'

if strcmp(User, 'Tremblay')
    root    = ['/Users/tremblay/Dropbox/Teaching/ESYS-300/'] 
elseif strcmp(User, 'Planat')
    root    = ['./']
elseif strcmp(User, 'Student')
    root    = ['C:\Users\titan\OneDrive\McGill\ESYS_300\lab3\']       % [ENTER your own path]
end

MatfilePath = [root 'Matfiles'] ; 
DataPath    = [root 'Data'] ;
PlotPath    = [root 'Plot'] ;

addpath(MatfilePath,DataPath, PlotPath) 


%%
% NAME:        

GivenName   = ['Evan'] ;                   % ENTER given name
FamilyName  = ['Henderson'] ;                  % ENTER family name

disp(' ')
disp(['Name: ', GivenName, ' ', FamilyName])    % Print Full Name
disp(' ')
%%
disp(' ')
disp('Do a scatter plot of Global Jan-Feb-Mar Mean temperature anomaly')
disp('as a function of the Global Yearly Mean temperature anomaly from 1880 to 2015')
disp(' ')
%
% 0 winter JFM, 1 spring AMJ, 2 summer JAS, 3 fall OND
iter_num = 100; %control number of iterations
ssn_num = 4;    %control seasonal devisions, currently only works for 4 seasons
doplot = false;
percent = 0.8;
divide = 1/percent;
stds = zeros([ssn_num iter_num]);
for season = 1:ssn_num
for iterate = 1:iter_num
    disp(iterate)
season_names = ["Winter (JFM)" "Spring (AMJ)" "Summer (JAS)" "Fall (OND)"];
season_name = season_names(season);
% Load GISS-Temperature data: 1880 to today

%ncdisp("gistemp1200_GHCNv4_ERSSTv5.nc")

time = ncread("gistemp1200_GHCNv4_ERSSTv5.nc",'time');
lat = ncread("gistemp1200_GHCNv4_ERSSTv5.nc",'lat');
lon = ncread("gistemp1200_GHCNv4_ERSSTv5.nc",'lon');
time_bnds = ncread("gistemp1200_GHCNv4_ERSSTv5.nc",'time_bnds');
tempanomaly = ncread("gistemp1200_GHCNv4_ERSSTv5.nc",'tempanomaly');
%%
% Calculate the means 
tempanomaly_eqpac = tempanomaly(1:40,43:48,:); %focus on equatorial pacific
tempanomaly = tempanomaly(1:40,43:48,:);

mean_mat_eqpac = mean(tempanomaly_eqpac,"omitmissing"); %mean the lons
pre_mean_mat_eqpac = squeeze(mean_mat_eqpac);   %squeeze out the singleton dimension
monthly_mean_mat_eqpac = mean(pre_mean_mat_eqpac,"omitmissing"); %mean the lats

mean_mat = mean(tempanomaly,"omitmissing"); %mean the lons
pre_mean_mat = squeeze(mean_mat);   %squeeze out the singleton dimension
monthly_mean_mat = mean(pre_mean_mat,"omitmissing");      %mean the lats
month_number = 1:1728;
%%
monthly_mean_mat_eqpac = [monthly_mean_mat_eqpac,NaN,NaN,NaN,NaN,NaN];
monthly_mean_mat_eqpac = monthly_mean_mat_eqpac((season-1)*3 + 1:end-(12-(season-1)*3));

monthly_mean_mat = [monthly_mean_mat,NaN,NaN,NaN,NaN,NaN];
monthly_mean_mat = monthly_mean_mat((season-1)*3 + 1:end-(12-(season-1)*3));
yearly_mean_mat = reshape(monthly_mean_mat,12,143);
yearly_mean_mat = squeeze(mean(yearly_mean_mat,"omitmissing"));
year_number = 1880:2022;
%%
%do random selection stuff
N = length(year_number);
P = randperm(N);
cut = floorDiv(N,divide);
P = P(1:cut);
%%
% and seasonal means (do not show your answer)
season_add = [1,2,3];
season_add = season_add + (season-1)*3;
seasonID = [];
for i = 0:142
    seasonID = [seasonID, (season_add+(12*i))];
end    
%%
seasonal_means_eqpac = monthly_mean_mat_eqpac(seasonID);
seasonal_means_eqpac = reshape(seasonal_means_eqpac,3,143);
seasonal_means_eqpac = squeeze(mean(seasonal_means_eqpac,"omitmissing"));

seasonal_means = monthly_mean_mat(seasonID);
seasonal_means = reshape(seasonal_means,3,143);
seasonal_means = squeeze(mean(seasonal_means,"omitmissing"));

% Plot the yearly mean temperature time series just to be safe (figure(1))
if doplot
figure(1)
plot(year_number,yearly_mean_mat)
title('Yearly Global Mean Temperature')
xlabel('Year')
ylabel('Global Mean Temperature Anomaly (K)')


%%
% Do the scatter plot
figure(1+season)
plot(year_number,seasonal_means_eqpac,'-')
hold on
plot(year_number,yearly_mean_mat,'-')
hold on
title(['Global '+season_name+' Seasonal Mean temperature anomaly as a function of the Global Yearly Mean temperature anomaly from 1880 to 2015'])
xlabel(['Global '+season_name+' Seasonal Mean Temperature Anomaly (K)'])
ylabel('Global Yearly Mean Temperature Anomaly (K)')
end
%%
%
% Pretend that we are on April 1, 2016.
% Use the data from 1880 to December 2015 to calculate the best linear fit
% Use JFM of 2016 as a predictor for the 2016 yearly mean projection.
% Use observed 2016 yearly mean to compare your projection with reality
% We will use 2017 to repeat the exercise.
%

%
%don't remove!!! end has been removed so data no goes until 2022
%now using this section to select my 80 percent of data
seasonal_means_p = seasonal_means_eqpac(P);
yearly_mean_p = yearly_mean_mat(P);

%polyfit
p = polyfit(seasonal_means_p,yearly_mean_p,1);

%start with normal plot with line of best fit
%figure(3)
%plot(seasonal_means,yearly_mean_mat,'.')
%title('Global Jan-Feb-Mar Mean temperature anomaly as a function of the Global Yearly Mean temperature anomaly from 1880 to 2015')
%xlabel('Global Jan-Feb-Mar Mean Temperature Anomaly (K)')
%('Global Yearly Mean Temperature Anomaly (K)')
%hold on
%x1=linspace(-1,2,50);
y1= @(x) p(1)*x + p(2);
%plot(x1,y1(x1),'--')


%find stdev

difs = [];
for i = 1:length(seasonal_means_p)
    difs = [difs, (yearly_mean_p(i)-y1(seasonal_means_p(i)))^2];
end
standdev = sqrt(sum(difs)/length(seasonal_means_p));
stds(season,iterate) = standdev;
yup = @(x) p(1)*x + p(2)+2*standdev;    %upper bound on 95%
ybt = @(x) p(1)*x + p(2)-2*standdev;    %lower bound on 95%
if doplot
figure(5+season)
plot(seasonal_means_p,yearly_mean_p,'.b','DisplayName','Data')
title(['Global Yearly Mean Temperature Anomaly as a function of the '+season_name+' Seasonal Mean Temperature Anomaly from 1880 to 2015'])
xlabel(['Global '+season_name+' Seasonal Mean Temperature Anomaly (K)'])
ylabel('Global Yearly Mean Temperature Anomaly (K)')
hold on
x1=linspace(-1,2,50);
y1= @(x) p(1)*x + p(2);
plot(x1,y1(x1),'--k','DisplayName','Line of Best Fit')
hold on
plot(x1,yup(x1),'--k','DisplayName','95% Certainty Interval')
hold on
plot(x1,ybt(x1),'--k','DisplayName','95% Certainty Interval')
hold on
%plot predeictions
%plot(seasonal_means(end-6), y1(seasonal_means(end-6)),'r+','DisplayName','2016 Prediction')
hold on
%plot(seasonal_means(end-5), y1(seasonal_means(end-5)),'r+','DisplayName','2017 Prediction')
%plot reality
%plot(seasonal_means(end-6), yearly_mean_mat(end-6),'g^','DisplayName','2016 Outcome')
hold on
%plot(seasonal_means(end-5), yearly_mean_mat(end-5),'g^','DisplayName','2017 Outcome')

legend({'Data','Line of Best Fit','95% Certainty Interval','95% Certainty Interval'},'Location','southeastoutside')
end
end
end
%%
%process stds
mean_stds = squeeze(mean(stds,2));
std_stds = std(stds,0,2);
%%
%plot results
disp(["Iterations: "+num2str(iter_num)])
figure(1000)
y = mean_stds;
x = 1:4;
bar(y)
set(gca, 'XTickLabel',{"Winter (JFM)", "Spring (AMJ)", "Summer (JAS)", "Fall (OND)"})
title(["Standard Deviation Associated with Each Season's Preditictive Power. Iterations:" num2str(iter_num)])
ylabel('Standard Deviation (K)')

hold on

er = errorbar(x,y,std_stds,std_stds);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off