c = @cmu.colors;
%creating a data set

%loading the data 
load('Pileup.mat')

%LacZ gene position on the chromosome 
Lac = 365000;

%The span of the chromosome starting at the DSB 
%to fit the model, for example 100 kb
Total = 100000;

% Select the region to the left or to the right of the break for the data
% set with 1 Chi in the Chi array 
data_left= wrev(WholePileup5322hits(Lac-Total: Lac)); 
data_right = WholePileup5322hits(Lac :Lac + Total); 
sequence_left = WholePileupSequence5322(Lac -Total: Lac ); 
sequence_right = WholePileupSequence5322(Lac: Lac  + Total); 

%%
%Find the positions of the Chi sites
%..on the orgin side
Chi = 'GCTGGTGG';
ChiPositions= wrev(Total - strfind(sequence_left', Chi));
%.. and on the terminus side
ChiRev = 'CCACCAGC';
ChiPositionsRev = strfind(sequence_right', ChiRev);

%%
%Binning of the raw data,
%the size of the bin is 
w=250;

%The binned data 
data_right_binned = data_bin(data_right,w);
data_left_binned = data_bin(data_left,w);

%Correct the positions of Chi sites for binning
ChiLacRight = round(ChiPositionsRev./w);
ChiLacLeft = round(ChiPositions./w);

%Saving the data sets
save('data_right_binned.mat', 'data_right_binned')
save('ChiLacRight.mat', 'ChiLacRight')
save('data_left_binned.mat', 'data_left_binned')
save('ChiLacLeft.mat', 'ChiLacLeft')

%%
%plotting the binned data on the right side of DSB 
figure, hold on
plot(data_left_binned ,'Color', c('navy blue'),'LineWidth',1.1);
plot(0 ,0,'diamond','MarkerSize',15,'MarkerFaceColor',	c( 'raspberry'),'MarkerEdgeColor',	'k')
plot(ChiLacLeft,0,'o','MarkerSize', 15,'MarkerFaceColor', c( 'amber'),'MarkerEdgeColor',	'k')
%and on the left side 
figure, hold on
plot(data_right_binned,'Color', c('navy blue'),'LineWidth',1.1);
plot(0 ,0,'diamond','MarkerSize',15,'MarkerFaceColor',	c( 'raspberry'),'MarkerEdgeColor',	'k')
plot(ChiLacRight,0,'o','MarkerSize', 15,'MarkerFaceColor',	c( 'amber'),'MarkerEdgeColor',	'k')

%%
% Subtract the background

background_region = 2000;
data = backgroundSubtraction(data_left_binned, background_region,w);

%% Grid Sampling - Simple Model
tic 
%Specify the grid for each of the three parameters
ps_range = [0.02 :0.001:0.03];
pchi_range = [0.0:0.1:1];
pm_range = [0.0:0.1:1];
  
[L,T,J,K, ps_max,pchi_max,pm_max] = ...
RecBCDgridSamplingSimple(ps_range,pchi_range,pm_range,data, ChiLacLeft)

toc
%%
%PLot the predition of the simple model together with the raw data and the smoothed data 

RecBCDplot([ps_max,pchi_max,pm_max], data, ChiLacLeft,w);
print -depsc -tiff -r1200 fitLeftLac

%%
%%Compute BIC score for the optimal model

[BIC_simple] = RecBCDdeterministic_BIC([ps_max,pchi_max,pm_max],data_left_binned,ChiLacLeft)

%% Parameter uncertainty - MCMC, simple model

theta_0 =[ps_max,pchi_max,pm_max];
epsilon_0 = [0.001,0.01,0.01];
[Theta, acceptance] = MC_RecBCDsimple(data, ChiLacLeft, 10000,theta_0,epsilon_0 )                                                                                                                                                                                               


figure, 
 subplot(1,3,1)
 legend('ps', 'pchi1','pm1' )
boxplot(Theta(1,:)/250,'outliersize',1)
legend('ps')
subplot(1,3,2)
boxplot(Theta(2,:),'outliersize',1)
legend('pchi')
subplot(1,3,3)
boxplot(Theta(3,:),'outliersize',1)

legend('pm' )

% print -depsc -tiff -r1200 5chiConfInt

%% Grid Sampling - Mixed Model

tic,

%Grid
ps_range = [0.02: 0.002: 0.03];
pchi_range1 = [0.2: 0.02: 0.4];
pm_range1 = [0.8: 0.02: 1];
pchi_range2 = [0.7: 0.02: 1];
pm_range2 =[0.4: 0.02: 0.6];
a_range = [0.5: 0.02: 0.6];


[log_likelihood,T,J,K,JJ,KK,N, ps_max,pchi_max1,pm_max1,pchi_max2,pm_max2,a_max]= ...
     RecBCDgridSamplingMixed(ps_range,pchi_range1,pm_range1,pchi_range2,pm_range2,a_range, data, ChiLacLeft);
%% 
%PLot the predition of the mixed model together with the raw data and the smoothed data 
RecBCDplotMixed([ps_max(1),pchi_max1(1),pm_max1(1),pchi_max2(1),pm_max2(1),a_max(1)], data, ChiLacLeft,w);

%%
%%Compute BIC score for the optimal model
[BIC_mixed] = RecBCDdeterministicMixed_BIC([ps_max,pchi_max1,pm_max1,pchi_max2,pm_max2,a_max],data_left_binned,ChiLacLeft)

%%
%Compare the BIC scores of the simple and mixed model 

delta_BIC = BIC_mixed - BIC_simple

% I fthe difference between the score is larger than 7 there is a strong
% preference towards mixed model 

%%
