%Este es un script para replicar la actividad neuronal de un experimento 
close all; clear all; clc;
stimulationTime = 50; %in ms
deltaT= 0.01;
t=0:deltaT:stimulationTime;
%t=30:deltaT:stimulationTime;
tiempo=length(t);
 
%% Specify the external current
changeTimes = [0]; % in ms
%currentLevels = 20;%[6.231]; %change this to see effect of different currents on nA (suggested values: 3,20,50,1000)
currentLevels = randn(1,tiempo)+2.231;%13;%[6.231];  %<<<<<----------------------------- BASAL 
%currentLevels = randn(1,tiempo)+30;%13;%[6.231]; 

I(1:numel(t)) = currentLevels;
%I(15000:22400)= randn(1,length (I(15000:22400)))+25;%50; %<<<<<<<<<<<<<<---------------VIBRISA



%%   Constant parameter
% stimulationTime = 50; %in ms
% deltaT= 0.01;
% t=0:deltaT:stimulationTime;
%All of these can be found in table 3 art: Hodking Huxley 1952
gk= [20:50];
gk= gk(randperm(length(gk)));
gbar_K=gk(1);%36;

gNa= [80:160];
gNa= gNa(randperm(length(gNa)));
gbar_Na=gNa(1);%120;

gbar_Ca= 3;

%gbar_syn = 6; %<<<<<<<<<<<------Revisar
A= 0.1;
gsyn= 0.30;
Wsyn = normrnd(-10,10);
gbar_syn = (A*Wsyn+1)*gsyn; %<<<<<<<<<<<------Revisar

g_L=0.3;

E_K=-88; 
E_Na = 60; 
E_Ca = 120;
E_L=-54.387; %%%%%%<-----   REVISAR 
E_syn= 0;
C=1;

%%   Set the initial states 
%rat
V=-70; %baseline voltage 
alpha_n = (0.01*(V+55))/(1-exp(-(V+55)/10)); %equation 12 K
beta_n = 0.125*exp(-(V+65)/80); %Equation 13 K

alpha_m = (0.1*(V+40)/(1-exp(-(V+40)/10)));
beta_m = 4*exp(-(V+65)/20); %Equation 21   Na

alpha_h = 0.07*exp(-(V+65)/20); %Equation 23   Na
beta_h = 1/(1+exp(-(V+35)/10)); %Equation 24  Na


n(1) = alpha_n/(alpha_n+beta_n); %equation 9
m(1) = alpha_m/(alpha_m+beta_m); %equation 18
h(1) = alpha_h/(alpha_h+beta_h); %equation 18
M(1) = 1/(1+exp(-(V(1)+57)/6.2));
H(1) = 1/(1+exp((V(1)+81)/4));
tM = 0.612+1/(exp(-(V(1)+132)/16.7)+exp((V(1)+16.8)/18.2));
tH = 28+exp(-(V(1)+22)/10.5);

S=100;

for  i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time 
    
    S1= randperm(S); %<- ruido 
    S2= randperm(S);%<- ruido 


%     S1 = randi([90 100],1,1);
%     S2 = randi([90 100],1,1);

%      S1 = randi([0 90],1,1); %cierre 
%      S2 = randi([0 70],1,1); %apertura



    alpha_n(i) = (0.01*(V(i)+55))/(1-exp(-(V(i)+55)/10)); %equation 12 K
    beta_n(i) = 0.125*exp(-(V(i)+65)/80); %Equation 13 K

    alpha_m(i) = (0.1*(V(i)+40)/(1-exp(-(V(i)+40)/10)));
    beta_m(i) = 4*exp(-(V(i)+65)/20); %Equation 21   Na

    alpha_h(i) = 0.07*exp(-(V(i)+65)/20); %Equation 23   Na
    beta_h(i) = 1/(1+exp(-(V(i)+35)/10)); %Equation 24  Na

     %   Calculate the currrents
    I_Na = (m(i)^3)*gbar_Na*h(i)*(V(i)-E_Na); %equations 3 and 14
    I_Na2(i+1) = (m(i)^3)*gbar_Na*h(i)*(V(i)-E_Na); 

    I_K = (n(i)^4)*gbar_K*(V(i)-E_K);%Equations 4 and 6
    I_K2(i+1) = (n(i)^4)*gbar_K*(V(i)-E_K);


     % With Ca
    tM = 0.612+1/(exp(-(V(i)+132)/16.7)+exp((V(1)+16.8)/18.2));
    if V(i)<-80
        tH = exp((V(i)+467)/66.6);
    else
        tH = 28+exp(-(V(1)+22)/10.5);
    end

    M(i+1) = (tM^-1)*(-M(i)+M(1));
    H(i+1) = (tH^-1)*(-H(i)+H(1));
    
    
    
    
    I_Ca = gbar_Ca*(M(i)^2)*H(i)*(V(i)-E_Ca);
    I_Ca2(i+1) = gbar_Ca*(M(i)^2)*H(i)*(V(i)-E_Ca);
    
    %with synapsis
    a = 0.33; %ms mM
    B = 0.2; %ms
    if i<= 300;
        rj = (1-exp(-(a*i)));
    else 
        rj=(1-exp(-(a*i)))*exp(-(B*(i-3)));

       
    end 
    rj2(i+1)=rj;
    I_syn = gbar_syn*rj*(V(i)-E_syn);
    I_syn2(i+1)= gbar_syn*rj*(V(i)-E_syn);

% normal current     
%     I_L = g_L*(V(i)-E_L); %Equations 5
%     I_ion = I(i)-I_K-I_Na-I_L;
%     I_ion2(i+1) = I(i)-I_K-I_Na-I_L;


% % With Ca current 
%     I_L = g_L*(V(i)-E_L); %Equations 5
%     I_ion = I(i)-I_K-I_Na-I_L-I_Ca;
%     I_ion2(i+1) = I(i)-I_K-I_Na-I_L-I_Ca;


% With synapsis  current 
    I_L = g_L*(V(i)-E_L); %Equations 5
    I_ion = I(i)-I_K-I_Na-I_L-I_Ca-I_syn;
    I_ion2(i+1) = I(i)-I_K-I_Na-I_L-I_Ca-I_syn;    
    
    
    %   Calculate the derivates usion Euler first order approximation,    States open-close 
     V(i+1) = V(i) + deltaT*I_ion/C;
    n(i+1) = n(i) + deltaT*((S2(1)/100)*alpha_n(i)*(1-n(i)) - (S1(1)/100)*beta_n(i)*n(i)); %Equation 7
    m(i+1) = m(i) + deltaT*((S2(1)/100)*alpha_m(i)*(1-m(i)) - (S1(1)/100)*beta_m(i)*m(i)); %Equation 15
    h(i+1) = h(i) + deltaT*((S2(1)/100)*alpha_h(i)*(1-h(i)) - (S1(1)/100)*beta_h(i)*h(i)); %Equation 16
end 


plot(t,V)
title('voltage')
figure
plot(t,rj2)
title('fraction of bound receptors described')
figure 
plot(t,I_syn2)
title('I_syn')
