function experimento4
close all 
%variables globales 
f1= figure;
f2= figure;
f3= figure;
f4= figure;
f5= figure;
f6= figure;
f7= figure;
f8= figure;
f9= figure;
f10= figure;





% Zona  con minutos TRATAMIENTO
N= input('Escribe los minutos de estimulacion:');
% High frequency current stimulation 
A = 0.5; %asintota inferior 
K = 2.5; %asintota superior 
B = 1.05; % tasa de crecimiento 
vi = 1;
Q = 70; %potencial de membrana en reposo 
C= 1; 
E= 20; %input('Escribe las neuronas control  del experimento:');  % mismo valor que while 
Y = A+((K-A)/((C+Q*exp(-B*N))^(1/vi)));



contador=1;

stimulationTime = 50; %in ms
deltaT= 0.01;
t=0:deltaT:stimulationTime;
t3= 0:deltaT:50;
%t=30:deltaT:stimulationTime;
tiempo=length(t);
t2= 0.01:deltaT:50;%length(t)-1;
 
%% Specify the external current
changeTimes = [0]; % in ms
%currentLevels = 20;%[6.231]; %change this to see effect of different currents on nA (suggested values: 3,20,50,1000)
currentLevels =9.47979;% 0.2*(randn(1,tiempo));%+2.231;%13;%[6.231];  %<<<<<----------------------------- BASAL 
%currentLevels = randn(1,tiempo)+30;%13;%[6.231]; 




I(1:numel(t)) = currentLevels;
%I(1:400)=-100;
%I(15000:22400)= randn(1,length (I(15000:22400)))+25;%50; %<<<<<<<<<<<<<<---------------VIBRISA
% figure(f1) 
% plot(t,I)
% title('corriente sobre E')
% xlabel('ms')
% ylabel('nA')
% %xlim([50, 100])


%% iniciamos union de NEURONAS 

while contador<=20 %%<<<<<<<<---------------nivel 
    %% NEURONA EXCITADORA <<<<<<<<<<<<<<<<<<<<<<<<<---------------------------
neuronE= neuronaE(I);
neuronE>30;
p=neuronE;
[p,z]=findpeaks(neuronE,t);
p2=z(p>20); %<<<<<<------------este es RASTER COPIAR Y PEGAR
x= [1:20];
x2=x(randperm(length(x)));
p2=p2+x2(1,1);
% p3=p2(p2>74);
% p3=p3(p3<300);


figure(f2)
hold on
plot(p2,contador,'.');
%figure
% plot(t,neuronE,"m")
 title('Raster E')
 xlim([0, 50])
 ylim([0, 20])
% ylim([-100, 50])

% liberacion de vesiculas 
p3 = Poisson(t);
% figure 
% plot(p3)
% title('probabilidad de liberacion')




% corriente exctiadora 
I_E= 10*diff(neuronE); %hacemos esto porque la diferencia es pequeña y necesitamos emular la corriente de la neurona
I_E4= I_E;%(5000:10000);
figure(f3)

I_E2=abs(I_E4);
I_E3= ((I_E2+I_E4)/2).*(1+p3(:,:));%1.4; %<<<<<<<----- PPE
%I_E3=((I_E-I_E2)/2).*(1+p3(:,:));%2.0;  %<<<<<<<----- PPI
% hold on 
% plot(t2,I_E3)
%xlim([3.62, 35.89])
% title('corriente sobre FS')
% xlabel('ms')
% ylabel('nA')

if contador>1
    I_E3=I_E3(1:length(I_SST3))+I_FS3+I_SST3;
else 
    I_E3= I_E3;
end 



%% NEURONA FS <<<<<<<<<<<<<<<<<<<<<<<<<---------------------------
neuronFS= neuronaFS(I_E3);
neuronFS>30;
p=neuronFS;
[p,z]=findpeaks(neuronFS,t);
p2=z(p>20); %<<<<<<------------este es RASTER COPIAR Y PEGAR
x= [1:20];
x2=x(randperm(length(x)));
p2=p2+x2(1,1);
% p3=p2(p2>74);
% p3=p3(p3<300);





figure(f4)
hold on
plot(p2,contador,'.');
%plot(t3,neuronFS,"r")
title('Raster FS')
xlim([0, 50])
ylim([0, 20])
%ylim([-100, 50])

                % corriente inhibitoria de fs 
                I_FS = 10*diff(neuronFS); %hacemos esto porque la diferencia es pequeña y necesitamos emular la corriente de la neurona
%                 figure(f5)
%                 hold on
                %plot(t2,I_FS)
                I_FS2= abs(I_FS);
                %plot(t2,I_FS2)
                I_FS3= ((I_FS-I_FS2)/2).*(100+p3(:,:));%2.0;  %<<<<<<<----- PPI
%                 plot(t2,I_FS3)
%                 title('corriente aplicada sobre E y FS')
%                 xlabel('ms')
%                 ylabel('nA')



%% NEURONA SST <<<<<<<<<<<<<<<<<<<<<<<<<---------------------------
 %I_ISST= I_E3;
figure(f6)
plot(t2,I_E3) 
title('corriente sobre SST')

neuronSST= neuronaSST(I_E3);
neuronSST>30;
p=neuronSST;
[p,z]=findpeaks(neuronSST,t);
p2=z(p>20); %<<<<<<------------este es RASTER COPIAR Y PEGAR
x= [1:20];
x2=x(randperm(length(x)));
p2=p2+x2(1,1);




figure(f7)
hold on
plot(p2,contador,'.');
%plot(t3,neuronSST,"g")
title('Raster SST')
xlim([0, 50])
ylim([0, 20])
%ylim([-100, 50])


                % corriente inhibitoria de SST 
                I_SST = 10*diff(neuronSST); %hacemos esto porque la diferencia es pequeña y necesitamos emular la corriente de la neurona
%                 figure (f8)
                hold on
                %plot(t2,I_FS)
                I_SST2= abs(I_SST);
                %plot(t2,I_FS2)
                I_SST3= ((I_SST-I_SST2)/2).*(100+p3(:,:));%2.0;  %<<<<<<<----- PPI
%                 plot(t2,I_SST3)
%                 title('corriente aplicada')
%                 xlabel('ms')
%                 ylabel('nA')


%if contador >2
 %% NEURONA VIP <<<<<<<<<<<<<<<<<<<<<<<<<---------------------------   
I2= I_FS3+I_SST3;
figure(f9)
plot(t2,I2) 
title('corriente sobre VIP')                

neuronVIP= neuronaVIP(I2);
neuronVIP>30;
p=neuronVIP;
[p,z]=findpeaks(neuronVIP,t);
p2=z(p>20); %<<<<<<------------este es RASTER COPIAR Y PEGAR
x= [1:20];
x2=x(randperm(length(x)));
p2=p2+x2(1,1);




figure(f10)
hold on
plot(p2,contador,'.')
%plot(t3,neuronVIP,"b")
title('Raster VIP')
xlim([0, 50])
ylim([0, 20])
%ylim([-100, 50])

                
                
                
I= I_FS3+I_SST3+I(1:length(I_SST3));
%I_E3= 
contador= contador+1;
end 
'terminado :D'

end 

function p3 = Poisson(t)
%observaciones 
T= 0.05;
t=[0:0.001:1.4];
N= 50;
n= N*t./T.*exp(-t./T);
% figure
% plot(t,n)
% 
% ylabel('eventos (liberacion de vesiculas)')
% xlabel ('time s')
% xlim([0 0.05])

%% ajuste de liberacion de vesiculas por poisson 

X= n;%[0:100];
for i=1:numel(X)-1
x=round(X(i));%70; % numeros de veces que ocurre un evento durante un intervalo definido ESTE EJEMPLO t 1.4s
mu= mean(n);% probabilidad de ocurrencia es la misma para cualesquiera de los intervalos de igual longitud
P= (exp(-mu)*mu^x)/factorial(round(x)); %funcion de probabilidad de poisson 
P2(i+1)=P;
P4(i+1)=sum(P2);  % funcion de densidad de probabilidad  de poisson 
end 

%opcional 
p3=smooth(P2);
% figure
% plot(t,p3)
% ylabel('probabilidad de eventos')
% xlabel ('time s')
% xlim([0 0.05])

% figure
% plot(t,P4)
% ylabel('probabilidad acumulada')
% xlabel ('time s')
% xlim([0 0.05])

end 

function VE = neuronaE(I)
%%   Constant parameter
% Zona  con minutos TRATAMIENTO
N= 10;%input('Escribe los minutos de estimulacion:');
% High frequency current stimulation 
A = 0.5; %asintota inferior 
K = 2.5; %asintota superior 
B = 1.05; % tasa de crecimiento 
vi = 1;
Q = 70; %potencial de membrana en reposo 
C= 1; 
E= 20; %input('Escribe las neuronas control  del experimento:');  % mismo valor que while 
Y = A+((K-A)/((C+Q*exp(-B*N))^(1/vi)));








stimulationTime = 50; %in ms
deltaT= 0.01;
t=0:deltaT:stimulationTime;
%All of these can be found in table 3 art: Hodking Huxley 1952
gk= [20:50];
gk= gk(randperm(length(gk)));
%gbar_K=gk(1);%36;
gbar_K=36;

gNa= [80:160];
gNa= gNa(randperm(length(gNa)));
%gbar_Na=gNa(1);%120;
gbar_Na=120;

gbar_Ca= 3;

gbar_syn = 6; %<<<<<<<<<<<------Revisar

g_L=0.3;

E_K=-88; 
E_Na = 60; 
E_Ca = 120;
E_L=-54.387; %%%%%%<-----   REVISAR 
E_syn= 0;
C=1;

%%   Set the initial states 
%rat
VE=-70; %baseline voltage 
alpha_n = (0.01*(VE+55))/(1-exp(-(VE+55)/10)); %equation 12 K
beta_n = 0.125*exp(-(VE+65)/80); %Equation 13 K

alpha_m = Y*(0.1*(VE+40)/(1-exp(-(VE+40)/10)));
beta_m = Y*(4*exp(-(VE+65)/20)); %Equation 21   Na

alpha_h = Y*(0.07*exp(-(VE+65)/20)); %Equation 23   Na
beta_h = 1/(1+exp(-(VE+35)/10)); %Equation 24  Na


n(1) = alpha_n/(alpha_n+beta_n); %equation 9
m(1) = alpha_m/(alpha_m+beta_m); %equation 18
h(1) = alpha_h/(alpha_h+beta_h); %equation 18
M(1) = 1/(1+exp(-(VE(1)+57)/6.2));
H(1) = 1/(1+exp((VE(1)+81)/4));
tM = 0.612+1/(exp(-(VE(1)+132)/16.7)+exp((VE(1)+16.8)/18.2));
tH = 28+exp(-(VE(1)+22)/10.5);

S=100;

for  i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time 
    
    S1=100; %randperm(S); %<- ruido 
    S2=100;% randperm(S);%<- ruido 


%     S1 = randi([90 100],1,1);
%     S2 = randi([90 100],1,1);

%      S1 = randi([0 90],1,1); %cierre 
%      S2 = randi([0 70],1,1); %apertura



    alpha_n(i) = (0.01*(VE(i)+55))/(1-exp(-(VE(i)+55)/10)); %equation 12 K
    beta_n(i) = 0.125*exp(-(VE(i)+65)/80); %Equation 13 K

    alpha_m(i) = (0.1*(VE(i)+40)/(1-exp(-(VE(i)+40)/10)));
    beta_m(i) = 4*exp(-(VE(i)+65)/20); %Equation 21   Na

    alpha_h(i) = 0.07*exp(-(VE(i)+65)/20); %Equation 23   Na
    beta_h(i) = 1/(1+exp(-(VE(i)+35)/10)); %Equation 24  Na

     %   Calculate the currrents
    I_Na = (m(i)^3)*gbar_Na*h(i)*(VE(i)-E_Na); %equations 3 and 14
    I_Na2(i+1) = (m(i)^3)*gbar_Na*h(i)*(VE(i)-E_Na); 

    I_K = (n(i)^4)*gbar_K*(VE(i)-E_K);%Equations 4 and 6
    I_K2(i+1) = (n(i)^4)*gbar_K*(VE(i)-E_K);


     % With Ca
    tM = 0.612+1/(exp(-(VE(i)+132)/16.7)+exp((VE(1)+16.8)/18.2));
    if VE(i)<-80
        tH = exp((VE(i)+467)/66.6);
    else
        tH = 28+exp(-(VE(1)+22)/10.5);
    end

    M(i+1) = (tM^-1)*(-M(i)+M(1));
    H(i+1) = (tH^-1)*(-H(i)+H(1));
    
    
    
    
    I_Ca = gbar_Ca*(M(i)^2)*H(i)*(VE(i)-E_Ca);
    I_Ca2(i+1) = gbar_Ca*(M(i)^2)*H(i)*(VE(i)-E_Ca);
    
    %with synapsis
    a = 0.33; %ms mM
    B = 0.2; %ms
    if i<= 300;
        rj = (1-exp(-(a*i)));
    else 
        rj=(1-exp(-(a*i)))*exp(-(B*(i-3)));


    end 
    I_syn = gbar_syn*rj*(VE(i)-E_syn);
    I_syn2(i+1)= I_syn;
% normal current     
%     I_L = g_L*(V(i)-E_L); %Equations 5
%     I_ion = I(i)-I_K-I_Na-I_L;
%     I_ion2(i+1) = I(i)-I_K-I_Na-I_L;


% % With Ca current 
%     I_L = g_L*(V(i)-E_L); %Equations 5
%     I_ion = I(i)-I_K-I_Na-I_L-I_Ca;
%     I_ion2(i+1) = I(i)-I_K-I_Na-I_L-I_Ca;


% With synapsis  current 
    I_L = g_L*(VE(i)-E_L); %Equations 5
    I_ion = I(i)-I_K-I_Na-I_L-I_Ca-I_syn;
    I_ion2(i+1) = I(i)-I_K-I_Na-I_L-I_Ca-I_syn;    
    
    
    %   Calculate the derivates usion Euler first order approximation,    States open-close 
     VE(i+1) = VE(i) + deltaT*I_ion/C;
    n(i+1) = n(i) + deltaT*((S2(1)/100)*alpha_n(i)*(1-n(i)) - (S1(1)/100)*beta_n(i)*n(i)); %Equation 7
    m(i+1) = m(i) + deltaT*((S2(1)/100)*alpha_m(i)*(1-m(i)) - (S1(1)/100)*beta_m(i)*m(i)); %Equation 15
    h(i+1) = h(i) + deltaT*((S2(1)/100)*alpha_h(i)*(1-h(i)) - (S1(1)/100)*beta_h(i)*h(i)); %Equation 16

% plot(I_syn)
   
 end 

% figure
% plot(I_syn2)
% title('synapse E')
end 

function VSST = neuronaSST(I_E3)
%%   Constant parameter
% Zona  con minutos TRATAMIENTO
N= 10;%input('Escribe los minutos de estimulacion:');
% High frequency current stimulation 
A = 0.5; %asintota inferior 
K = 2.5; %asintota superior 
B = 1.05; % tasa de crecimiento 
vi = 1;
Q = 70; %potencial de membrana en reposo 
C= 1; 
E= 20; %input('Escribe las neuronas control  del experimento:');  % mismo valor que while 
Y = A+((K-A)/((C+Q*exp(-B*N))^(1/vi)));




stimulationTime = 50; %in ms
deltaT= 0.01;
t=0:deltaT:stimulationTime;
%All of these can be found in table 3 art: Hodking Huxley 1952
gk= [20:50];
gk= gk(randperm(length(gk)));
%gbar_K=gk(1);%36;
gbar_K=36;

gNa= [80:160];
gNa= gNa(randperm(length(gNa)));
%gbar_Na=gNa(1);%120;
gbar_Na=120;

gbar_Ca= 3;

gbar_syn = 6; %<<<<<<<<<<<------Revisar

g_L=0.3;

E_K=-88; 
E_Na = 60; 
E_Ca = 120;
E_L=-54.387; %%%%%%<-----   REVISAR 
E_syn= 0;
C=1;

%%   Set the initial states 
%rat
VSST=-70; %baseline voltage 
alpha_n = (0.01*(VSST+55))/(1-exp(-(VSST+55)/10)); %equation 12 K
beta_n = 0.125*exp(-(VSST+65)/80); %Equation 13 K

alpha_m = Y*(0.1*(VSST+40)/(1-exp(-(VSST+40)/10)));
beta_m =Y*(4*exp(-(VSST+65)/20)); %Equation 21   Na

alpha_h = Y*(0.07*exp(-(VSST+65)/20)); %Equation 23   Na
beta_h = 1/(1+exp(-(VSST+35)/10)); %Equation 24  Na


n(1) = alpha_n/(alpha_n+beta_n); %equation 9
m(1) = alpha_m/(alpha_m+beta_m); %equation 18
h(1) = alpha_h/(alpha_h+beta_h); %equation 18
M(1) = 1/(1+exp(-(VSST(1)+57)/6.2));
H(1) = 1/(1+exp((VSST(1)+81)/4));
tM = 0.612+1/(exp(-(VSST(1)+132)/16.7)+exp((VSST(1)+16.8)/18.2));
tH = 28+exp(-(VSST(1)+22)/10.5);

S=100;

for  i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time 
    
    S1= 100;%randperm(S); %<- ruido 
    S2= 100; %randperm(S);%<- ruido 


%     S1 = randi([90 100],1,1);
%     S2 = randi([90 100],1,1);

%      S1 = randi([0 90],1,1); %cierre 
%      S2 = randi([0 70],1,1); %apertura



    alpha_n(i) = (0.01*(VSST(i)+55))/(1-exp(-(VSST(i)+55)/10)); %equation 12 K
    beta_n(i) = 0.125*exp(-(VSST(i)+65)/80); %Equation 13 K

    alpha_m(i) = (0.1*(VSST(i)+40)/(1-exp(-(VSST(i)+40)/10)));
    beta_m(i) = 4*exp(-(VSST(i)+65)/20); %Equation 21   Na

    alpha_h(i) = 0.07*exp(-(VSST(i)+65)/20); %Equation 23   Na
    beta_h(i) = 1/(1+exp(-(VSST(i)+35)/10)); %Equation 24  Na

     %   Calculate the currrents
    I_Na = (m(i)^3)*gbar_Na*h(i)*(VSST(i)-E_Na); %equations 3 and 14
    I_Na2(i+1) = (m(i)^3)*gbar_Na*h(i)*(VSST(i)-E_Na); 

    I_K = (n(i)^4)*gbar_K*(VSST(i)-E_K);%Equations 4 and 6
    I_K2(i+1) = (n(i)^4)*gbar_K*(VSST(i)-E_K);


     % With Ca
    tM = 0.612+1/(exp(-(VSST(i)+132)/16.7)+exp((VSST(1)+16.8)/18.2));
    if VSST(i)<-80
        tH = exp((VSST(i)+467)/66.6);
    else
        tH = 28+exp(-(VSST(1)+22)/10.5);
    end

    M(i+1) = (tM^-1)*(-M(i)+M(1));
    H(i+1) = (tH^-1)*(-H(i)+H(1));
    
    
    
    
    I_Ca = gbar_Ca*(M(i)^2)*H(i)*(VSST(i)-E_Ca);
    I_Ca2(i+1) = gbar_Ca*(M(i)^2)*H(i)*(VSST(i)-E_Ca);
    
    %with synapsis
    a = 0.33; %ms mM
    B = 0.2; %ms
    if i<= 300;
        rj = (1-exp(-(a*i)));
    else 
        rj=(1-exp(-(a*i)))*exp(-(B*(i-3)));


    end 
    I_syn = gbar_syn*rj*(VSST(i)-E_syn);
    I_syn2(i+1)= I_syn;
% normal current     
%     I_L = g_L*(V(i)-E_L); %Equations 5
%     I_ion = I(i)-I_K-I_Na-I_L;
%     I_ion2(i+1) = I(i)-I_K-I_Na-I_L;


% % With Ca current 
%     I_L = g_L*(V(i)-E_L); %Equations 5
%     I_ion = I(i)-I_K-I_Na-I_L-I_Ca;
%     I_ion2(i+1) = I(i)-I_K-I_Na-I_L-I_Ca;


% With synapsis  current 
    I_L = g_L*(VSST(i)-E_L); %Equations 5
    I_ion = I_E3(i)-I_K-I_Na-I_L-I_Ca-I_syn;
    I_ion2(i+1) = I_E3(i)-I_K-I_Na-I_L-I_Ca-I_syn;    
    
    
    %   Calculate the derivates usion Euler first order approximation,    States open-close 
     VSST(i+1) = VSST(i) + deltaT*I_ion/C;
    n(i+1) = n(i) + deltaT*((S2(1)/100)*alpha_n(i)*(1-n(i)) - (S1(1)/100)*beta_n(i)*n(i)); %Equation 7
    m(i+1) = m(i) + deltaT*((S2(1)/100)*alpha_m(i)*(1-m(i)) - (S1(1)/100)*beta_m(i)*m(i)); %Equation 15
    h(i+1) = h(i) + deltaT*((S2(1)/100)*alpha_h(i)*(1-h(i)) - (S1(1)/100)*beta_h(i)*h(i)); %Equation 16

% plot(I_syn)
   
 end 

% figure
% plot(I_syn2)
% title('synapse ssp')
end 

function VFS = neuronaFS(I_E3)
%%   Constant parameter
% Zona  con minutos TRATAMIENTO
N=10;% input('Escribe los minutos de estimulacion:');
% High frequency current stimulation 
A = 0.5; %asintota inferior 
K = 2.5; %asintota superior 
B = 1.05; % tasa de crecimiento 
vi = 1;
Q = 70; %potencial de membrana en reposo 
C= 1; 
E= 20; %input('Escribe las neuronas control  del experimento:');  % mismo valor que while 
Y = A+((K-A)/((C+Q*exp(-B*N))^(1/vi)));











stimulationTime = 50; %in ms
deltaT= 0.01;
t=0:deltaT:stimulationTime;
%All of these can be found in table 3 art: Hodking Huxley 1952
gk= [20:50];
gk= gk(randperm(length(gk)));
%gbar_K=gk(1);%36;
gbar_K=10;

gNa= [80:160];
gNa= gNa(randperm(length(gNa)));
%gbar_Na=gNa(1);%120;
gbar_Na=90;

gbar_Ca= 3;

gbar_syn = 6; %<<<<<<<<<<<------Revisar

g_L=0.3;

E_K=-88; 
E_Na = 60; 
E_Ca = 120;
E_L=-54.387; %%%%%%<-----   REVISAR 
E_syn= 0;
C=1;

%%   Set the initial states 
%rat
VFS=-70; %baseline voltage 
alpha_n = (0.01*(VFS+55))/(1-exp(-(VFS+55)/10)); %equation 12 K
beta_n = 0.125*exp(-(VFS+65)/80); %Equation 13 K

alpha_m = Y*(0.1*(VFS+40)/(1-exp(-(VFS+40)/10)));
beta_m = Y*(4*exp(-(VFS+65)/20)); %Equation 21   Na

alpha_h = Y*(0.07*exp(-(VFS+65)/20)); %Equation 23   Na
beta_h = 1/(1+exp(-(VFS+35)/10)); %Equation 24  Na


n(1) = alpha_n/(alpha_n+beta_n); %equation 9
m(1) = alpha_m/(alpha_m+beta_m); %equation 18
h(1) = alpha_h/(alpha_h+beta_h); %equation 18
M(1) = 1/(1+exp(-(VFS(1)+57)/6.2));
H(1) = 1/(1+exp((VFS(1)+81)/4));
tM = 0.612+1/(exp(-(VFS(1)+132)/16.7)+exp((VFS(1)+16.8)/18.2));
tH = 28+exp(-(VFS(1)+22)/10.5);

S=100;

for  i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time 
    
    S1= 100;%randperm(S); %<- ruido 
    S2= 100; %randperm(S);%<- ruido 


%     S1 = randi([90 100],1,1);
%     S2 = randi([90 100],1,1);

%      S1 = randi([0 90],1,1); %cierre 
%      S2 = randi([0 70],1,1); %apertura



    alpha_n(i) = (0.01*(VFS(i)+55))/(1-exp(-(VFS(i)+55)/10)); %equation 12 K
    beta_n(i) = 0.125*exp(-(VFS(i)+65)/80); %Equation 13 K

    alpha_m(i) = (0.1*(VFS(i)+40)/(1-exp(-(VFS(i)+40)/10)));
    beta_m(i) = 4*exp(-(VFS(i)+65)/20); %Equation 21   Na

    alpha_h(i) = 0.07*exp(-(VFS(i)+65)/20); %Equation 23   Na
    beta_h(i) = 1/(1+exp(-(VFS(i)+35)/10)); %Equation 24  Na

     %   Calculate the currrents
    I_Na = (m(i)^3)*gbar_Na*h(i)*(VFS(i)-E_Na); %equations 3 and 14
    I_Na2(i+1) = (m(i)^3)*gbar_Na*h(i)*(VFS(i)-E_Na); 

    I_K = (n(i)^4)*gbar_K*(VFS(i)-E_K);%Equations 4 and 6
    I_K2(i+1) = (n(i)^4)*gbar_K*(VFS(i)-E_K);


     % With Ca
    tM = 0.612+1/(exp(-(VFS(i)+132)/16.7)+exp((VFS(1)+16.8)/18.2));
    if VFS(i)<-80
        tH = exp((VFS(i)+467)/66.6);
    else
        tH = 28+exp(-(VFS(1)+22)/10.5);
    end

    M(i+1) = (tM^-1)*(-M(i)+M(1));
    H(i+1) = (tH^-1)*(-H(i)+H(1));
    
    
    
    
    I_Ca = gbar_Ca*(M(i)^2)*H(i)*(VFS(i)-E_Ca);
    I_Ca2(i+1) = gbar_Ca*(M(i)^2)*H(i)*(VFS(i)-E_Ca);
    
    %with synapsis
    a = 0.33; %ms mM
    B = 0.2; %ms
    if i<= 300;
        rj = (1-exp(-(a*i)));
    else 
        rj=(1-exp(-(a*i)))*exp(-(B*(i-3)));


    end 
    I_syn = gbar_syn*rj*(VFS(i)-E_syn);
    I_syn2(i+1)= I_syn;
% normal current     
%     I_L = g_L*(V(i)-E_L); %Equations 5
%     I_ion = I(i)-I_K-I_Na-I_L;
%     I_ion2(i+1) = I(i)-I_K-I_Na-I_L;


% % With Ca current 
%     I_L = g_L*(V(i)-E_L); %Equations 5
%     I_ion = I(i)-I_K-I_Na-I_L-I_Ca;
%     I_ion2(i+1) = I(i)-I_K-I_Na-I_L-I_Ca;


% With synapsis  current 
    I_L = g_L*(VFS(i)-E_L); %Equations 5
    I_ion = I_E3(i)-I_K-I_Na-I_L-I_Ca-I_syn;
    I_ion2(i+1) = I_E3(i)-I_K-I_Na-I_L-I_Ca-I_syn;    
    
    
    %   Calculate the derivates usion Euler first order approximation,    States open-close 
     VFS(i+1) = VFS(i) + deltaT*I_ion/C;
    n(i+1) = n(i) + deltaT*((S2(1)/100)*alpha_n(i)*(1-n(i)) - (S1(1)/100)*beta_n(i)*n(i)); %Equation 7
    m(i+1) = m(i) + deltaT*((S2(1)/100)*alpha_m(i)*(1-m(i)) - (S1(1)/100)*beta_m(i)*m(i)); %Equation 15
    h(i+1) = h(i) + deltaT*((S2(1)/100)*alpha_h(i)*(1-h(i)) - (S1(1)/100)*beta_h(i)*h(i)); %Equation 16

% plot(I_syn)
   
 end 

% figure
% plot(I_syn2)
% title('synapse ssp')
end 

function VVIP = neuronaVIP(I2)
%%   Constant parameter
% Zona  con minutos TRATAMIENTO
N=10;% input('Escribe los minutos de estimulacion:');
% High frequency current stimulation 
A = 0.5; %asintota inferior 
K = 2.5; %asintota superior 
B = 1.05; % tasa de crecimiento 
vi = 1;
Q = 70; %potencial de membrana en reposo 
C= 1; 
E= 20; %input('Escribe las neuronas control  del experimento:');  % mismo valor que while 
Y = A+((K-A)/((C+Q*exp(-B*N))^(1/vi)));












stimulationTime = 50; %in ms
deltaT= 0.01;
t=0:deltaT:stimulationTime;
%All of these can be found in table 3 art: Hodking Huxley 1952
gk= [20:50];
gk= gk(randperm(length(gk)));
%gbar_K=gk(1);%36;
gbar_K=36;

gNa= [80:160];
gNa= gNa(randperm(length(gNa)));
%gbar_Na=gNa(1);%120;
gbar_Na=120;

gbar_Ca= 3;

gbar_syn = 6; %<<<<<<<<<<<------Revisar

g_L=0.3;

E_K=-88; 
E_Na = 60; 
E_Ca = 120;
E_L=-54.387; %%%%%%<-----   REVISAR 
E_syn= 0;
C=1;

%%   Set the initial states 
%rat
VVIP=-70; %baseline voltage 
alpha_n = (0.01*(VVIP+55))/(1-exp(-(VVIP+55)/10)); %equation 12 K
beta_n = 0.125*exp(-(VVIP+65)/80); %Equation 13 K

alpha_m = Y*(0.1*(VVIP+40)/(1-exp(-(VVIP+40)/10)));
beta_m = Y*(4*exp(-(VVIP+65)/20)); %Equation 21   Na

alpha_h = 0.07*exp(-(VVIP+65)/20); %Equation 23   Na
beta_h = 1/(1+exp(-(VVIP+35)/10)); %Equation 24  Na


n(1) = alpha_n/(alpha_n+beta_n); %equation 9
m(1) = alpha_m/(alpha_m+beta_m); %equation 18
h(1) = alpha_h/(alpha_h+beta_h); %equation 18
M(1) = 1/(1+exp(-(VVIP(1)+57)/6.2));
H(1) = 1/(1+exp((VVIP(1)+81)/4));
tM = 0.612+1/(exp(-(VVIP(1)+132)/16.7)+exp((VVIP(1)+16.8)/18.2));
tH = 28+exp(-(VVIP(1)+22)/10.5);

S=100;

for  i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time 
    
    S1= 100;%randperm(S); %<- ruido 
    S2= 100; %randperm(S);%<- ruido 


%     S1 = randi([90 100],1,1);
%     S2 = randi([90 100],1,1);

%      S1 = randi([0 90],1,1); %cierre 
%      S2 = randi([0 70],1,1); %apertura



    alpha_n(i) = (0.01*(VVIP(i)+55))/(1-exp(-(VVIP(i)+55)/10)); %equation 12 K
    beta_n(i) = 0.125*exp(-(VVIP(i)+65)/80); %Equation 13 K

    alpha_m(i) = (0.1*(VVIP(i)+40)/(1-exp(-(VVIP(i)+40)/10)));
    beta_m(i) = 4*exp(-(VVIP(i)+65)/20); %Equation 21   Na

    alpha_h(i) = 0.07*exp(-(VVIP(i)+65)/20); %Equation 23   Na
    beta_h(i) = 1/(1+exp(-(VVIP(i)+35)/10)); %Equation 24  Na

     %   Calculate the currrents
    I_Na = (m(i)^3)*gbar_Na*h(i)*(VVIP(i)-E_Na); %equations 3 and 14
    I_Na2(i+1) = (m(i)^3)*gbar_Na*h(i)*(VVIP(i)-E_Na); 

    I_K = (n(i)^4)*gbar_K*(VVIP(i)-E_K);%Equations 4 and 6
    I_K2(i+1) = (n(i)^4)*gbar_K*(VVIP(i)-E_K);


     % With Ca
    tM = 0.612+1/(exp(-(VVIP(i)+132)/16.7)+exp((VVIP(1)+16.8)/18.2));
    if VVIP(i)<-80
        tH = exp((VVIP(i)+467)/66.6);
    else
        tH = 28+exp(-(VVIP(1)+22)/10.5);
    end

    M(i+1) = (tM^-1)*(-M(i)+M(1));
    H(i+1) = (tH^-1)*(-H(i)+H(1));
    
    
    
    
    I_Ca = gbar_Ca*(M(i)^2)*H(i)*(VVIP(i)-E_Ca);
    I_Ca2(i+1) = gbar_Ca*(M(i)^2)*H(i)*(VVIP(i)-E_Ca);
    
    %with synapsis
    a = 0.33; %ms mM
    B = 0.2; %ms
    if i<= 300;
        rj = (1-exp(-(a*i)));
    else 
        rj=(1-exp(-(a*i)))*exp(-(B*(i-3)));


    end 
    I_syn = gbar_syn*rj*(VVIP(i)-E_syn);
    I_syn2(i+1)= I_syn;
% normal current     
%     I_L = g_L*(V(i)-E_L); %Equations 5
%     I_ion = I(i)-I_K-I_Na-I_L;
%     I_ion2(i+1) = I(i)-I_K-I_Na-I_L;


% % With Ca current 
%     I_L = g_L*(V(i)-E_L); %Equations 5
%     I_ion = I(i)-I_K-I_Na-I_L-I_Ca;
%     I_ion2(i+1) = I(i)-I_K-I_Na-I_L-I_Ca;


% With synapsis  current 
    I_L = g_L*(VVIP(i)-E_L); %Equations 5
    I_ion = I2(i)-I_K-I_Na-I_L-I_Ca-I_syn;
    I_ion2(i+1) = I2(i)-I_K-I_Na-I_L-I_Ca-I_syn;    
    
    
    %   Calculate the derivates usion Euler first order approximation,    States open-close 
     VVIP(i+1) = VVIP(i) + deltaT*I_ion/C;
    n(i+1) = n(i) + deltaT*((S2(1)/100)*alpha_n(i)*(1-n(i)) - (S1(1)/100)*beta_n(i)*n(i)); %Equation 7
    m(i+1) = m(i) + deltaT*((S2(1)/100)*alpha_m(i)*(1-m(i)) - (S1(1)/100)*beta_m(i)*m(i)); %Equation 15
    h(i+1) = h(i) + deltaT*((S2(1)/100)*alpha_h(i)*(1-h(i)) - (S1(1)/100)*beta_h(i)*h(i)); %Equation 16

% plot(I_syn)
   
 end 

% figure
% plot(I_syn2)
% title('synapse ssp')
end 

