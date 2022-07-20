clc
clear
% % Code to run milfoil growth model from Miller et al. 2011 (Original Model: Herb et
% % al. 2006). This does not include Weevils. Diana White (Last modified:
% % June 2017)
% %----------------------------------------------------------

% % preallocation of functions
conv1=1/(100);% to go from 1/m to 1/cm
conv2=1/100^2;% to go from 1/m^2 to 1/cm^2
conv3=1/100^3;
  
  timestep=0.5;
  % % discretization
  dt=1; % time step is one day
  N=270;  %270; % run simulation over N days
  T1 = 1:dt:N; % discretization



  I0=zeros(1,numel(T1)); %Daily irradiance levels
  T=zeros(1,numel(T1));  %Daily temperature
  W=zeros(1,numel(T1));  %Daily milfoil growth
  lambda=zeros(1,numel(T1));%Daily loss due to respiration
  Id=zeros(1,numel(T1));
  Im=zeros(1,numel(T1));
  h=zeros(1,numel(T1));
  Id2=zeros(1,numel(T1));
  Im2=zeros(1,numel(T1));
  
  
  % % Model parameters taken from Herb 2006 and 
  U0=0.1;
  thetag=1.06;
 
  Tb= 10;   %origina;10;
  Kwt= 0.5;  %original: 0.5;
  km= .006;  % original: 0.006;
  Ps=50;
% d=0.5; % water depth in meters
 %d=0.4; % water depth in meters (Eric's model)
  d=2;
  k1=90;

  lambda0=0.01;
  thetar=1.072;
  delta=0.0;
  
  
  
% % original description  
%   DL=12;
%   for n=0:DL
%   I0(n)=I0p*sin(pi*n/DL)*sin(pi*n/DL);
%   end
%   for n=DL:24
%   I0(n)=0;
%   end

% % new description of light IO and temp T change over 360 days (used in
% % original Herb paper).


 
for n=1:numel(T1)
%I0(1,n) = 90;    
 I0(1,n)=50*sin(pi*(n+15)/360);
 T(1,n)=25*sin(pi*(n-20)/360)*sin(pi*(n-20)/360);
%T(1,n)=15;

lambda(1,n)=lambda0*thetar^(T(1,n)-Tb);
end
  W(1,1)=5*conv2;   %% conversion made here
  
  
%%
  for n=1:numel(T1) %days
%   h=0.5;

  h(1,n)=min(d,(W(1,n)/conv2)/Ps);

  %h(1,n)=min(d,(W(1,n)/conv2)/Ps);

 % h(1,n)=d;
      
  Im(1,n)= I0(1,n)*exp(+Kwt*d-Kwt*h(1,n));
  Id(1,n)= I0(1,n)*exp(-Kwt*d-(km/conv2*W(1,n))); %%%conversion made here
  
%   Im2(1,n)= Im(1,n)*conv2;
%   Id2(1,n)= Id(1,n)*conv2;
  
  % writing out the differential equation for watermilfoil growth using
  % simple Explicit Euler Scheme.
  
  W(1,n+1) = W(1,n) + dt*(U0*thetag^(T(1,n)-Tb)*W(1,n))./((Kwt*h(1,n)+(km/conv2*W(1,n)))).*log((k1+Im(1,n))./(k1+Id(1,n)))-dt*(lambda(1,n) + delta).*W(1,n);
  
  
  % % conversion made in km term
  
  end

  
 


 %% To plot solution of max days Nmax
%  figure(1)
%  hold on
%  plot(h(1,:)*100,'r')
 figure(1)
 hold on
  N2=1:numel(T1)+1;
  plot(N2,W(1,:),'r')
 title(['Biomass over a Season (' num2str(d*100) 'cm depth)'])
    xlabel('Day') 
    ylabel('Biomass (g/cm^2)') 
 set(gca,'FontSize',18)
    
figure (2)
hold on
h_cm=h*100;
  N2=1:numel(T1)+1;
  plot(h_cm,'r')
 title(['Plant Height over a Season (' num2str(d*100) 'cm depth)'])
    xlabel('Day') 
    ylabel('Height (cm)')
    axis([0 270 0 230])
 set(gca,'FontSize',18)



% set(gca,'FontSize',18)