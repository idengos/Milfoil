close all
clear
clc
% % Code to run milfoil growth model from Miller et al. 2011 (Original Model: Herb et
% % al. 2006). Extension includes diffusion in space in 1D (assume straight "river"). Diana White & Isabel Dengos (Last modified:
% % June 2020).

% % Update on may 30th 2020. 2D diffusion code which includes a
% % "break"/stop of the program if stability condition not met.

% % here we include reaction (plant growth) and diffusion of stolins (EWM
% % root system
% %----------------------------------------------------------


% % Stability condition

% % For stability of the 2D diffusion equation, we must have D <1/4*(dx^2/dt)



% % preallocation of functions
  %dx=2;
  conv1=1/(100);% to go from 1/m to 1/cm
  conv2=1/100^2;% to go from 1/m^2 to 1/cm^2
  conv3=1/100^3;
  Nmax=225*2; % run simulation over a year period
  Xmax=500; % horizontal domain size
  Ymax=500; % vertical domain size
  %I0=zeros(Xmax,Ymax); %Daily irradiance levels
  %T=zeros(Xmax,Ymax);  %Daily temperature
  Wold=zeros(Xmax,Ymax);  %Daily milfoil growth
  Wnew=zeros(Xmax,Ymax);
  S=zeros(Xmax,Ymax);  %Daily stolon growth
  lambda=zeros(Xmax,Ymax);%Daily loss due to respiration
  Id=zeros(Xmax,Ymax); %Irradiance at lake bottom
  Im=zeros(Xmax,Ymax); %Irradiance at depth of plant top
  h=zeros(Xmax,Ymax);  %height of plant
  W0=5*conv2; % initial condition for milfoil
  
  % % Model parameters taken from Herb 2006 and 
  U0=0.1;
  %U0=0.5;
  thetag=1.06;
 
  Tb=10;
  Kwt=0.5;
  km=0.006;
  Ps=50;
  d=2; % water depth in meters
  k1=90;

  lambda0=0.01;
  thetar=1.072;
  delta=0.0;
  
  dt= 0.5; % time step is one day  (I made this 1/2 day, so I run the simulation twice as long)
  dx = 2.0;% Spatial step is 2 cm
  D = 0.1; % diffusion coefficient is in cm^2/day (max is 1.9 in terms of dt and dx)
  %c0= 0.4; % 1/2 max. Make this bigger to allow more time for milfoil to grow (it's growing too fast if not)
  cstar= 0.3; % invasion speed of the moving front (largest 1 cm/day) - assume this is the largest value (let's use a mean instead of the max)
  r= cstar.^2 /(2*D);
  %r=1; % intrinsic growth rate of the stolons. Here, cstar = sqrt(2*D*r), this is the asymptotic wave speed.
  k=1; %carrying capacity of the stolons is 1 (initial s can't be bigger than this)
 
  
  %pulling factors
  pulling_day=90;
  stolon_percent_left=1; %the percentage of stolons left after pulling based on what was already there
  
  DStabiliity = 0.25*(dx^2/dt); % D must be smaller than this
   if D > DStabiliity
     disp('D too large')
     return
 end
  
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

% % Model will  be more realistic if actual data used.
%% Initailization



% preallocation of W (milfoil) and S (stolons)

 Wold(:,:) = W0;  %must define this an non-zero everywhere!
 S(:,:) = 0;
%    
%initialize stolons 
    i=1;
    while(i~=10)
    S(215,201:230) = 0.5;
   
    S(215+i,201+i:230-i)=0.5;
   
    S(215-i,201+i:230-i)=0.5;
    i=i+1;
    end

% S(250:450,250:450)=0.5;  big box


%% Time and 2D space Integration
n=1;
% step through time loop
while(n < Nmax)
 
I0=50*sin(pi*(n*dt+15)/360);  % surface irradiance is constant along entire surface (1D for now)
T=25*sin(pi*(n*dt-20)/360)*sin(pi*(n*dt-20)/360); %water temperature is constant
lambda=lambda0*thetar.^(T-Tb);

%Height
h(:,:)=min(d,(Wold(:,:)/conv2)/Ps);
%Irradiance at each time
  Im(:,:)= I0*exp(+Kwt*d-Kwt*h(:,:));
  Id(:,:)= I0*exp(-Kwt*d-(km/conv2*Wold(:,:)));
 


% % S Numerical Intergration (diffusion eq order: verticle, horizontal,
% normal diagonal, backward diagonal)
S(2:499,2:499) = S(2:499,2:499) + dt*r*S(2:499,2:499).*(1-(S(2:499,2:499))) + dt/dx^2*D*(S(3:500,2:499)-2*S(2:499,2:499)+S(1:498,2:499))+dt/dx^2*D*(S(2:499,3:500)-2*S(2:499,2:499)+S(2:499,1:498));

% S No-flux boundary conditions
S(500,:) = S(499,:);
S(1,:) = S(2,:); 

S(2:499,500) = S(2:499,499);
S(2:499,1) = S(2:499,2);



% %Hill functional response
for i = 1:500
    for j=1:500
        
        if(S(i,j)<10^(-5))
            Wnew(i,j) = Wold(i,j); % no change if Stolon density very low
        else

% W Numerical Integration
% % change only if stolon density big enough
Wnew(i,j) = Wold(i,j) + S(i,j).*dt*(U0*thetag^(T-Tb)*Wold(i,j))./((Kwt*h(i,j)+(km/conv2*Wold(i,j)))).*log((k1+Im(i,j))./(k1+Id(i,j)))-dt*(lambda + delta).*Wold(i,j);
        end
    end
end


if(n==1)
 i_space_w_milfoil=0; 
 
   for i = 1:500
    for j=1:500
        
        if(Wnew(i,j)>W0)
           
            i_space_w_milfoil= i_space_w_milfoil+1;
        else
            
        end
    end
   
   end
end

if(n==pulling_day*2) % mid season distribution
    figure(5)
    mesh(Wnew)
    title(['Biomass Distribution before pulling (day ' num2str(pulling_day) ')']) 
    zlabel('Biomass (g/cm^2)')
    axis([0 500 0 500 0 0.035 0 0.035])
    colorbar
    set(gca,'FontSize',18)
  
    figure(6)
    mesh(S)
    title(['Stolon Distribution before pulling (day ' num2str(pulling_day) ')'])
    axis([0 500 0 500 0 1 0 1])
    colorbar
    set(gca,'FontSize',18)
    
 
%pull
     %down the middle of patch  
   invaded_before = numel(find(S>10^(-5)))
   inv_before_area=invaded_before*4
   total_mil_before =0;
    for i = 1:500
    for j=1:500
        
        if(S(i,j)>10^(-5))
            total_mil_before = total_mil_before +Wnew(i,j); % no change if Stolon density very low
            
        else
            
        end
    end
    end 
   before_grid_ave= total_mil_before/(500*500) 
   before_patch_ave= total_mil_before/invaded_before
   %down the middle of patch
     %make half milfoil grid for area to pull
    half_graph_W=Wnew(:,1:215);
     %find out how many spaces have milfoil and organize into cooridnates
     [row1,column1]=find(half_graph_W>W0);
    
     
     squares_denuted=numel(row1)
     
    
    
   %removal
   

    for i=1:squares_denuted
        
        S(row1(i),column1(i))= S(row1(i),column1(i))*stolon_percent_left;
    
        i=i+1;
    end

    Wnew(:,1:215)=W0;
       


invaded_after = numel(find(S>10^(-5)))
    
      
%     figure(7)
%     mesh(Wnew)
%     title(['Biomass Distribution after pulling (day ' num2str(pulling_day) ')'])
%     zlabel('Biomass (g/cm^2)')
%     axis([0 500 0 500 0 0.035 0 0.035])
%     set(gca,'FontSize',18)
  
    figure(8)
    mesh(S)
    title(['Stolon Distribution after pulling (day ' num2str(pulling_day) ')'])
    axis([0 500 0 500 0 1 0 1])
    set(gca,'FontSize',18)
 end

n=n+1;
Wold=Wnew;
end

% final distributions
    figure(9)
    mesh(Wnew)
    title(['End of Season Biomass (D=' num2str(D) ' cm^2/day)']) 
    zlabel('Biomass (g/cm^2)')
    axis([0 500 0 500 0 0.035 0 0.035])
    colorbar
    set(gca,'FontSize',18)
  
    figure(10)
    mesh(S)
    title(['End of Season Stolon Distribution (D=' num2str(D) ' cm^2/day)'])
    axis([0 500 0 500 0 1 0 1])
    colorbar
    set(gca,'FontSize',18)

% 
% 
 %% averaging over space (not including the initlizd milfoil
 W_total=0;
 space_w_milfoil=0;
 
 for i = 1:500
    for j=1:500
        
        if(S(i,j)>10^(-5))
            W_total = W_total +Wold(i,j); % no change if Stolon density very low
            space_w_milfoil= space_w_milfoil+1;
        else
            W_total= W_total;
        end
    end
end
W_ave_in_total_grid= W_total/(Xmax*Ymax) 
W_ave_within_patch= W_total/space_w_milfoil

% returns maximum value of milfoil growth in a single patch in the domain
% (should match the ODE model when D is very low, since r is very high and
% S will reach it's max of 1.

highest= max(max(max(Wnew,[],500)))
space_w_milfoil
per= (space_w_milfoil/(500*500))*100
area_cov= space_w_milfoil*4
