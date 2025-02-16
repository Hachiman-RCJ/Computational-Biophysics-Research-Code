% Traffic jam code draft

%Parameter values that can be changed
NMT=500; % Number of MT present in the axon
N=6000; % Number of consecuive time steps
M=250; % Number of motors
P=zeros(NMT, 1); % Polarity of MTs
dt=0.1; % Time step (sec)
vnot=1.0; % Motor speed (microns/sec)
probplus=0.0; % Probability of MT plus end leading
Laxon=(1000); % Length of axon (microns)
koff=0.1; % Rate of switching to another MT
R=zeros(N,4,M); % column 1: time, column 2: the MT # which the motor is attached to
% column 3: orientation of the MT that the motor is attached to
%column 4: the position of the motor in the axon
ksizestep=8; %the size of each step taken from kinesin motors (microns)

% Randomize MT polarity
% Plus leading=2 and minus leading =3 and still floating around = 1
% the relation operator will decide if it is plus or minus leading

for k=1:NMT
    if(rand<probplus)
        P(k,1)=2; % MT is plus leading
        Mant=find(P==2);
        SizeMant=numel(Mant);
        randi(SizeMant);
    else
        P(k,1)=3; % Mt is minus leading
        Mret=find(P==3);
        SizeMret=numel(Mret);
        randi(SizeMret);
    end


end 


% give motor proteins an initial attachment
% Also give each motor a random position in the axon
for A=1:1
    for B=1:M
        R(A,2,B)=randi(numel(P)); %choosing random MT location
        R(A,3,B)=P(R(A,2,B),1); % the orientation of the random choosen MT
        R(A,4,B)=randi(Laxon); % randomly assign an initial position to motor



    end 
end


% The chance of an axon region to have a polarity flaw
% such that the region is has mostly minus end out MTs

Pregion=round(Laxon/3); % proximal region distance
Mregion=round(Laxon/3)*2; % middle region distance
Dregion=round(Laxon/3)*3; % distal region distance
ProximalArray=zeros(M,1); % MTs in the proximal region
MiddleArray=zeros(M,1); % MTs in the middle region
DistalArray=zeros(M,1); % Mts in the distal region
ProximalPopulation = zeros(M,1); % MT population in the proximal region
MiddlePopulation = zeros(M,1); % MT population in the middle region
DistalPopulation = zeros(M,1); % MT population in the distal region



% which MTs are located in each of the three regions
for A=1:1
    for B =1:M
 if (R(1,4,B)<=Pregion) % checking which MTs are in the proximal region
     ProximalArray(B,1)=R(A,2,B); % storing which MTs are in the proximal region
    
 end
    end
end


for A=1:1
    for B =1:M
 if (Pregion<=R(1,4,B)) && (R(1,4,B)<= Mregion) % checking which MTs are in the middle region
MiddleArray(B,1)=R(A,2,B); % storing which MTs are in the middle region
 end
    end
end


for A=1:1
    for B =1:M
 if (Mregion<=R(1,4,B)) && (R(1,4,B)<= Dregion) % checking which MTs are in the distal region
DistalArray(B,1)=R(A,2,B); % storing which Mts are in the distal region
 end
    end
end


% Polarity flaw located in the middle region

for E=1:numel(MiddleArray)
    if MiddleArray(E,1)>0
        R(1,3,E)=3;
    end
end

% population from all three regions and the average overall of

% Middle region
MiddlePopulation = MiddleArray;
for H=1:numel(MiddleArray)
    for I=1:numel(MiddlePopulation)
        if I <= numel(MiddlePopulation)
             I;
        numel(MiddlePopulation);
            if MiddlePopulation(I,1)==0
            MiddlePopulation(I,:)=[];
            
            end
        end
     
    end
end

Mpop=numel(MiddlePopulation);


% Proximal region
ProximalPopulation = ProximalArray;
for H=1:numel(ProximalPopulation)
    for I=1:numel(ProximalPopulation)
        if I <= numel(ProximalPopulation)
            if ProximalPopulation(I,1)==0
            ProximalPopulation(I,:)=[];
            
        end
    end
end
end

Ppop=numel(ProximalPopulation);

% Distal region
DistalPopulation = DistalArray;
for H=1:numel(DistalPopulation)
    for I=1:numel(DistalPopulation)
        if I <= numel(DistalPopulation)
            if DistalPopulation(I,1)==0
            DistalPopulation(I,:)=[];
            
        end
    end
end
end

Dpop=numel(DistalPopulation);

%Average of all three regions

Averagepop=(Mpop+Ppop+Dpop)/3;


%The following for loop will help calculate the location of each motor
%protein at every time step when they are attaching to different MT


for i=1:N % Consecutive steps
    for j=1:M % cycling all motors
            if rand<koff*dt
                R(i+1,2,j)= randi(numel(P)); % motor j attaches to a random MT
                R(i+1,3,j) = P(R(i+1,2,j),1); % the orientation of MT which is attached to motor j
            else
                R(i+1,2,j)=R(i,2,j); % motor protein stays where it was
            R(i+1,3,j) = R(i,3,j); % motor protein remains floating in the abyss...
                
            end

        if R(i+1,3,j)== 2 % Motor on plus end out MT
           
            R(i+1,4,j) = R(i,4,j)+vnot*dt;
           
        else
             R(i+1,4,j) = R(i,4,j)-vnot*dt;
            
        end
        
        if  R(i+1,4,j)>Laxon % motor reaches axonal terminal
             R(i+1,4,j)=Laxon; % stays at this location
        end
        
        
        if  R(i+1,4,j)<0
             R(i+1,4,j)=0;
        end

        w=1;
    
        % This keeps track of how many iterations were completed

         R(i+1,1,j)=R(i,1,j)+1; % time array
    end

% this part of the code is comparing two motor proteins location (on the same MT)
% and position in the axon. Then we will implement the jamming condition

    for r=1:M-1
            for z=r+1:M
                if(R(i+1,2,r)==R(i,2,z))
                    if (R(i+1,4,r)-R(i,4,z)) <= ksizestep
                        R(i+1,2,j)=R(i,2,j);
                        R(i+1,4,r)=R(i,4,r);

                    end
                end
            end
        end
    end


%find the population in each region at the final timestep for each motor
%recording the population of each motor at the last timestep

% The chance of an axon region to have a polarity flaw
% such that the region is has mostly minus end out MTs

FProximalArray=zeros(M,1); % MTs in the proximal region
FMiddleArray=zeros(M,1); % MTs in the middle region
FDistalArray=zeros(M,1); % Mts in the distal region
FProximalPopulation = zeros(M,1); % MT population in the proximal region
FMiddlePopulation = zeros(M,1); % MT population in the middle region
FDistalPopulation = zeros(M,1); % MT population in the distal region


% which MTs are located in each of the three regions
for A=N:N
    for B =1:M
 if (R(N,4,B)<=Pregion) % checking which MTs are in the proximal region
     FProximalArray(B,1)=R(A,2,B); % storing which MTs are in the proximal region
    
 end
    end
end


for A=N:N
    for B =1:M
 if (Pregion<=R(N,4,B)) && (R(N,4,B)<= Mregion) % checking which MTs are in the middle region
FMiddleArray(B,1)=R(A,2,B); % storing which MTs are in the middle region
 end
    end
end


for A=N:N
    for B =1:M
 if (Mregion<=R(N,4,B)) && (R(N,4,B)<= Dregion+1) % checking which MTs are in the distal region
FDistalArray(B,1)=R(N,2,B); % storing which Mts are in the distal region
 end
    end
end



% population from all three regions and the average overall of

% Middle region
FMiddlePopulation = FMiddleArray;
for H=1:numel(FMiddleArray)
    for I=1:numel(FMiddlePopulation)
        if I <= numel(FMiddlePopulation)
             I;
        numel(FMiddlePopulation);
            if FMiddlePopulation(I,1)==0
            FMiddlePopulation(I,:)=[];
            
            end
        end
     
    end
end

FMpop=numel(FMiddlePopulation);


% Proximal region
FProximalPopulation = FProximalArray;
for H=1:numel(FProximalPopulation)
    for I=1:numel(FProximalPopulation)
        if I <= numel(FProximalPopulation)
            if FProximalPopulation(I,1)==0
            FProximalPopulation(I,:)=[];
            
        end
    end
end
end

FPpop=numel(FProximalPopulation);

% Distal region
FDistalPopulation = FDistalArray;
for H=1:numel(FDistalPopulation)
    for I=1:numel(FDistalPopulation)
        if I <= numel(FDistalPopulation)
            if FDistalPopulation(I,1)==0
            FDistalPopulation(I,:)=[];
            
        end
    end
end
end

FDpop=numel(FDistalPopulation);

%Average of all three regions

FAveragepop=(FMpop+FPpop+FDpop)/3;



% Calculate fraction of motors that have traveled # nm along the axon
count=0;
for C=1:M
    if(R(N,4,C)>=Laxon)
        count=count+1;
    end
end

fractionsuccess=count/numel(R(N,4,:))


% initial and final population timestep

figure
xbar=categorical({'Proximal','Middle','Distal','Overall'});
xbar=reordercats(xbar,{'Proximal','Middle','Distal','Overall'});
ybar=[Ppop,Mpop,Dpop,Averagepop];
bar(xbar,ybar)
ylabel('Motor Population in each region')
title('Initial Motor Distribution')

figure
xbar=categorical({'Proximal','Middle','Distal','Overall'});
xbar=reordercats(xbar,{'Proximal','Middle','Distal','Overall'});
ybar=[FPpop,FMpop,FDpop,FAveragepop];
bar(xbar,ybar)
ylabel('Motor Population in each region')
title('Final Motor Distribution')




% Plot trajectories
figure 
for D=1:M
    plot(R(:,1,D),R(:,4,D))
    hold on
end


% Plot frequency plot

figure
hist(R(N,4,:))
