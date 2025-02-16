%MotorAxonModel_1D_05072023



%{
This program is focused on modeling the transport of motor proteins through
the neuronal axons. We will be utilizing experimental parameters to simulate our
traffic model. 

%}

%Author: Roy Cruz Jr
%Date: May 7 2023

%Version 8: Motor Axon Model


clear all;
close all;


%Parameter values that can be changed
NMT=300; % Number of Microtubule (MT) present in the axon
N=10; % Number of consecuive time steps
M=100; % Number of motors
dt=0.1; % Time step (sec)
vnot=1.0; % Motor speed (microns/sec)
probplus=.4; % Probability of MT plus end leading
Laxon=100; % Length of axon (microns)
koff=.2; % Rate of switching to another MT





R=zeros(N,7,M);

%{
The array R will have a total of five columns. Each column will either store
data with respect to time(1), (2)MT numbeer which a motor protein is attached to,
(3)orientation of the MT that has a motor protein attached to it,(4) position of
the MT based on the axon,(5) the position of the motor protein based on the
MT it is attached to, (6)and number of binding sites available, and (7)
represents the position of the motor protein based on axon length.
%}


ksizestep=.008; %the size of each step taken from kinesin motors (microns)

% orientation of all MT's and there positon based on axon length. The last
% column is just telling us what number MT it is.

P=zeros(NMT,1);
iv=zeros(NMT,1);
%iv(:,3)=transpose(1:NMT);
 
% randomize MT polarity and positon based on the length of the axon

for i=1:NMT
    if (rand<probplus)
        P(i,1)=1; %MT is plus leading
        Mant=find(P==1);
        iv(i,1)=randi(Laxon);
    else
        P(i,1)=2; % MT is plus leading
        Mret=find(P==2);
        iv(i,1)=randi(Laxon);
    end
end



% Create an array that will store the MT length and number of binding sites

Lave=10;
mt_length=zeros(NMT,1);
mtnum_bs=zeros(NMT,1);


for j=1:NMT
    %mt_length(j,1)=randi(Lave);
    mt_length(j,1)=Lave;
    mtnum_bs(j,1)=round((mt_length(j,1))/ksizestep);
end

% Assign which MT will have Motor Protein (MP) attach to it. Also, assign a
% random position on the length of the MT attached to it.


mtwmp=zeros(M,1);
mp_loc=zeros(M,1);


for i=1:M
    %mtwmp(i,1)=randi(NMT); R(1,2,i)=i;
    R(1,2,i)=i;
    R(1,3,i)=P(i,1);
    R(1,5,i)=randi(mt_length(i,1));
    R(1,4,i)=iv(R(1,2,i),1);
    R(1,7,i)=R(1,4,i)+R(1,5,i);
    if (P(i,1)==2)
        R(1,6,i)=round((R(1,5,i))/ksizestep);
    else
        R(1,6,i)=round((mt_length(i,1)-R(1,5,i))/ksizestep);
    end
end


% create middle region polarity flaw
% This can be adjusted to either proximal, distal or all regions instead.
% Tracking MT by assigning them a number


Pregion=round(Laxon/3); % proximal region distance
Mregion=round(Laxon/3)*2; % middle region distance
Dregion=round(Laxon/3)*3; % distal region distance
ProximalArray=zeros(NMT,2); % MTs in the proximal region
MiddleArray=zeros(NMT,2); % MTs in the middle region
DistalArray=zeros(NMT,2); % Mts in the distal region
ProximalPopulation = zeros(M,1); % MT population in the proximal region
MiddlePopulation = zeros(M,1); % MT population in the middle region
DistalPopulation = zeros(M,1); % MT population in the distal region


for B=1:NMT
    if(iv(B,1)<=Pregion)
        ProximalArray(B,1)=iv(B,1);
        ProximalArray(B,2)=B;
    

    elseif(Pregion<iv(B,1) && iv(B,1)<=Mregion)
        MiddleArray(B,1)=iv(B,1);
        MiddleArray(B,2)=B;
   

    elseif(Mregion<iv(B,1) && iv(B,1)<=Dregion+1)
        DistalArray(B,1)=iv(B,1);
        DistalArray(B,2)=B;
    end
end

% write a input command that allows user to chooose where the polarity flaw
% is located


disp('Choose a region: Proximal=1, Middle=2, Distal=3 or none=4')
Region=input('Choose A Value ')

for E=1:round(M)
    if(Region==1)
        if(ProximalArray(E,1)==iv(E,1))
            R(1,3,E)=2;
          
        end
    elseif(Region==2)
        if(MiddleArray(E,1)==iv(E,1))
            R(1,3,E)=2;
            
        end
    elseif(Region==3)
        if(DistalArray(E,1)==iv(E,1))
            R(1,3,E)=2;
           
        end
    end
end


% new array that saves the information about where the MT are located at

Population1=nonzeros(ProximalArray(:,2));
Population2=nonzeros(MiddleArray(:,2));
Population3=nonzeros(DistalArray(:,2));


disp('Choose a jamming condition: Scenario_A=1, Scenario_B=2, Scenario_C=3 or Scenario_D=4  ')
jamming=input('Choose A Value: ')
mp_switch=0;
countJams=0;
countingswitches=0;

% teehee motor protein go brrrrrrr (Calving wrote this - Roy)
for i=1:N-1
    for j=1:M
        if( R(i,7,j)<Laxon)
            if(R(i,6,j)>0)
                if(rand<koff*dt)
                    countingswitches=countingswitches+1;
                    if(R(i,7,j)<=Pregion)
                        for A=1:length(Population1)
                            if((R(i,7,j)-iv(Population1(A,1),1))<=10 || R(i,2,j)~=Population1(A,1))

                        mp_switch=1;
                        R(i+1,2,j)=Population1(A,1);
                        R(i+1,3,j)=P(R(i+1,2,j),1);
                        R(i+1,4,j)=iv(R(i+1,2,j),1);
                        %R(i+1,5,j)=rand*5;
                        R(i+1,5,j)=abs(randi(mt_length(R(i+1,2,j),1)));
                        R(i+1,7,j)=R(i+1,4,j)+R(i+1,5,j);

                        if(R(i+1,3,j)==1)
                            R(i+1,6,j)=round((mt_length(R(i+1,2,j),1)-R(i+1,5,j))/ksizestep);
                        else
                            R(i+1,6,j)=round(R(i+1,5,j)/ksizestep);
                        end
                        break;    
                            end
                        end
                        

                    elseif(Pregion<R(1,7,5) && R(1,7,5)<=Mregion)
                        for B=1:length(Population2)
                            if((R(i,7,j)-iv(Population2(B,1),1))<=10 && R(i,2,j)~=Population2(B,1))


                        mp_switch=1;
                        R(i+1,2,j)=Population2(B,1);
                        R(i+1,3,j)=P(R(i+1,2,j),1);
                        R(i+1,4,j)=iv(R(i+1,2,j),1);
                        %R(i+1,5,j)=rand*5;
                        R(i+1,5,j)=abs(randi(mt_length(R(i+1,2,j),1)));
                        R(i+1,7,j)=R(i+1,4,j)+R(i+1,5,j);

                        if(R(i+1,3,j)==1)
                            R(i+1,6,j)=round((mt_length(R(i+1,2,j),1)-R(i+1,5,j))/ksizestep);
                        else
                            R(i+1,6,j)=round(R(i+1,5,j)/ksizestep);
                        end
                        break;
                            end
                        end

                    elseif(Mregion<R(i,7,j) && R(i,7,j)<=Dregion+1)
                        for C=1:length(Population3)
                            if((R(i,7,j)-iv(Population3(C,1),1))<=10 && R(i,2,j)~=Population3(C,1))

                        mp_switch=1;
                        R(i+1,2,j)=Population3(C,1);
                        R(i+1,3,j)=P(R(i+1,2,j),1);
                        R(i+1,4,j)=iv(R(i+1,2,j),1);
                        %R(i+1,5,j)=rand*5;
                        R(i+1,5,j)=abs(randi(mt_length(R(i+1,2,j),1)));
                        R(i+1,7,j)=R(i+1,4,j)+R(i+1,5,j);

                        if(R(i+1,3,j)==1)
                            R(i+1,6,j)=round((mt_length(R(i+1,2,j),1)-R(i+1,5,j))/ksizestep);
                        else
                            R(i+1,6,j)=round(R(i+1,5,j)/ksizestep);
                        end
                        break;
                            end
                        end

                    end
                end
            end
        end

         if (R(i,3,j)==1 && mp_switch==0 &&  R(i,7,j)<Laxon && R(i,6,j)>0 && R(i,5,j)<mt_length(R(i,2,j),1) && mt_length(R(i,2,j),1)-R(i,5,j)>=ksizestep)
            R(i+1,2,j)=R(i,2,j);
            R(i+1,3,j)=R(i,3,j);
            R(i+1,4,j)=R(i,4,j);
            R(i+1,5,j)=R(i,5,j)+ksizestep;
            R(i+1,6,j)=R(i,6,j)-1;
            R(i+1,7,j)=R(i,7,j)+ksizestep;


        elseif (R(i,3,j)==2 && mp_switch==0 &&  R(i,7,j)<Laxon && R(i,6,j)>0 && R(i,5,j)<mt_length(R(i,2,j),1) && mt_length(R(i,2,j),1)-R(i,5,j)>=ksizestep)
            R(i+1,2,j)=R(i,2,j);
            R(i+1,3,j)=R(i,3,j);
            R(i+1,4,j)=R(i,4,j);
            R(i+1,5,j)=R(i,5,j)-ksizestep;
            R(i+1,6,j)=R(i,6,j)-1;
            R(i+1,7,j)=R(i,7,j)-ksizestep;

        end



        if (R(i,6,j)<=0  && R(i,7,j)<=Laxon || (mt_length(R(i,2,j),1)-R(i,5,j))>=ksizestep)
%             disp('HERE')
            if(R(i,7,j)<=Pregion)
                for A=1:length(Population1)
                    if((R(i,7,j)-iv(Population1(A,1),1))<=5 && R(i,2,j)~=Population1(A,1))

                R(i+1,2,j)=Population1(A,1);
                R(i+1,3,j)=P(R(i+1,2,j),1);
                R(i+1,4,j)=iv(R(i+1,2,j),1);
                %R(i+1,5,j)=rand*5;
                R(i+1,5,j)=abs(randi(mt_length(R(i+1,2,j),1)));
                R(i+1,7,j)=R(i+1,4,j)+R(i+1,5,j);

                if(R(i+1,3,j)==1)
                    R(i+1,6,j)=round((mt_length(R(i+1,2,j),1)-R(i+1,5,j))/ksizestep);
                else
                    R(i+1,6,j)=round(R(i+1,5,j)/ksizestep);
                end
                break
                    end
                end

            elseif(Pregion<R(i,7,j) && R(i,7,j)<=Mregion)
                for B=1:length(Population2)
                    if((R(i,7,j)-iv(Population2(B,1),1))<=5 && R(i,2,j)~=Population2(B,1))

                R(i+1,2,j)=Population2(B,1);
                R(i+1,3,j)=P(R(i+1,2,j),1);
                R(i+1,4,j)=iv(R(i+1,2,j),1);
                %R(i+1,5,j)=rand*5;
                R(i+1,5,j)=abs(randi(mt_length(R(i+1,2,j),1)));
                R(i+1,7,j)=R(i+1,4,j)+R(i+1,5,j);

                if(R(i+1,3,j)==1)
                    R(i+1,6,j)=round((mt_length(R(i+1,2,j),1)-R(i+1,5,j))/ksizestep);
                else
                    R(i+1,6,j)=round(R(i+1,5,j)/ksizestep);
                end
                break
                    end
                end

            elseif(Mregion<R(i,7,j) && R(i,7,j)<=Dregion+1)
                for C=1:length(Population3)
                    if((R(i,7,j)-iv(Population3(C,1),1))<=5 && R(i,2,j)~=Population3(C,1))

                R(i+1,2,j)=Population3(C,1);
                R(i+1,3,j)=P(R(i+1,2,j),1);
                R(i+1,4,j)=iv(R(i+1,2,j),1);
                %R(i+1,5,j)=rand*5;
                R(i+1,5,j)=abs(randi(mt_length(R(i+1,2,j),1)));
                R(i+1,7,j)=R(i+1,4,j)+R(i+1,5,j);

                if(R(i+1,3,j)==1)
                    R(i+1,6,j)=round((mt_length(R(i+1,2,j),1)-R(i+1,5,j))/ksizestep);
                else
                    R(i+1,6,j)=round(R(i+1,5,j)/ksizestep);
                end
                break
                    end
                end
            end


        end

         R(i+1,1,j)=R(i,1,j)+1; % time array
         %countingswitches=countingswitches+1;


         %replace all negative positions and number of biding sites
         %available with a zero instead.

         % An alternative if statement could be if(R(i+1,7,j)>=Laxon)

         if (R(i+1,7,j)>=Laxon)
           R(i+1,2,j)=R(i,2,j);
            R(i+1,3,j)=R(i,3,j);
            R(i+1,4,j)=R(i,4,j);
            R(i+1,5,j)=R(i,5,j);
            R(i+1,6,j)=0;
            R(i+1,7,j)=Laxon;

         end
         

        if (R(i+1,7,j)<=0)
            R(i+1,7,j)=0;
        end


    end



        % This part of the code is comparing two motor proteins location (on
        % the same MT) and position in the axon.  Then we will implement the
        % jamming condition


        for r=1:M-1
            for z=r+1:M
                if(R(i+1,2,r)==R(i,2,z) && R(i+1,5,r)-R(i,5,z)<=ksizestep)
                    if(jamming==1)
                        R(i,2,z)=0; % previous mp is floating
                        R(i,3,z)=0;
                        R(i,4,z)=0;
                        R(i,5,z)=0;
                        R(i,6,z)=0;
                        R(i,7,z)=R(i,7,z);

                        countJams=countJams+1;

                    elseif (jamming==2)
                        R(i+1,2,r)=0; % new mp floating
                        R(i+1,3,r)=0;
                        R(i+1,4,r)=0;
                        R(i+1,5,r)=0;
                        R(i+1,6,r)=0;
                        R(i+1,7,r)=R(i+1,7,r);

                        countJams=countJams+1;

                    elseif(jamming==3)
                        R(i,2,z)=0; % previous mp floating
                        R(i,3,z)=0;
                        R(i,4,z)=0;
                        R(i,5,z)=0;
                        R(i,6,z)=0;
                        R(i,7,z)=R(i,7,z);

                        R(i+1,2,r)=0; % new mp floating
                        R(i+1,3,r)=0;
                        R(i+1,4,r)=0;
                        R(i+1,5,r)=0;
                        R(i+1,6,r)=0;
                        R(i+1,7,r)=R(i+1,7,r);

                        countJams=countJams+1;

                    elseif(jamming==4)
                        R(i+1,2,r)=R(i,2,r); % prevent move from happening
                        R(i+1,3,r)=R(i,3,r);
                        R(i+1,4,r)=R(i,4,r);
                        
                        R(i+1,6,r)=R(i,6,r)-1;
                        
                        if (R(i+1,2,r)==1)
                            R(i+1,5,r)=R(i,5,r)+ksizestep;
                            R(i+1,7,r)=R(i,7,r)+ksizestep;
                        elseif(R(i+1,2,r)==2)
                            R(i+1,5,r)=R(i,5,r)-ksizestep;
                            R(i+1,7,r)=R(i,7,r)-ksizestep;
                        end

                        countJams=countJams+1;

                    end
                end
            end
        end
end

%show the number of jams that occur

countJams
            
% plot the results




% Calculate fraction of motors that have traveled # nm along the axon
count=0;
for C=1:M
    if(R(N,7,C)>=Laxon)
        count=count+1;
    end
end

fractionsuccess=count/length(R(N,7,:))*100






%Plot trajectories

for i=1:(M)-1
    for j=1:N-1
        CM=jet(M);

        plot(R(j,1,i),R(j,7,i),'.','MarkerSize',10,'Color',CM(i,:));
        plot(R(j,1,i),R(j,7,i),'.','MarkerSize',10,'Color','k');
        xlabel('time')
        ylabel('MP position (microns)')
        hold on;
        grid on;


        %Traching trajectory of each particle at each time step

        figure(i)

        dots=plot(R(j,1,i),R(j,7,i),'.','MarkerSize',12,'Color',CM(i,:));
        dots=plot(R(j,1,i),R(j,7,i),'.','MarkerSize',12,'Color','k');
        axis manual;

        p.Xdata=R(j,1,i);
        p.Ydata=R(j,7,1);
        drawnow;
        axis([0,N,0,Laxon]);

        line([R(j,1,i),R(j+1,1,i)],[R(j,7,i),R(j+1,7,i)],'color',CM(i,:));
        line([R(j,1,i),R(j+1,1,i)],[R(j,7,i),R(j+1,7,i)],'color','k');

    end
end
        

% % Plot frequency plot
% 
% figure
% hist(R(N,7,:))
% xlabel('time')
% ylabel('MP position (microns)')
% 







% need to check why sometimes it jumps so far away, start by checking the
% region distances first

% make similar figures as in the axonal code, in order to show the
% distribution of motor proteins in each region from the initial and final
% datapoints.


% make a figure Fraction success Vs. length of axon



%need to fix the jamming condition by adding fixing the 3,4,5,6, and 7th
%column i think. I am not sure.



% need to fix the if(rand<koff* dt) because an alternative could be better.


%check out the processivity of the kinesin motors when they are detachment


% add a new if condition 


% debug the code a bit more, and work on the poster


