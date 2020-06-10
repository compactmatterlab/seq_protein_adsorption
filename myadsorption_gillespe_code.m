tic
clear all
timepoints=(2:1:300);
p=0.5;
nSites=10;
timescale=1;

nStart = 1; %defines number of starting bonds
pMinus = p; % defines probability of breaking a bond
% pMinus = repmat(p,nSites,1)+0.04*randn(nSites,1);
pAdd = 1-p; %defines probability of adding a bond

samplesize = 1000000; %initialize how many proteins will be simulated
tTime = zeros(1,samplesize); %initialize total time matrix
%n = [];

for runs = 1:samplesize %runs through each protein that is to be simulated
    nBinding = nStart; %(re)sets # of bound sites to # of starting bound sites
    Time = 0; %initializing step counter
    %N = 0;
    while nBinding > 0 && Time < 100000000 %checks to see if protein still has at least 1 bond 
        X = rand; %creates random number to determine the action of the protein
        if nBinding < nSites %check to see if protein is fully bound to surface
            if X <= pMinus %check to see if random # within boundaries to break a bond
                nBinding = nBinding - 1; %breaks a bond
            elseif X > pMinus && X <= (pMinus+pAdd) %check to see if random # falls within boundary to create a bond
                nBinding = nBinding + 1; %creates a bond
            
            end
        else %if protein is in fully bound state 
            if X <= pMinus %check to see if random # falls within boundary to break bond
                nBinding = nBinding - 1; %breaks a bond
            elseif X > pMinus && X <= (pMinus+pAdd) %check to see if random # doesnt fall in "do nothing" boundary 
                nBinding = nBinding-1; %breaks a bond
            end
        end
        Time = Time + 1; %Protein step counter
        %N = N+1;
    end
    tTime(runs) = sum(exprnd(timescale,Time,1)); %multiplies # of steps taken by the time each step takes creating matrix of total residence time of proteins
end


A = timepoints; %sets A equal to same timepoints as the experimental data
n = zeros(1,length(A)); %intialize "# of proteins left at each time point" matrix
for i = 1:length(A) %runs through each timepoint
    L = tTime > A(i); %returns logical (1 or 0) depending on if value of tTime is greater than current timepoint
    M = find(L == 1); %finds where value of tTime was greater than current timepoint
    n(i) = length(M); %returns number of proteins that stuck for a time (tTime) longer than current timepoint (A(i))
end
A1 = A;
n=n'; %/n(1);
CRT = n/samplesize;


% Line = fit(A1',CRT,'power1')

% Line2 = fit(A1',CRT,'exp1')
% 
% Line3 = fit(A1',CRT,'exp2')


% hold on
% plot(Line(A1'),'b')
figure
plot(A1,CRT,'k-')

% plot(Line2(A),'r')
toc

for i = 1:100
    sampleno = randi(samplesize,1,20000);
    xtTime = tTime(sampleno(:));
    sampleavg(i)=mean(xtTime);
    samplestc(i)=std(xtTime);
end

figure
hist(sampleavg)