function [CRT] = myadsorption_gillespe_opt(timepoints,experdata,timescale,p,nSites)

%nSites = 800;
nStart = 1; %defines number of starting bonds
pMinus = p; % defines probability of breaking a bond
pAdd = 1-p; %defines probability of adding a bond

samplesize = 1000000; %initialize how many proteins will be simulated
tTime = zeros(1,samplesize); %initialize total time matrix
%n = [];

for runs = 1:samplesize; %runs through each protein that is to be simulated
    nBinding = nStart; %(re)sets # of bound sites to # of starting bound sites
    Time = 0; %initializing step counter
    %N = 0;
    while nBinding > 0 && Time < 50000 %checks to see if protein still has at least 1 bond 
        X = rand; %creates random number to determine the action of the protein
        if nBinding < nSites; %check to see if protein is fully bound to surface
            if X <= pMinus; %check to see if random # within boundaries to break a bond
                nBinding = nBinding - 1; %breaks a bond
            elseif X > pMinus && X <= (pMinus+pAdd); %check to see if random # falls within boundary to create a bond
                nBinding = nBinding + 1; %creates a bond
            
            end
        else %if protein is in fully bound state 
            if X <= pMinus; %check to see if random # falls within boundary to break bond
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

f = experdata(1:end);
A = timepoints; %sets A equal to same timepoints as the experimental data
n = zeros(1,length(A)); %intialize "# of proteins left at each time point" matrix
for i = 1:length(A); %runs through each timepoint
    L = tTime > A(i); %returns logical (1 or 0) depending on if value of tTime is greater than current timepoint
    M = find(L == 1); %finds where value of tTime was greater than current timepoint
    n(i) = length(M); %returns number of proteins that stuck for a time (tTime) longer than current timepoint (A(i))
end
A1 = A;
n=n';

%section to find the number of events lasting longer than 1.25 to use to
%normalize the simulated data to match the normalization of experimental
%data
normalize = n(1)/f(1);
%end normalization section
CRT = n/normalize;


%Line = fit(A1,n1,'power1')

%Line2 = fit(A1,n,'exp1')

%Line3 = fit(A1,n1,'exp2')

%plot(Line(A),'b')
%plot(Line2(A),'y')
%plot(Line3(A),'g')