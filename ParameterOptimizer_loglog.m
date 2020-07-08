clear all;
clc;

load 0415_45pMFb_ManualCountsSurvivalFunctions.mat; %load experimental data
%load BSA647.mat; %load experimental data

exper_data1 = total_survFunc(:,2); %input experimental survival data values
timepoints1 = total_survFunc(:,1); %input experimental data timepoints

exper_data = exper_data1(2:length(exper_data1)-1); % If needed
timepoints = timepoints1(2:length(timepoints1)-1); % If needed
clearvars -except timepoints exper_data;

p = .4; %initialize probabilities (pMinus & pAdd)
timescale = 4; %initialize time each step takes
nSites = 5; %initialize number of possible binding sites
R2 = 1000000; %initialize coefficient of determination (used to determine fit quality)
[CRT] = myadsorption_gillespe(timepoints,exper_data, timescale, p, nSites);%runs the simulation
% R2 = (sum((log(CRT) - log(exper_data)).^2)) %calculates value of R2 for current fit no weight
R2 = (sum((abs(log(exper_data).^-1)).*(log(CRT) - log(exper_data)).^2)) %calculates value of R2 for current fit with weights
success=0;
fails = 0;
R2new = 10;
% tStep = .01;
% SitesStep = 2;
% pStep = .01;
avgjump=15; % average percentage change
jump = exprnd(avgjump); % actual percentage change

while fails < 100 %sets threshold to determine if data is a good enough fit
    
    timescalenew=timescale;
    pnew=p;
    nSitesnew=nSites;
    
    A = [timescale;p;nSites;fails;R2]
    
    RANDOM = rand;
    
    if RANDOM < .33
        RAND2 = rand;
        if RAND2 < .5
            timescalenew = (1+jump/100)*timescale;
        else
            timescalenew = (1-jump/100)*timescale;
        end
        if timescalenew<0.01
            timescalenew=0.01;
        end
        if timescalenew>2
            timescalenew=2;
        end
        [CRT] = myadsorption_gillespe_opt(timepoints,exper_data, timescale, p, nSites);%runs the simulation
        
    elseif RANDOM > .33 &&  RANDOM < .66
        RAND2 = rand;
        if RAND2 < .5
            nSitesnew = round((1+jump/100)*nSites);
        else
            nSitesnew = round((1-jump/100)*nSites);
        end
        if nSitesnew<5
            nSitesnew=5;
        end
        if nSitesnew>45
            nSitesnew=45;
        end
        
        [CRT] = myadsorption_gillespe_opt(timepoints,exper_data, timescale, p, nSites);%runs the simulation
        
    elseif RANDOM > .66
        RAND2 = rand;
        if RAND2 < .5
            pnew = (1+jump/100)*p;
        else
            pnew = (1-jump/100)*p;
        end
        if pnew<0.46
            pnew=0.46;
        end
        if pnew>0.54
            pnew=0.54;
        end
        
        [CRT] = myadsorption_gillespe_opt(timepoints,exper_data, timescale, p, nSites);%runs the simulation
    end
    %     R2new = (sum((log(CRT) - log(exper_data)).^2));
    R2new = (sum((abs(log(exper_data).^-1)).*(log(CRT) - log(exper_data)).^2));
    R2new
    %     loglog(timepoints, CRT, timepoints, exper_data)
    %     legend('simulated', 'experimental')
    if R2new<R2 || isinf(R2new)
        timescale=timescalenew;
        p=pnew;
        nSites=nSitesnew;
        R2=R2new;
        success=success+1;
        fails = 0;
    elseif isinf(R2new)
        fails = 0;
    else
        success=0;
        fails = fails + 1;
    end
    
    if success>10
        avgjump = avgjump*2;
    end
    if fails>40 && avgjump > 1
        fails=0;
        avgjump=avgjump/2;
    end
    jump = exprnd(avgjump); % actual percentage change
    
end
    R2
    loglog(timepoints, CRT, timepoints, exper_data)