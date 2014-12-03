%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
% 
% This function just plots the cluster statistics and returns the handle of
% the produced figure.
%
function [handle] = plot_ClustStats(tevnts,Tmax,pars_tru,Nexo,gens,NumGens,h,liftim)
    handle = figure;

        subplot(2,2,1)

            plot(tevnts,1:length(tevnts)); hold on;
            plot(0:Tmax,pars_tru(1)*(0:Tmax),'k--'); grid on; hold off;

            xlabel('t - Time'); 
            ylabel('N_t(t) - Counting process');
            legend(['Complete process with \Lambda_{emp}/\Lambda_{tru} = '...
                ,num2str(length(tevnts)/Tmax/pars_tru(1)*(1-pars_tru(2)))],...
                ['Exogeneous process with \mu_{emp}/\mu_{tru} = '...
                ,num2str(Nexo/Tmax/pars_tru(1))],...
                'Location','NorthWest');

        subplot(2,2,3)

            plot(tevnts,gens,'co'); 

            grid on;        
            xlabel('t - Time'); 
            ylabel('k(t) - Generation of event as proxy for current clusters size');
            legend(['\mu = ',num2str(pars_tru(1)),'. n = ',num2str(pars_tru(2)),'.']) %Just wanted to have that printed somewhere in the plots.

        subplot(2,2,2)

            semilogy(0:length(NumGens)-2,NumGens(2,1:end-1),'r.-'); %Exclude last entry of NumGens because it is always zero (k_{max}+1 = 0).

            grid on; 
            xlabel('k - Generation'); 
            ylabel('N_k(k) - Number of events per generation'); %Equivalent: h=hist(gens,max(gens)+1); %Count how many events occured per k-Generation

        subplot(2,2,4)
            h = hist(h,max(h));

            loglog(h,'mo'); grid on; 

            xlabel('s - Size of cluster'); 
            ylabel('N_s(s) - Number of occurrence of size'); 
            legend(['Lifetime of largest cluster: ',num2str(liftim)],'Location','NorthEast')
end
