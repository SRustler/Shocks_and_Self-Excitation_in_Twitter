%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
% 
% Notes: 
% - This function only works for EXPONENTIAL kernels since only for
%   this one the analytical expression tau is known: R_{exp}(t) ~ exp[t *
%   (n-1)/tau] such that tau [in Phi(t) ~ exp(t/tau)] can be determined.
%
% Input:
% - N: number of events
% - Nexo: number of exogeneous events
% - Tmax: maximal time of simulation
% - b: characteristic time from the fit ~exp(b*x)
%
% Output:
% - pars_clu = [mu n tau] analogous to the pars needed for simulating the
%   data in the first place.
%
function [pars_clu] = clu_estim(N,Nexo,Tmax,b)
    pars_clu(1) = Nexo/Tmax; %mu
    pars_clu(2) = (N-Nexo)/N; %Calculate via n = Nendo/N
    pars_clu(3) = (pars_clu(2)-1)/b; %b is from the fit of R(t) = a*exp(b*x) = 1/tau * exp(-(1-n)/tau * t). 
%     pars_clu(3) = (1-pars_clu(2))*tren; %Recall: tren := tau_{ren} = tau/(1-n)
end
