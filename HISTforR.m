%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
% 
% Notes:
% - This function obtains the correctly scaled resolvent R(t)
% - The exponential fit is performed as described in
%   http://www.matrixlab-examples.com/exponential-regression.html
%
% Input:
% - rwt: Renormalized response times
% - binwidth: /equiv Delta_t ~ error in ~R(s)!
% - n: branching ratio for correct scaling of R(t)
%
% Output:
% - [R,t] = Correctly scaled R(t) in t-resolution binwidth, where t is a row vector of
%   equidistant bar-centers of the histogram
%
function [R,t] = HISTforR(rwt,binwidth,n)
    bins = max(rwt)/binwidth; %#bins should be at least 10*max(rwt) for binwidth to be small enough as it is the Delta_t for the num_laplace(). Delta_t must be very small (<~-0.1) for ~R(s) to good enough.
    rwt = rwt(rwt~=0); %Clean rwt from its zero entries. They correspond to exo-events waiting for themselves which does not make sense and distorts the following histogram.
    [R,t] = hist(rwt,bins); 
    K = sum(R)*binwidth*(1-n); %Constant of rescaling to 1/(1-n) because Int[R_tru(t),{t,0,inf}]=1/(1-n) s.t. 1/K Int[R_emp(t),{t,0,inf}] != 1/(1-n)
    R = R/K; %Rescaling
end
