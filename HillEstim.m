%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
%
% Notes:
% - This function calculates the PL exponent theta via the Hill estimator
%   in dependence on the short time cutoff c.
% - The output is one element shorter than the input because the function
%   starts with a short time cutoff smaller than the largest time as
%   opposed to the largest time itself.
% - In order to determine the right theta and c, they have to be plotted
%   against each other and a region of convergence must be found.
% 
% Input:
% - X: Random variable that is PL distributed. Note it does not necessarily
%   have to follow a PDF. Further exponents are best recovered if the
%   distribution in the loglog-plot is a line on the entire range. E.g. a
%   c-shift would introduce a bending around t=c.
% 
% Output:
% - theta_hat: Estimated PL exponent versus short time cutoff
% - c_hat: Tested short time cutoff
% 
function [theta_hat,c_N] = HillEstim(X)
    X = sort(X,'descend'); 	
    num = length(X);
    theta_hat = zeros(num,1);
    c_N = zeros(num,1);
    for ic = 2:num
        c_N(ic) = X(ic);
        theta_hat(ic) = mean(log(X(1:ic)/c_N(ic)));
        theta_hat(ic) = 1/theta_hat(ic);
    end
    theta_hat=theta_hat(2:end);
    c_N=c_N(2:end);
end
