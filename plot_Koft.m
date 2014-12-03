%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
% 
% This function takes the K-matrix produced by get_rwt() to produce an
% average K(t),i.e. K of t, that is compared to the theoretically predicted
% K(t)-curve. 
% 
% Input:
% - K: Matrix with each row corresponding to the independent clusters and
%   the n-th column corresponding to the longest time it took for the k-th
%   generation to occur since immigration.
% - pars: Parameters for plotting the theoretical K(t)-curve.
% 
% Output:
% - Kavg: K-matrix that is averaged over all rows.
% 
% Notes:
% - K(t)_emp and K(t)_theo won't perfectly agree on their scale probably
% because in the derivation many approximations have been made. What counts
% is that the slopes, i.e. the power law exponent theta, are the same.
% 
function [Kavg,Kerr] = plot_Koft(K,pars)
    KMAX = size(K,2);
    Kavg = zeros(1,KMAX);
    Kerr = zeros(1,KMAX);
    for j = 1:KMAX %for each column calculate average time
        Kavg(j) = mean(nonzeros(K(:,j)));
        Kerr(j) = std(nonzeros(K(:,j)))/(sqrt(numel(nonzeros(K(:,j)))));
    end
    figure; grid on; hold on;
        h1 = loglog(Kavg+1.96*Kerr,1:length(Kavg),'c-');
        h2 = loglog(Kavg-1.96*Kerr,1:length(Kavg),'c-');
        h3 = loglog(Kavg,1:length(Kavg),'b-');
        h4 = loglog((pars(3)*pars(4))/(1-pars(4))*(1:KMAX).^(1/pars(4)),1:KMAX,'g-'); 
%                     y1=Kavg+1.96*Kvar;                 %Create first curve
%                     y2=Kavg-1.96*Kvar;                 %Create second curve
%                     xx=[1:KMAX,fliplr(1:KMAX)];        %Create continuous xx value array for plotting
%                     yy=[y1,fliplr(y2)];                %Create yy-values for out and then back
%                     fill(xx,yy,[0.5 0.5 0.5]);         %Fill with grey
        ylabel('log[K(t)] - Current highest generation');
        xlabel('log(t) - Time since immigration');
        legend([h3,h4,h1],{'Empirical results','Theoretical prediction','0.95-confidence bounds'});
        title(['pars = [',num2str(pars(1)),', ',num2str(pars(2)),', ', num2str(pars(3)),', ', num2str(pars(4)),']']);
        set(gca,'XScale','log','YScale','log');
end
