%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
% 
% Notes:
% - This function estimates the PDF of sampled data via an AKDE whose
%   bandwidth is varied depending on the value of the sampled data point 
%   (posAKDE). Smaller values result in smaller bandwidths. 
% - Hence, this function is only suitable for strictly decaying PDFs!
% - This function takes reasonably long for length(X)=10^5 (few minutes) 
%   and quite long for 10^6 (~10min).
% - The default value of k (coefficient for bandwidth) needs to be adapted
%   to the unit of X. E.g. for days it should be 0.5, for hours 12.
% 
% Input:
% - X: Sampled data of which we want to find the underlying pdf.
% - xi: Points at which we want to calculate the pdf. Note: They do not
%   have to be equispaced!
% - k: Factor for bandwidth that is per default set to 0.5*X(i). That is
%   depending on the position of the sampled data point its bandwidth will be
%   chosen.
% - BW: Adaptive bandwidth calculated in advance as
%     inds = Rhisto ~= 0; 
%     K = Thisto(inds); 
%     BW = ones(lentgh(K));
%     BW(2:end) = K(2:end) - K(1:end-1); 
%     BW(1) = BW(2);
%   Or: BW(i) = BW(find(abs(rwt(i)-T)<binwidth)) (tbc from notes)
% Output:
% - fi(xi): Point pairs for the estimated PDF of the sampled data.
% 
function [fi,xi] = posAKDE(X,xi,k,BW)
    
    xi = [-fliplr(xi) xi];
    [X,inds] = sort(X); %Sort sampled data for BW-calculation
    XX = length(X);
    fi = zeros(1,length(xi));
    
    if nargin <4
        BW = X;
    end
    
    if nargin >3
        BW = BW(inds); %Sort BW in the same way X was sorted
    end
    
    if nargin <3
        k = 0.5; %k should be 0.5 if X is in days and 12 if it is in hours!
    end
    
    if length(BW)~=length(X)
        error('Bandwidths and random variable must be of same lenght!')
    end
    
    for i = 1:XX
%         BW(i) = k*50*(X(i)-X(i-1))+k; %This line is quite arbitrary!
        a = normpdf(xi,X(i),k*BW(i));
%         if isnan(sum(a))
%             i
%             x = X(i)
%             norm = sum(a)
%         end
        fi = fi + 1/XX*a; %y = normpdf(x,mu,sigma)
        if mod(i/XX*100,1) 
            disp(['AKDE: ',num2str(i/XX*100),' %']) %Print progress
        end
    end
    %Post-reflection:
    r = fliplr(fi(xi < 0)); %ri-values to be reflected and added to the respective positive ri.
    pos = find(xi>0); pos = pos(1); %Position of the first positive ri-value.
    fi(pos:pos+length(r)-1)= fi(pos:pos+length(r)-1) + r; %Only update the first ri3 entries.
    fi(1:pos-1) = []; 
    xi(1:pos-1) = []; 
    
    beep on; beep; beep off; %Beep at the end of this lengthy computation.
end
   
        
