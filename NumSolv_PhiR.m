%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
%
% Notes:
% - The quality of the numerical transformation will heavily depend on the
%   integration method (trapezoidal is better than rectangular), the
%   branching ratio and the binwidth (a higher n or bw amplifies the error).
% - Multimodality will be carried over to Phi such that this waviness will touch
%   negative Phi too at the long tail. Reducing the error means reducing
%   this waviness.
% 
% Input:
% - Ri: Correctly scaled resolvent!!!
% - t: Positive time stamps in R(t) that need to be equispaced!
% - n: Branching ratio used for the transformation to Phi.
% - method: Integration method. Either 'rect' for rectangular or 'trap' for
%   trapezoidal or 'simp' for the composite Simpson's rule (quadratic).
% Output:
% - time: Computation time for the entire routine
% - Phi: Numerically calculated bare kernel that should be normalized to 1.
%   Test via sum(Phi)/deltat=1?
% 
function [time, Phi] = NumSolv_PhiR(R,t,n,method)
    StartNum = tic;

    deltat = t(2)-t(1);
    Phi = zeros(size(t));
    
    Phi(1) = R(1); %Zeroth step: Phi(t0) = R(t0). Note that Ri itself is correctly scaled. 
    %Note that this equality only holds at t=0! If R(1) is actually at t>0 
    %this error will be kept for the entire calculation!

    if strcmp(method,'rect')
        for i = 2:length(t)
            sum = 0;
            for j = 1:i-1
                sum = sum + Phi(j) * R(i-j) * deltat;
            end
            Phi(i) = R(i) - n * sum;
            if mod(i,100) == 0
                disp(['NumSolv: ',num2str(i/length(t)*100,'%.2f'),' %']);
            end
        end
    end  
    
    if strcmp(method,'trap')
        for i = 2:length(t)
            sum = 0;
            for j = 1:i-2
                sum = sum + Phi(j) * R(i-j) * deltat;
            end
            Phi(i) = R(i) - n*( 0.5*Phi(1)*R(i-1)*deltat + sum + 0.5*Phi(i-1)*R(1)*deltat ); %The last term does not make sense since it is always R(0)!
            if mod(i,100) == 0
                disp(['NumSolv: ',num2str(i/length(t)*100,'%.2f'),' %']);
            end
        end
    end
   
    if strcmp(method,'simp')
        for i = 2:length(t)
            sum2 = 0;
            sum4 = 0;
            for j = 1:(i/2-1)
                sum2 = sum2 + 2 * Phi(2*j) * R(i-2*j) * deltat;
            end
            for j = 1:(i/2) 
                sum4 = sum4 + 4 * Phi(2*j-1) * R(i-2*j+1) * deltat;
            end
            Phi(i) = R(i) - n/3*( Phi(1)*R(i-1)*deltat + sum2 + sum4 + Phi(i-1)*R(1)*deltat ); %The last term does not make sense since it is always R(0)!
            if mod(i,100) == 0
                disp(['NumSolv: ',num2str(i/length(t)*100,'%.2f'),' %']);
            end
        end
    end
    
    time = toc(StartNum); 
    beep on; beep; beep off; %Beep at the end of this lengthy computation.
end
