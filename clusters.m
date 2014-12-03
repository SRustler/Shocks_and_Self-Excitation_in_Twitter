%% Initializations
clear; clc; %userpath(pwd) %Changes startup folder to pwd.
Print = [1 0 0]; %Vector of booleans for printing: cluster stats, RMK density estimations, RMK stats

%Initialized variables:
StartSim = tic;
Tmax = 5000000;
kernel = 'multi-exp'; %'single-exp' 'multi-exp' In approximation of a PL-kernel
EXP = 0;
POWSUM = 0;
switch kernel
    case 'single-exp'
        EXP = 1; %Make simple boolean to enter later if-routines most quickly.
        H = hawkes('exp');
        pars_tru = [0.1 0.7 2]; %[mu n tau] where \Phi(t) = 1/tau * exp(-t/tau).
    case 'multi-exp'   %Same as powsum but allowing increasing exps, too.
        POWSUM = 0; %Make simple boolean to enter later if-routines most quickly.
        H = hawkes('powsum'); 
        pars_tru = [0.1 0.7 0.1 0.8]; %[mu n c theta] where \Phi(t) = \Phi(t) = \theta*c^\theta / (t+c)^(1+\theta). Choose lower c(>!0) and higher theta (<!1) for shorter memory.
    otherwise
        error('Warning: No correct kernel selected!');
end

%Simulations:
[tevnts,~,NumGens,gens,cid] = H.fastsim(pars_tru, Tmax); %gens and cid are vectors of length(tevnts). Entries of the former are the generation of the event (e.g. k=0 for immigrants, k=2 for grand-daughters). Entries of the latter are the cluster index of the event (e.g. max(cid) will be the latest occuring cluster independ of all other clusters with a different index.).
[rwt,K,t0] = get_rwt(cid,tevnts,gens); %Calculate renormalized waiting times and retrieve time stamps of immigrants
maxrwtgen = gens(max(rwt)==rwt); %Find highest generation of the most long-lived cluster.
rwt=rwt(rwt>0); %The entries of rwt=0 have to be removed because rwt \in ]0;inf[ and otherwise fit of R(t) would be spoiled.

disp(['Simulation took ', num2str(toc(StartSim)), ' seconds.']);

%True kernel functions:
if EXP
    phi_tru = @(t) exp(-t/pars_tru(3))/pars_tru(3);
    R_tru = @(t) exp(-t/pars_tru(3)*(1-pars_tru(2)))/pars_tru(3); %Analytical solution via Laplace method.
end 
if POWSUM %This can take very long.
    phi_ideal = @(t) pars_tru(4)*pars_tru(3)^pars_tru(4) * (t+pars_tru(3)).^-(1+pars_tru(4)); %Exact form of PL tb approximated (c-shifted Omori). Goodness of approximation can be tested in test_kernels.m.
    phi_tru = @(t) phi_powsum(pars_tru, t); %Approximation of the above form as it was used in Hardiman, Bercot, Bouchaud 2013. Note that we use the subscript 'tru' because this kernel is also used for the simulation!
%     [tt,R_tru,msg] = get_Rpowsum(pars_tru,Tmax/100); %This is the analytical solution of phi_tru obtained via Laplace transformation (more details in function comments). 
%     fprintf(msg); %Display performance of get_Rpowsum().
    % Or load R_tru(tt) instead:
    load('R_tru_mu1n04c4theta09.mat','R_tru','tt');
end 

%Life time of largest cluster:
[h,liftim] = get_liftim(cid,tevnts); %Note that liftim does NOT have to equal the largest response time!

disp(['Initialization completed in ',num2str(toc(StartSim)),' seconds.'])

%% Estimation of renormalized kernel R(t)
StartEstimR = tic;

% Estimate R(t) as histogram:
binwidth = 0.1; 
[Ri,ti] = HISTforR(rwt,binwidth,(length(tevnts)-length(t0))/length(tevnts)); %The third argument is the empirical branching ratio used for the correct scaling of R.

% Estimate R(t) via different kernels for KDE: 
bw1 = 1; bw2 = 0.9; %Bandwidth for KDE. Note: The smaller bandwidth the better the fit at the head but also overfitting at the tail.
[ri1,xi] = ksdensity(rwt,ti,'kernel',@exppdf,'bandwidth',bw1); %Use exponential PDFs as underlying kernel
[ri2,~] = ksdensity(rwt,ti,'kernel',@posnormpdf,'bandwidth',bw1); %Use positive half-Gaussian as underlying kernel
[ri3,~] = ksdensity(rwt,ti,'kernel','normal','bandwidth',bw2); %At rwt-values, i.e. X not "y"-values, normzd kernels are summed up to form a pdf. X will then be points at which f=pdf(X), i.r. ri1=pdf(xi). Note that we take X from the histogram because it does not matter which specific X values we take.
[ri4,xi4] = ksdensity(rwt,[fliplr(-ti) ti],'kernel','normal','bandwidth',bw2); %This turns out to be the best KDE for our data.
    r = fliplr(ri4(xi4 < 0)); %ri-values to be reflected and added to the respective positive ri.
    pos = find(xi4>0); pos = pos(1); %Position of the first positive ri-value.
    ri4(pos:pos+length(r)-1)= ri4(pos:pos+length(r)-1) + r; %Only update the first ri4 entries.
    ri4(1:pos-1) = []; 
    xi4(1:pos-1) = [];
    
% Estimate R(t) via AKDE with normal kernel and post-reflection:  
[R_akde,~,~] = posAKDE(rwt,ti,0.5); %Note that R_akde is normalized to 1.
if EXP %For analytical EXP-case still employ fitting.
    R_tru = R_tru(ti); %Turn function into simple vector such that it is like in the POWSUM case.
    tt = ti;
    fo = fitoptions('Method','NonlinearLeastSquare'); %This is already the default fitting method. Potentially change.
    Rf_akde = fit(ti',transpose(R_akde/(1-((length(tevnts)-length(t0))/length(tevnts)))),'exp1',fo); %Perform fit on best-KDE'd R(t).
    ab = coeffvalues(Rf_akde); %Retrieve coeffs of f(x) = a*exp(b*x) OR f(x) = Sum[a_i*exp(-b_i*x),{i,1,15}]
    conf_ab = confint(Rf_akde); %Retrieve confidence bounds of a and b.
    Rf_akde = @(t) ab(1)*exp(ab(2)*t);  
    Rf_akde_lower = @(t) conf_ab(1,1)*exp(conf_ab(1,2)*t);
    Rf_akde_upper = @(t) conf_ab(2,1)*exp(conf_ab(2,2)*t);
end

disp(['Estimating the renormalized kernel took ', num2str(toc(StartEstimR)), ' seconds.']);

%% Calculation of bare kernel Phi(t)
StartEstimPhi = tic;

%Estimation via MLE and cluster knowledge [=knowing R(t)] (Latter only works for EXP)
[pars_mle,~] = H.estim(tevnts, Tmax); %Parameters estimated via MLE. ???Note: estim relies on gif_powsum which is not yet updated.
if EXP
    phi_mle = @(t) exp(-t/pars_mle(3))/pars_mle(3); %[Result to be plotted]  
    [pars_clu] = clu_estim(length(tevnts),length(t0),Tmax,ab(2)); %Parameters estimated from knowledge of cluster and analytical solution of Laplace trafo
%     phi_clu = exp(-ti./pars_clu(3))/pars_clu(3); %[Result to be plotted]
    [~,phi_clu] = NumSolv_PhiR(ri4,ti,pars_clu(2),'simp'); %[Result to be plotted]  %??? This still takes too long!
%     ff = fit(ti',phi_clu','exp1',fo);
%     aa=coeffvalues(ff);
%     pars_clu1 = pars_clu;
%     [pars_clu] = clu_estim(length(tevnts),length(t0),Tmax,aa(2)); %Parameters estimated from knowledge of cluster and analytical solution of Laplace trafo
end 
if POWSUM
    phi_mle = @(t) pars_mle(4)*pars_mle(3)^pars_mle(4) * (t+pars_mle(3)).^-(1+pars_mle(4)); %[Result to be plotted]
    [alpha,c,~] = plfit(rwt); %Test this!???
    [pars_clu] = [length(tevnts)/Tmax (length(tevnts)-length(t0))/length(tevnts) c alpha-1]; %Note that Clauset's alpha = theta+1.
    %Note: phi_clu(1:3) only exists for the an analytical solution of the single-exp case of R(t). As soon as R(t) is multi-exp, the analytical solution becomes too elaborate (confirmed with gettingPhiFromRasExp.nb where 'InverseLaplaceTransform[\Philap[s, 1, 5, 10, 4, 2, 1, n], s,t]' becomes very ugly. Note that here we used a triple sum with known/fitted coeffs.)
    [~,phi_clu] = NumSolv_PhiR(R_akde,ti,pars_clu(2),'trap'); %[Result to be plotted]  
end

disp(['Calculating the bare kernel took ', num2str(toc(StartEstimPhi)), ' seconds.']);

%% Plots

if Print(1) % Cluster statistics
    fig1 = plot_ClustStats(tevnts,Tmax,pars_tru,length(t0),gens,NumGens,h,liftim);
end

if Print(2) % Comparing different estimations of R(t)
    fig2 = figure; 
        hold on; grid on;

        plot(ti,Ri,'yo','MarkerSize',10); 
        plot(ti,R_akde/(1-pars_clu(2)),'g.');
        plot(xi,ri1/(1-pars_clu(2)),'c-'); %Scale ri1 to correct value: Recall that kappa = 1/(1-n).
        plot(xi,ri2/(1-pars_clu(2)),'b-'); 
        plot(xi,ri3/(1-pars_clu(2)),'r-'); 
        plot(xi4,ri4/(1-pars_clu(2)),'m-'); 

        set(gca,'YScale','log');
        xlabel('t - Waiting time');
        ylabel('log[R(t)] - Renormalized kernel');
        xlim([0 max(rwt)]); 
        ylim([0.0000001 1]); 
        legend(['Histogram with binwidth ',num2str(binwidth)],...
            'Adaptive KDE with normal kernel and post-reflection',...
            ['KDE with exponential kernel and bandwidth ',num2str(bw1)],...
            ['KDE with positive-normal kernel and bandwidth ',num2str(bw1)],...
            ['KDE with normal kernel and bandwidth ',num2str(bw2)],...
            ['KDE with normal kernel and post-reflection and bandwidth ',num2str(bw2)]);

        if POWSUM
            set(gca,'XScale','log'); 
            xlabel('log(t) - Waiting time');
            xlim([ti(1)-ti(1)/2 max(rwt)]);
        end
        hold off;
end
            
if Print(3) % Comparing different estimations of Phi(t)
    fig3 = figure; 
    
        subplot(2,1,1);
            hold on; grid on;
            
            set(gca,'YScale','log'); 
            ylabel('log[R(t)] - Renormalized memory kernel');
%             ylim([0.0001 1]);
            
            if EXP
                
                h3 = plot(ti,Rf_akde_upper(ti),'Color',[0.5 0.5 0.5]); %Grey bounds
                h4 = plot(ti,Rf_akde_lower(ti),'LineStyle','None'); %Grey bounds
                    y1=Rf_akde_lower(ti);                      %Create first curve
                    y2=Rf_akde_upper(ti);                      %Create second curve
                    xx=[ti,fliplr(ti)];                   %Create continuous xx value array for plotting
                    yy=[y1,fliplr(y2)];                 %Create yy-values for out and then back
                    fill(xx,yy,[0.5 0.5 0.5]);%Fill with grey
                h2 = plot(ti,Rf_akde(ti),'k--');
                h1 = plot(ti,R_akde/(1-pars_clu(2)),'bo','MarkerSize',5); %Plot best KDE-method from fig2, i.e. "AKDE, normal with post-reflection"
                h5 = plot(tt,R_tru,'g-'); %tt will different and shorter than X as it is calculated in INVLAP() in get_Rpowsum(). Matching the range of X will dramatically increase the computation time!
                
                xlim([0 max(rwt)]);
                legend([h1 h2 h3 h5],{...
                    'Renormalized kernel R(t) estimated via AKDE',...
                    ['Fitted R(t) with characteristic decay time ',num2str(-1/ab(2)),' = \tau/(1-n)'],...
                    '0.95-confidence bounds',... %Label for h3 and h4 combined 
                    ['Theoretical R(t) = 1/\tau_{tru}*exp(-(1-n_{tru})/\tau_{tru}*t), where \tau_{tru}/(1-n_{tru}) = ',num2str(pars_tru(3)/(1-pars_tru(2)))]});          
            end
            if POWSUM
                h1 = plot(ti,R_akde/(1-pars_clu(2)),'bo','MarkerSize',5); %Plot best KDE-method from fig2, i.e. "AKDE, normal with post-reflection"
                h5 = plot(tt,R_tru,'g-'); %tt will different and shorter than X as it is calculated in INVLAP() in get_Rpowsum(). Matching the range of X will dramatically increase the computation time!
                
                set(gca,'XScale','log'); hold off;
                xlim([ti(1)-ti(1)/2 max(rwt)]);
                legend([h1 h5],{...
                    'Renormalized kernel R(t) estimated via AKDE',...
                    ['Theoretical R(t|pars) with pars = [\mu n c \theta] = [',num2str(pars_tru),']']});   %???Analytical expression?       
            end
            
        subplot(2,1,2);
            hold on; grid on;
            
            plot(ti,R_akde/(1-pars_clu(2)),'Color',[0.5  0.5  1],'Marker','.','LineStyle','None');
%             plot((ti(1):length(R_tru))/length(R_tru)*ti(end),R_tru,'g--');
            plot(tt,R_tru,'g--');
            plot(ti,phi_tru(ti),'g-','LineWidth',3); %True Phi(t) used to generate the data in the first place.
            plot(ti,phi_mle(ti),'c-','LineWidth',2);
            plot(ti,phi_clu,'mo','MarkerSize',5); %Phi(t) obtained from cluster knowledge (For EXP analytically, for POWSUM via NumSolv_PhiR)
            
            set(gca,'YScale','log'); hold off;
            xlabel('t - Waiting time');
            ylabel('log[\Phi(t)] - Bare memory kernel');
            xlim([0 max(rwt)]);
%             ylim([0.0001 1]);
            legend(...
                'Renormalized kernel R(t) estimated via AKDE',...
                'R_{tru}(t) - Theoretical renormalized memory kernel',...
                ['\Phi_{tru}(t|\tau_{tru} = ',num2str(pars_tru(3)),') - True generating kernel'],...
                ['\Phi(t|\tau_{MLE} = ',num2str(pars_mle(3)),') - Obtained via MLE'],...
                ['\Phi_{clu}(t|\tau_{clu} = ',num2str(pars_clu(3)),') - Obtained numerically from R(t)']);
            
            if POWSUM
                set(gca,'XScale','log'); hold off;
                xlabel('log(t) - Waiting time'); 
                xlim([ti(1)-ti(1)/2 max(rwt)]);
                legend(...
                    'Renormalized kernel R(t) estimated via AKDE',...
                    'R_{tru}(t) - Theoretical renormalized memory kernel',...
                    ['\Phi_{tru}(t|\theta_{tru} = ',num2str(pars_tru(4)),') - True generating kernel'],...
                    ['\Phi(t|\theta_{MLE} = ',num2str(pars_mle(4)),') - Obtained via MLE'],...
                    ['\Phi_{clu}(t|\theta_{clu} = ',num2str(pars_clu(4)),') - Obtained numerically from R(t)']);
            end
end

disp('PARS: true       MLE   cluster'); %Print little table of parameters for quick comparison
disp(horzcat(pars_tru',pars_mle',pars_clu')); 
disp(['Total time elapsed: ',num2str(toc(StartSim)),' seconds.'])
disp(['The most long-lived cluster lived ',num2str(max(rwt)),' time units and had ',num2str(maxrwtgen),' generations.']) %The k-value here should be \leq max(gens)
disp(['The highest-generation cluster lived ',num2str(liftim),' time units and had ',num2str(max(gens)),' generations.'])
