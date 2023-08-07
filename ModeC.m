function main() 
    % Table 1
    % Main Simulation Parameters 
    % transmission Power of the BS 
    P_BS = 30;                              % W
    % circuit Power of the BS 
    P_c_BS = 10;                            % W
    % max Transmission Power UEs 
    P_max_D = 0.25;                        % W 
    % circuit Power of UEs
    P_c_D = 0.1;                            % W
    % channel bandwidth
    W = 1e7 ;                          % Hz 
    noiseFigure = 10^0.7;                    % dB
    % Thermal Noise Density 
    N0 =  10^(-20.4) ;                 %W/Hz = dBm/Hz
    %distance between UE1 and BS
    d1 = 0.1;    % km
    %distance between UE2 and BS
    d2 = 0.2;    % km
  
    alpha = 3.76;
    theta = 0.9;
    %path Loss for Mode B
%     P_L_B1 =  pathLossD2DLink(eta,distance); 
%     P_L_B2 = pathLossD2DLink(1-eta,distance);
%     SNR_coef = signalToNoiseRatioCoefficient();
      
        
%-----------------------------------simulation Mode A----------------------------
   %brake interval  of transmission power range (0,0.25] into m itnervals 
    m = 200;
    n = 100000;
    % generate 100000 samples
    % channel coefficient, exponentially distributed 
     f10 = - log(1 - rand(1,n))/2;
     f20 = - log(1 - rand(1,n))/2;
%      f10 = normrnd(0,1,1,n) + 1i * normrnd(0,1,1,n);
%      f10 = abs(sqrt(1/2) *f10).^2;
%      f20 = normrnd(0,1,1,n) + 1i * normrnd(0,1,1,n);
%      f20 = abs(sqrt(1/2) *f20).^2;

% build vectors of EE and SE for simulation and analytical cases
    for i = 1 : m  
        P = P_max_D * i / m;
        average_eeC_sim(i) = mean(energyEfficiencyModeC(f10, f20, P,d1,d2))/10^6;
        average_seC_sim(i) = mean(spectralEfficiencyModeC(f10, f20, P,d1,d2));
          for k = 1:4
              average_eeC_apr_N(k,i) = energyEfficiencyApproximateModeC(theta, P,k*2-1,d1,d2)/10^6;
              average_seC_apr_N(k,i) = spectralEfficiencyApproximateModeC(theta, P, k*2-1,d1,d2);
          end
        
     end
% ------------------------ Energy efficiency subplot-------------------------
  
    P = linspace( P_max_D/m, P_max_D, m);
    figure;subplot(2,2,1); 
    plot(P, average_eeC_sim,'b-');
    hold on;
    for k = 1:4
         plot(P, average_eeC_apr_N(k,:),'LineStyle','--','color',rand(1,3));
    end
    hold off;  
    legend('Through BS : Simulation','Through BS : Analytical N=1','Through BS : Analytical N=3','Through BS : Analytical N=5','Through BS : Analytical N=7');
    xlabel('UE Transmittion Power (W)');
%     axis([0 0.25 180 360]);
    ylabel('Energy Efficiency (Mbits/Joule)');
    title('Average energy efficiency') 

% ------------------------ Spectral efficiency subplot-------------------------
   
     subplot(2,2,2);
     plot(P, average_seC_sim,'b-');
     hold on;
     for k = 1:4
         plot(P, average_seC_apr_N(k,:),'LineStyle','--','color',rand(1,3));
     end    
     hold off;
     legend('Through BS : Simulation','Through BS : Analytical N=1','Through BS : Analytical N=3','Through BS : Analytical N=5','Through BS : Analytical N=7');
     xlabel('UE Transmittion Power (W)');
     ylabel('Spectral Efficiency (bits/s/Hz)');
     title('Average spectral efficiency');

% ------------------------ Average SNR subplot-------------------------

    subplot(2,2,3);
    plot(P,  pow2db( signalToNoiseRatioModeC1(theta, P, d1,d2)) ,P,   pow2db( signalToNoiseRatioModeC2(theta,P,d1,d2))) ;
    legend('Through BS : SNR-C1','Through BS : SNR-C2');
    xlabel('UE Transmittion Power (W)');
    ylabel('Average Received SNR (dB)');
    title('Signal to noise ratio');

% ------------------------ Energy efficiency vs Spectral efficiency subplot-------------------------

    subplot(2,2,4);
    plot(average_seC_sim, average_eeC_sim,'-');
    hold on;
    for k = 1:4
         plot(average_seC_apr_N(k,:), average_eeC_apr_N(k,:),'LineStyle','--','color',rand(1,3));
    end    
    hold off;
    legend('Through BS : Simulation','Through BS : Analytical N=1','Through BS : Analytical N=3','Through BS : Analytical N=5','Through BS : Analytical N=7');
    xlabel('Spectral Efficiency (bits/s/Hz)');
    ylabel('Energy Efficiency (Mbits/Joule)');
    title('Energy efficiency versus spectral efficiency');
    
   
% ---------------------------FUNCTIONS----------------------------------------    

    % SNR at UE1 
    function snrC1 = signalToNoiseRatioModeC1(theta, P_D1,d1,d2)
         v_C1 =  (d2/d1)^alpha .* (P_BS + P_D1) ./ P_D2(P_D1,d1,d2);
         snrC1 = theta * mu_C(d1) ./ (1 + v_C1);
    end

    % SNR at UE2 
    function snrC2 = signalToNoiseRatioModeC2(theta, P_D1,d1,d2)
         v_C2 =  (d1/d2)^alpha .* (P_BS + P_D2(P_D1,d1,d2)) ./ P_D1;
         snrC2 = theta * mu_C(d2) ./ (1 + v_C2);

    end
    % return the path loss for Cellular  link regarding to the distance d. 
    function pathLoss = pathLossCellularLink(d)   
         pathLoss =10^12.81 * d^3.76 * noiseFigure; 
    end
    %     returns the mu_C value regarding to the distance d
    function mu = mu_C(d) 
        P_L_C =  pathLossCellularLink(d) ; 
        mu = P_BS /( W * N0 * P_L_C);
    end
    %   transmission power at the UE2. 
    function pd2 = P_D2(P_D1, d1,d2)
        pd2 = P_D1 .* (d2/d1)^alpha;
        pd2 = min(P_max_D, pd2); 
    end

    %transmission rate at UE1 (simulation)
    function rC1 = transmissionRateC1(f10,f20, P_D1, d1,d2)
        v_C1 =  (d2/d1)^alpha * (P_BS + P_D1) ./ P_D2(P_D1,d1,d2);
        rC1 = W / 2 * log2 ( 1 + mu_C(d1) * f10 .* f20 / (v_C1*f10 + f20));
    end

    %transmission rate at UE2 (simulation)
    function rC2 = transmissionRateC2(f10,f20, P_D1, d1,d2)
        v_C2 =  (d1/d2)^alpha * (P_BS + P_D2(P_D1,d1,d2)) ./ P_D1;
        rC2 = W / 2 * log2 ( 1 + mu_C(d2) .* f10 .*f20 / (v_C2*f10 + f20));
        
     end
    
    % returns vector of n elements of EE for each f10, f20 with given
    % transmission power P_D_1
    function eeC = energyEfficiencyModeC(f10, f20, P_D1,d1,d2)
        P_total_C =  P_BS + P_D1 + P_D2(P_D1,d1,d2) + P_c_BS + 2 * P_c_D; 
        eeC = (transmissionRateC1(f10,f20, P_D1, d1,d2) + (transmissionRateC2(f10,f20, P_D1, d1,d2)))/ P_total_C;
          
    end
     % returns vector of n elements of SE for each f10, f20 with given transmission power P_D1 
    function seC = spectralEfficiencyModeC(f10, f20, P_D1, d1,d2)
        seC = (transmissionRateC1(f10,f20, P_D1, d1,d2) +  transmissionRateC2(f10,f20, P_D1, d1,d2)) / W; 
    end

    %returns approximate value of  SE for given theta, UE1 transmission power  P_D1, dstance  d1(d2) from UE1(UE2) to BS,  and
    %number of Taylor series members N used in approximation
    function seC = spectralEfficiencyApproximateModeC( theta, P_D1, N, d1,d2)
        SNR_C1  = signalToNoiseRatioModeC1(theta, P_D1, d1,d2);
        TN_C1 = T_N(equation_solver(N),N) - T_N(1/SNR_C1,N);
        SNR_C2  = signalToNoiseRatioModeC2(theta, P_D1, d1,d2);
        TN_C2 = T_N(equation_solver(N),N) - T_N(1/SNR_C2,N);
        seC =  (exp(1/SNR_C1) * TN_C1 + exp(1/SNR_C2) * TN_C2 )/ ( 2 * log(2));
    end

    %returns approximate value of  EE for given theta, UE1 transmission power  P_D1, dstance  d1(d2) from UE1(UE2) to BS,  and
    %number of Taylor series members N used in approximation
    function eeC = energyEfficiencyApproximateModeC(theta, P_D1, N, d1,d2)
        P_total_C =  P_BS + P_D1 + P_D2(P_D1,d1,d2) + P_c_BS + 2 * P_c_D;
        SNR_C1  = signalToNoiseRatioModeC1(theta, P_D1, d1,d2);
        TN_C1 = T_N(equation_solver(N),N) - T_N(1/SNR_C1,N);
        SNR_C2  = signalToNoiseRatioModeC2(theta, P_D1, d1,d2);
        TN_C2 = T_N(equation_solver(N),N) - T_N(1/SNR_C2,N);
        eeC =  W * (exp(1/SNR_C1) * TN_C1 + exp(1/SNR_C2) * TN_C2 )/ ( 2 * log(2) * P_total_C);
    end

    %  equation_solver(1)
    % finds the real root of the equation sum(first N members of Taylor series) = 0.
    function root = equation_solver(N)
        eq = [];
            for n = 0:N
              eq(n+1) = (-1)^(N-n) / factorial(N-n);
            end
        eq = roots(eq);
        root = eq(real(eq)>0&imag(eq)==0);
    end

    % finds the value of function T_N(t) in the point t.
    function TN = T_N(t,N)
        value = 0;
        for n = 1:N
           value = value + power(-1,n) * t^n / (n * factorial(n)); 
        end
        TN = value + log(t);
    end

    % n factorial, for integer n
    function fact = factorial(n)
        if ( n == 0 ) 
            fact = 1;
        else
            fact = prod(1:n);
        end
    end

  
end
