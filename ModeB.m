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
    %distance between UEs
    distance = 0.02;    % km
    eta = 0.7;
    theta =0.99;
    %path Loss for Mode B
%     P_L_B1 =  pathLossD2DLink(eta,distance); 
%     P_L_B2 = pathLossD2DLink(1-eta,distance);
%     SNR_coef = signalToNoiseRatioCoefficient();
      
        
%-----------------------------------simulation Mode A----------------------------
   %brake interval  of transmission power range (0,0.25] into m itnervals 
    m = 200;
    n = 10000;
    % generate 100000 samples
    % channel coefficient, exponentially distributed 
    g10 = - log(1 - rand(1,n))/2;
    g20 = - log(1 - rand(1,n))/2;
   % g10 = randn(1,n) + 1i * randn(1,n);
   % g20 = randn(1,n) + 1i * randn(1,n);
   % g10 = abs(sqrt(1/2) * g10).^2;
   % g20 = abs(sqrt(1/2) * g20).^2;

% build vectors of EE and SE for simulation and analytical cases
    for i = 1 : m  
        P = P_max_D * i / m;
        average_eeB_sim(i) = mean(energyEfficiencyModeB(eta, g10, g20, P))/10^6;
        average_seB_sim(i) = mean(spectralEfficiencyModeB(eta, g10, g20, P));
          for k = 1:4
              average_eeB_apr_N(k,i) = energyEfficiencyApproximateModeB(eta,theta, P,k*2-1)/10^6;
              average_seB_apr_N(k,i) = spectralEfficiencyApproximateModeB(eta, theta, P, k*2-1);
          end
        
     end
% ------------------------ Energy efficiency subplot-------------------------
  
    P = linspace( P_max_D/m, P_max_D, m);
    figure;subplot(2,2,1); 
    grid on;
    plot(P, average_eeB_sim,'b-');
    hold on;
    for k = 1:4
         plot(P, average_eeB_apr_N(k,:),'LineStyle','--','color',rand(1,3));
    end
    
    %finds optimal transmission power P* for simulation EE
    [M I] = max(average_eeB_sim);
     P_star = P_max_D*I/m
     mean(energyEfficiencyModeB(eta, g10, g20, P_star))/10^6
     mean(spectralEfficiencyModeB(eta, g10, g20, P_star))
     line([P_star P_star], [M-50 M + 50],'Color', 'blue','LineStyle','-');
    
    % finds optimal transmission power P* for analytical EE
    xB_star = bisectionMethod(@fB, 0.000001, 100, 0.00000001);
    P_star = 1 / xB_star
    line([P_star P_star], [M-50 M + 50],'Color', 'blue','LineStyle','--');
    hold off;
    legend('Multihop-D2D : Simulation','Multihop-D2D : Analytical N=1','Multihop-D2D : Analytical N=3','Multihop-D2D : Analytical N=5','Multihop-D2D : Analytical N=7');
    xlabel('UE Transmittion Power (W)');
    axis([0 0.25 180 330]);
    ylabel('Energy Efficiency (Mbits/Joule)');
    title('Average energy efficiency') 

% ------------------------ Spectral efficiency subplot-------------------------
   
     subplot(2,2,2);
     plot(P, average_seB_sim,'b-');
     hold on;
     for k = 1:4
         plot(P, average_seB_apr_N(k,:),'LineStyle','--','color',rand(1,3));
     end    
     hold off;
     legend('Multihop-D2D : Simulation','Multihop-D2D : Analytical N=1','Multihop-D2D : Analytical N=3','Multihop-D2D : Analytical N=5','Multihop-D2D : Analytical N=7');
     xlabel('UE Transmittion Power (W)');
     ylabel('Spectral Efficiency (bits/s/Hz)');
     title('Average spectral efficiency');

% ------------------------ Average SNR subplot-------------------------

    subplot(2,2,3);
    plot(P,  pow2db(mean(g10) * signalToNoiseRatioModeB1(theta,eta,P) + mean(g20) * signalToNoiseRatioModeB2(theta,eta,P)) ) ;
    xlabel('UE Transmittion Power (W)');
    ylabel('Average Received SNR (dB)');
    title('Signal to noise ratio');

% ------------------------ Energy efficiency vs Spectral efficiency subplot-------------------------

    subplot(2,2,4);
    plot(average_seB_sim, average_eeB_sim,'-');
    hold on;
    for k = 1:4
         plot(average_seB_apr_N(k,:), average_eeB_apr_N(k,:),'LineStyle','--','color',rand(1,3));
    end    
    hold off;
    legend('Multihop-D2D : Simulation','Multihop-D2D : Analytical N=1','Multihop-D2D : Analytical N=3','Multihop-D2D : Analytical N=5','Multihop-D2D : Analytical N=7');
    xlabel('Spectral Efficiency (bits/s/Hz)');
    ylabel('Energy Efficiency (Mbits/Joule)');
    title('Energy efficiency versus spectral efficiency');
      
% ----------------------------------- Plot EE(eta), SE(eta) -------------
       
    for j = 0 : 9  
        eta_ = 0.5 + j * 0.05;
        for i = 1 : m  
            P = P_max_D * i / m;
            eeB_sim(i) = mean(energyEfficiencyModeB(eta_, g10, g20, P))/10^6;
            seB_sim(i) = mean(spectralEfficiencyModeB(eta_, g10, g20, P));
            eeB_apr(i) = energyEfficiencyApproximateModeB(eta_,theta, P,7)/10^6;
            seB_apr(i) = spectralEfficiencyApproximateModeB(eta_, theta, P, 7);        
        end
%         find optimal transmission power
        [ee_eta_sim(j+1) P_star_sim]  = max (eeB_sim);
        [ee_eta_apr(j+1) P_star_apr] = max (eeB_apr);
        se_eta_sim(j+1) = seB_sim(P_star_sim);
        se_eta_apr(j+1) = seB_apr(P_star_apr);
    end
     eta_ = 0.5:0.05:0.95; 
     
     figure; plot(eta_, ee_eta_sim,'-', eta_, ee_eta_apr,'--o');
     legend('Multihop-D2D : Simulation','Multihop-D2D : Analytical');
     xlabel('\eta');
     ylabel('Energy Efficiency (Mbits/Joule)');
     title('Energy efficiency versus \eta');
     axis([0.5 0.95 250 330]);
     
     figure; plot(eta_, se_eta_sim,'-', eta_, se_eta_apr,'--o');
     legend('Multihop-D2D : Simulation','Multihop-D2D : Analytical');
     xlabel('\eta');
     ylabel('Spectral Efficiency (bits/s/Hz)');
     title('Spectral efficiency versus \eta');
     axis([0.5 0.95 8.5 11.5]);

% ---------------------------FUNCTIONS----------------------------------------    

    % SNR at UE1 
    function snrB1 = signalToNoiseRatioModeB1(theta, eta, P_D1)
        P_L_B1 =  pathLossD2DLink(eta,distance); 
        snrB1 = theta * P_D1 / (3 * W * N0 * P_L_B1);
    end

    % SNR at UE2 
    function snrB2 = signalToNoiseRatioModeB2(theta, eta, P_D1)
        P_L_B2 = pathLossD2DLink(1-eta,distance);
        b2 = ( (eta/(1-eta))^4 + 2 ) * W * N0 * P_L_B2  / theta;
        snrB2 =  P_D1 / b2;

    end

   
    % return the path loss for D2D  link regarding to the distance d. 
    function pathLoss = pathLossD2DLink(eta_, d)           
        pathLoss = (10^14.8) * (d * eta_)^4 * noiseFigure; 
    end

    % returns vector of n elements with EE for each h0 with given transmission power P
    function eeB = energyEfficiencyModeB(eta, g10, g20, P_D1)
        P_D2 = P_D1 * ((1-eta) /eta) ^ 4;
        P_D3 = P_D1;
        P_L_B1 =  pathLossD2DLink(eta,distance) ; 
        P_L_B2 = pathLossD2DLink(1-eta,distance);
        mu_B1 = P_D3 /( W * N0 * P_L_B1); 
        v_B1 = (P_D3 + P_D1) / P_D2 * ((1-eta)/eta) ^ 4;
        r_B1 = (W * log2(1 + (mu_B1 .* g10 .* g20)/(v_B1 .* g10 + g20) )  / 2);
        mu_B2 = P_D3 /( W * N0 * P_L_B2 ); 
        v_B2 = (P_D3 + P_D2) / P_D1 * (eta/(1-eta)) ^ 4;
        r_B2 = W * log2(1 + (mu_B2 .* g10 .* g20)/(v_B2 .* g20 + g10) )  / 2;
        P_total = P_D1 + P_D2 + P_D3 + 3 * P_c_D;
        eeB = (r_B1 + r_B2)/P_total;
    end
     % returns vector of n elements with SE for each |g10|^2, |g20|^2 with given transmission power P_D, and eta 
    function seB = spectralEfficiencyModeB(eta, g10, g20, P_D1)
        P_D2 = P_D1 * ((1-eta) /eta) ^ 4;
        P_D3 = P_D1;
        P_L_B1 =  pathLossD2DLink(eta,distance); 
        P_L_B2 = pathLossD2DLink(1-eta,distance);
        mu_B1 = P_D3 /( W * N0 * P_L_B1 ); 
        v_B1 = (P_D3 + P_D1) / (P_D2) * ((1-eta)/eta) ^ 4;
        r_B1 = W * log2(1 + (mu_B1 .* g10 .* g20)/(v_B1 .* g10 + g20) )  / 2;
        mu_B2 = P_D3 /( W * N0 * P_L_B2 ); 
        v_B2 = (P_D3 + P_D2) / (P_D1) * (eta/(1-eta)) ^ 4;
        r_B2 = W * log2(1 + (mu_B2 .* g10 .* g20)/(v_B2 .*g20 + g10) )  / 2;
        seB =(r_B1 + r_B2)/ W;
    end

    %returns approximate value of  SE for given eta, theta, transmission power P and
    %number of Taylor series members N used in approximation
    function seB = spectralEfficiencyApproximateModeB(eta, theta, P_D1, N)
        SNR_B1  = signalToNoiseRatioModeB1(theta,eta,P_D1);
        TN_B1 = T_N(equation_solver(N),N) - T_N(1/SNR_B1,N);
        SNR_B2  = signalToNoiseRatioModeB2(theta,eta,P_D1);
        TN_B2 = T_N(equation_solver(N),N) - T_N(1/SNR_B2,N);
        seB =  (exp(1/SNR_B1) * TN_B1 + exp(1/SNR_B2) * TN_B2 )/ ( 2 * log(2));
    end

    %returns approximate value of  EE for given transmission power P and
    %number of Taylor series members N used in approximation
    function eeB = energyEfficiencyApproximateModeB(eta,theta, P_D1, N)
        P_D2 = P_D1 * ((1-eta) /eta) ^ 4;
        P_D3 = P_D1;
        P_total = P_D1 +P_D2 +P_D3 + 3 * P_c_D;
        SNR_B1  = signalToNoiseRatioModeB1(theta,eta,P_D1);
        TN_B1 = T_N(equation_solver(N),N) - T_N(1/SNR_B1,N);
        SNR_B2  = signalToNoiseRatioModeB2(theta,eta,P_D1);
        TN_B2 = T_N(equation_solver(N),N) - T_N(1/SNR_B2,N);
        eeB =  W * (exp(1/SNR_B1) * TN_B1 + exp(1/SNR_B2) * TN_B2 )/ ( 2 * log(2) * P_total);
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

    %first derivative of the EE function 
    function yB = fB(eta, x)
        P_L_B1 = pathLossD2DLink(eta,distance); 
        P_L_B2 = pathLossD2DLink(1-eta,distance);
        b1 = 3 * W *N0 * P_L_B1  / theta; 
        b2 = ( (eta/(1-eta))^4 + 2 ) * W * N0 * P_L_B2  / theta;
        b3 = ( (1-eta)/eta)^4 + 2;
        yB =   (psi(x,b1) + psi(x,b2)) .* (b3 ./ x + 3 * P_c_D) + (phi(x,b1) + phi(x,b2)) * b3 ./ x.^2;
    end
    
    %auxilary function to find optimal transmission power 
    function phi = phi(x,b) 
        phi = exp(b .* x) * (b .* x -1 - log(b .* x));
    end

    %auxilary function to find optimal transmission power
    function psi = psi(x,b)
        psi = b * phi(x,b) + exp(b .* x) * (b - 1. / x);
    end 

    %finds sum of first n members of the Taylor Series 
    function taylorS = taylorSeries(t,n)
        taylorS = 1;
        for i = 1:n
            taylorS = taylorS + power(-t,i)/factorial(i) ;
        end
    end

    % n factorial, for integer n
    function fact = factorial(n)
        if ( n == 0 ) 
            fact = 1;
        else
            fact = prod(1:n);
        end
    end

    % Repeatedly bisects an interval and then selects a subinterval in which a root must lie for further processing.
    % Input: f – function , [a,b] -interval,  epsilon – accuracy.
    % Output:  estimation of root with accuracy epsilon.

    function root = bisectionMethod(fB, a, b, epsilon)
        % Check that that neither end-point is a root
        % and if f(a) and f(b) have the same sign, throw an exception.
        if ( fB(eta, a) == 0 )
            root = a;
            return;
        elseif ( fB(eta, b) == 0 )
            root = b;
            return;
        elseif ( fB(eta, a) * fB(eta, b) > 0 )
            error( 'f(a) and f(b) do not have opposite signs' );
        end

        % We will iterate until |a-b| <= epsilon)
        while (abs(a-b) > epsilon)
        % Find the mid-point
        c = (a + b)/2;
                   % Check if we found a root or whether or not
            % we should continue with:
            %          [a, c] if f(a) and f(c) have opposite signs, or
            %          [c, b] if f(c) and f(b) have opposite signs.
            if ( fB(eta,c) == 0 )
                root = c;
                return;
            elseif ( fB(eta, c)*fB(eta,a) < 0 )
                b = c;
            else
                a = c;
            end
        end
        root = (b + a)/2;
    end
end
