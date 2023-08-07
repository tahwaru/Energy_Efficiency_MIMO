function main() 
    % Table 1
    % Main Simulation Parameters 
    % transmission Power of the BS 
    clc
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
    distance = 0.02;                     % km    
    %path Loss for Mode A
    eta = 0.5;
    alpha = 3.76;
    theta = 0.999;
    thetaC = 0.8;
    epsilon = 0.0000001;
            
%-----------------------------------simulation Mode A----------------------------
    %brake interval  of transmission power range (0,0.25] into m itnervals
    m = 200;
    n = 100000;
    % generate 100000 samples 
%   channel coefficient, exponentially distributed 
    %h0 = exprnd(1,1,n);
    h0 = randn(1,n) + 1i * randn(1,n);
    h0 = abs(sqrt(1/2) * h0).^2;
    g10 = - log(1 - rand(1,n))/2;
    g20 = - log(1 - rand(1,n))/2;
    f10 = - log(1 - rand(1,n))/2;
    f20 = - log(1 - rand(1,n))/2;
  
%  Mode A,B,C build vectors of EE and SE for simulation and analytical cases
d1=0.3;
d2 = 0.3;

    for d = 1:10
        %mode A
        distance = d/100;
        P_L_A =  pathLossD2DLink(distance); 
        xA_star = bisectionMethodA(@fA, epsilon, 0.2, epsilon);
        P_A_star = W * N0 * P_L_A* noiseFigure / xA_star; 
        average_eeA_sim(d) = mean(energyEfficiencyModeA(h0,P_A_star))/10^6;
        average_seA_sim(d) = mean(spectralEfficiencyModeA(h0,P_A_star));
        average_eeA_apr_N(d) = energyEfficiencyApproximateModeA(P_A_star,7)/10^6;
        average_seA_apr_N(d) = spectralEfficiencyApproximateModeA(P_A_star,7);
        % mode B
        xB_star = bisectionMethodB(@fB, epsilon, 100, epsilon);
        P_B_star = 1 / xB_star; 
        average_eeB_sim(d) = mean(energyEfficiencyModeB(eta, g10, g20, P_B_star))/10^6;
        average_seB_sim(d) = mean(spectralEfficiencyModeB(eta, g10, g20, P_B_star));
        average_eeB_apr_N(d) = energyEfficiencyApproximateModeB(eta,theta, P_B_star,7)/10^6;
        average_seB_apr_N(d) = spectralEfficiencyApproximateModeB(eta, theta, P_B_star,7);
        %mode C
        P_C_star = P_max_D;
        average_eeC_sim(d) = mean(energyEfficiencyModeC(f10, f20, P_C_star,d1,d2))/10^6;
        average_seC_sim(d) = mean(spectralEfficiencyModeC(f10, f20, P_C_star,d1,d2));
        average_eeC_apr_N(d) = energyEfficiencyApproximateModeC(thetaC, P_C_star,7,d1,d2)/10^6;
        average_seC_apr_N(d) = spectralEfficiencyApproximateModeC(thetaC, P_C_star,7,d1,d2);
   end
   comparisonPlot(10,100,10);
   clearVar();
   distance = 0.08;   
   for d = 1:8
        %mode A
        d1 = d *0.05;
        d2 = d1;
        P_L_A =  pathLossD2DLink(distance); 
        xA_star = bisectionMethodA(@fA, epsilon, 0.2, epsilon);
        P_A_star = W * N0 * P_L_A* noiseFigure / xA_star; 
        average_eeA_sim(d) = mean(energyEfficiencyModeA(h0,P_A_star))/10^6;
        average_seA_sim(d) = mean(spectralEfficiencyModeA(h0,P_A_star));
        average_eeA_apr_N(d) = energyEfficiencyApproximateModeA(P_A_star,7)/10^6;
        average_seA_apr_N(d) = spectralEfficiencyApproximateModeA(P_A_star,7);
        % mode B
        xB_star = bisectionMethodB(@fB, epsilon, 100, epsilon);
        P_B_star = 1 / xB_star; 
        average_eeB_sim(d) = mean(energyEfficiencyModeB(eta, g10, g20, P_B_star))/10^6;
        average_seB_sim(d) = mean(spectralEfficiencyModeB(eta, g10, g20, P_B_star));
        average_eeB_apr_N(d) = energyEfficiencyApproximateModeB(eta,theta, P_B_star,7)/10^6;
        average_seB_apr_N(d) = spectralEfficiencyApproximateModeB(eta, theta, P_B_star,7);
        %mode C
        P_C_star = P_max_D;
        average_eeC_sim(d) = mean(energyEfficiencyModeC(f10, f20, P_C_star,d1,d2))/10^6;
        average_seC_sim(d) = mean(spectralEfficiencyModeC(f10, f20, P_C_star,d1,d2));
        average_eeC_apr_N(d) = energyEfficiencyApproximateModeC(thetaC, P_C_star,7,d1,d2)/10^6;
        average_seC_apr_N(d) = spectralEfficiencyApproximateModeC(thetaC, P_C_star,7,d1,d2);
   
    end
   comparisonPlot(50,400,8);

     
% ------------------------ Energy efficiency subplot-------------------------
   

    
% ---------------------------FUNCTIONS----------------------------------------    
    function clearVar()
        average_eeA_sim = [];
        average_seA_sim =[];
        average_eeA_apr_N =  [];
        average_seA_apr_N = []; 
        average_eeB_sim = [];
        average_seB_sim = [];
        average_eeB_apr_N =  [];
        average_seB_apr_N = []; 
        average_eeC_sim = [];
        average_seC_sim = [];
        average_eeC_apr_N =  [];
        average_seC_apr_N = []; 
    end
    function comparisonPlot(start_point, end_point, points)
         d = linspace( start_point, end_point,points);

        figure;subplot(1,2,1); 
        hold on;
        grid on;
        set(gca, 'YScale', 'log')
        plot(d, average_eeA_sim,'b-');
        plot(d, average_eeA_apr_N,'bo--');
        plot(d, average_eeB_sim,'r-');
        plot(d, average_eeB_apr_N,'ro--');
        plot(d, average_eeC_sim,'g-');
        plot(d, average_eeC_apr_N,'go--');
        hold off;
        legend('Direct-D2D : Simulation','Direct-D2D : Analytical N=7','Multihop-D2D : Simulation','Multihop-D2D : Analytical N=7','Through-BS Simulaion', 'Through-BS : Analytical N=7');
        xlabel('D2D distance between UE1 and UE2');
        ylabel('Energy Efficiency (Mbits/Joule)');
        % ------------------------ Spectral efficiency subplot------------------------%
         subplot(1,2,2);
          plot(d, average_seA_sim,'b-');
        hold on;
        plot(d, average_seA_apr_N,'bo--');
        plot(d, average_seB_sim,'r-');
        plot(d, average_seB_apr_N,'ro--');
        plot(d, average_seC_sim,'g-');
        plot(d, average_seC_apr_N,'go--');
        hold off;
        legend('Direct-D2D : Simulation','Direct-D2D : Analytical N=7','Multihop-D2D : Simulation','Multihop-D2D : Analytical N=7','Through-BS Simulaion', 'Through-BS : Analytical N=7');
         xlabel('D2D distance between UE1 and UE2');
         ylabel('Spectral Efficiency (bits/s/Hz)');
         title('Average spectral efficiency');
    end
    % SNR constant  
    function snrA = signalToNoiseRatioModeA(P)
         snrA = P/ ( W * N0 * P_L_A * noiseFigure);
    end

    % return the path loss for D2D  link regarding to the distance d. 
    function pathLoss = pathLossD2DLink(d)           
        pathLoss = (10^14.8) * d^4; 
    end
    function pathLoss = pathLossD2DLinkB(eta_, d)           
        pathLoss = (10^14.8) * (d * eta_)^4 *noiseFigure; 
    end
     % SNR at UE1 
    function snrB1 = signalToNoiseRatioModeB1(theta, eta, P_D1)
        P_L_B1 =  pathLossD2DLinkB(eta,distance); 
        snrB1 = theta * P_D1 / (3 * W * N0 * P_L_B1);
    end

    % SNR at UE2 
    function snrB2 = signalToNoiseRatioModeB2(theta, eta, P_D1)
        P_L_B2 = pathLossD2DLinkB(1-eta,distance);
        b2 = ( (eta/(1-eta))^4 + 2 ) * W * N0 * P_L_B2  / theta;
        snrB2 =  P_D1 / b2;

    end
  
    % returns vector of n elements with SE for each h0 with given transmission power P_D 
    function seA = spectralEfficiencyModeA(h,P_D)
        SNR = signalToNoiseRatioModeA(P_D);
        seA = log2(1 +  SNR .* h );
    end

    %returns approximate value of  SE for given transmission power P and
    %number of Taylor series members N used in approximation
    function seA = spectralEfficiencyApproximateModeA(P,N)
        SNR  = signalToNoiseRatioModeA(P);
        TN = (T_N(equation_solver(N),N) - T_N(1/SNR,N));
        seA = exp(1/SNR) * TN / log(2);
    end

 % returns vector of n element with EE for each h0 with given transmission power P 
    function eeA = energyEfficiencyModeA(h, P_D) 
        SNR  = signalToNoiseRatioModeA(P_D);
        eeA = W * log2(1 + h .* SNR ) / (2 * P_D + 2 * P_c_D);
    end

    %returns approximate value of  EE for given transmission power P and
    %number of Taylor series members N used in approximation
    function eeA = energyEfficiencyApproximateModeA(P,N)
        SNR  = signalToNoiseRatioModeA(P);
        TN = (T_N(equation_solver(N),N) - T_N(1/SNR,N));
        P = 2 * (P + P_c_D);
        eeA = exp(1./SNR) .* TN * W ;
        eeA = eeA / (log(2) * P);
    end

    % returns vector of n elements with EE for each h0 with given transmission power P
    function eeB = energyEfficiencyModeB(eta, g10, g20, P_D1)
        P_D2 = P_D1 * ((1-eta) /eta) ^ 4;
        P_D3 = P_D1;
        P_L_B1 =  pathLossD2DLinkB(eta,distance) ; 
        P_L_B2 = pathLossD2DLinkB(1-eta,distance);
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
        P_L_B1 =  pathLossD2DLinkB(eta,distance); 
        P_L_B2 = pathLossD2DLinkB(1-eta,distance);
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
    function seB = spectralEfficiencyApproximateModeB(eta, theta, P, N)
        SNR_B1  = signalToNoiseRatioModeB1(theta,eta,P);
        TN_B1 = T_N(equation_solver(N),N) - T_N(1/SNR_B1,N);
        SNR_B2  = signalToNoiseRatioModeB2(theta,eta,P);
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

%     first derivative of the function ln(EE) 
    function yA = fA(x)
        a= W * N0 * P_L_A * noiseFigure;
        yA = (x.^2 - x .* log(x) - 1 ) / (x.^2 - x .* log(x) - x) +( 2 * a) ./ (2 * P_c_D .* x .^ 2 + 2 * a .* x );
    end
      %first derivative of the EE function 
    function yB = fB(eta, x)
        P_L_B1 = pathLossD2DLinkB(eta,distance); 
        P_L_B2 = pathLossD2DLinkB(1-eta,distance);
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

     function root = bisectionMethodA(fA, a, b, epsilon)
        % Check that that neither end-point is a root
        % and if f(a) and f(b) have the same sign, throw an exception.

        if ( fA(a) == 0 )
            root = a;
            return;
        elseif ( fA(b) == 0 )
            root = b;
            return;
        elseif ( fA(a) * fA(b) > 0 )
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
            if ( fA(c) == 0 )
                root = c;
                return;
            elseif ( fA(c)*fA(a) < 0 )
                b = c;
            else
                a = c;
            end
        end
        root = (b + a)/2;
    end
    function root = bisectionMethodB(fB, a, b, epsilon)
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
