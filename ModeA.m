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
    P_L_A =  pathLossD2DLink(distance); 
    
        
%-----------------------------------simulation Mode A----------------------------
    %brake interval  of transmission power range (0,0.25] into m itnervals
    m = 200;
    n = 100000;
    % generate 100000 samples 
%   channel coefficient, exponentially distributed 
%      h0 = exprnd(1/2,1,n);

     h0 = randn(1,n) + 1i * randn(1,n);
     h0 = abs(sqrt(1/2) * h0).^2;
  
% build vectors of EE and SE for simulation and analytical cases
    for i = 1 : m  
        P = P_max_D * i / m;
        average_eeA_sim(i) = mean(energyEfficiencyModeA(h0,P))/10^6;
        average_seA_sim(i) = mean(spectralEfficiencyModeA(h0,P));
        for k = 1:4
            average_eeA_apr_N(k,i) = energyEfficiencyApproximateModeA(P,k*2-1)/10^6;
            average_seA_apr_N(k,i) = spectralEfficiencyApproximateModeA(P, k*2-1);
        end
        
     end
% ------------------------ Energy efficiency subplot-------------------------
    P = linspace( P_max_D/m, P_max_D, m);
    figure;subplot(2,2,1); 
    plot(P, average_eeA_sim,'-');
    hold on;
    for k = 1:4
         plot(P, average_eeA_apr_N(k,:),'LineStyle','--','color',rand(1,3));
    end
    
    %finds optimal transmission power P* for simulation EE
    [M I] = max(average_eeA_sim);
    P_star = P_max_D*I/m
    line([P_star P_star], [M-50 M + 50],'Color', 'blue','LineStyle','-');
    
    % finds optimal transmission power P* for analytical EE
    xA_star = bisectionMethod(@fA, 0.000001, 0.025, 0.00000001);
    P_star = W * N0 * P_L_A * noiseFigure/ xA_star; 
    P_star
    f= mean(energyEfficiencyModeA(h0,P_star))/10^6
    l = mean(spectralEfficiencyModeA(h0,P_star))
    line([P_star P_star], [M-50 M + 50],'Color', 'blue','LineStyle','--');
    hold off;
    legend('Direct-D2D : Simulation','Direct-D2D : Analytical N=1','Direct-D2D : Analytical N=3','Direct-D2D : Analytical N=5','Direct-D2D : Analytical N=7');
    xlabel('UE Transmittion Power (W)');
%     axis([0 0.25 180 390]);
    ylabel('Energy Efficiency (Mbits/Joule)');
    title('Average energy efficiency') 

% ------------------------ Spectral efficiency subplot------------------------%
     subplot(2,2,2);
     plot(P, average_seA_sim,'-');
     hold on;
     for k = 1:4
         plot(P, average_seA_apr_N(k,:),'LineStyle','--','color',rand(1,3));
     end    
     hold off;
     legend('Direct-D2D : Simulation','Direct-D2D : Analytical N=1','Direct-D2D : Analytical N=3','Direct-D2D : Analytical N=5','Direct-D2D : Analytical N=7');
     xlabel('UE Transmittion Power (W)');
     ylabel('Spectral Efficiency (bits/s/Hz)');
     title('Average spectral efficiency');

% ------------------------ Average SNR subplot-------------------------
    subplot(2,2,3);
    plot(P,  10 * log10( mean(h0) * signalToNoiseRatioModeA(P)));
    xlabel('UE Transmittion Power (W)');
    ylabel('Average Received SNR (dB)');
    title('Signal to noise ratio');

% ------------------------ Energy efficiency vs Spectral efficiency subplot-------------------------

    subplot(2,2,4);
    plot(average_seA_sim, average_eeA_sim,'-');
    hold on;
    for k = 1:4
         plot(average_seA_apr_N(k,:), average_eeA_apr_N(k,:),'LineStyle','--','color',rand(1,3));
    end    
    hold off;
    legend('Direct-D2D : Simulation','Direct-D2D : Analytical N=1','Direct-D2D : Analytical N=3','Direct-D2D : Analytical N=5','Direct-D2D : Analytical N=7');
    xlabel('Spectral Efficiency (bits/s/Hz)');
    ylabel('Energy Efficiency (Mbits/Joule)');
    title('Energy efficiency versus spectral efficiency')
      
    
% ---------------------------FUNCTIONS----------------------------------------    
    % SNR  
    function snrA = signalToNoiseRatioModeA(P)
         snrA = P/ ( W * N0 * P_L_A * noiseFigure);
    end

    % return the path loss for D2D  link regarding to the distance d. 
    function pathLoss = pathLossD2DLink(d)           
        pathLoss = (10^14.8) * d^4; 
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

    function root = bisectionMethod(fA, a, b, epsilon)
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
end
