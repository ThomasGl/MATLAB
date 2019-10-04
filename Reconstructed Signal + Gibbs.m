%% Signal - Reconstruction

%Step Function Definition

T=1;N=3;T1=generate_T1(3,1);t=linspace(-T1,T1,1000);
[x,y,s]=generate_vars(1,t);

subplot(3,1,1);plot(t,x);title('Cos(2\pi * x)');ylim([-1.5 1.5]);

%Dicontinued Cosine + Step Function

subplot(3,1,2);
plot(t,y,'r',t,s,'c');title('Previous + Step Function');ylim([-1.5 1.5]);
legend('y(t)','s(t)');title(legend, 'Wo(s)=1');

%Discontinous Reconstructed Signal

subplot(3,1,3);plot(t,FindFourierCoeff(t,s,100));title('Reconstruced Cosine');

%% Rectangular-Pulse Function Reconstruction

%Gibbs Phenomena Presentation, it happens for functions with discontinuetes
%approximated by Fourier Series

%Function Declaration
time = linspace(-4,4,2000);    
    
Original_Function = zeros(size(time));
Original_Function( time<=-2 ) = 0;
Original_Function( (time > -2) & (time < 2) ) = 1;
Original_Function( time>=2 ) = 0;

%Reconstruction Factors are declared on the vector below:
nVec = [1 2 4 8 16 64 128];
PlotFunctions(time, Original_Function, nVec);
        
function [x,y,s] = generate_vars(T,t)
    x=cos(2*pi*t/T);

    y(t<0)=-1;
    y(t>0)=1;

    s=x.*y;    
end

function [T1] = generate_T1(N,T)
    T1=(N+0.5)*T;
end
%The function below is the applied version of the Harmoic Analysis Theorem
%for function reconstruction.
function [fReconst] = FindFourierCoeff(time, fOrig, k)
    dx = time(2) - time(1);
    a0 = (1/pi) .* sum(fOrig.*dx);
        
    sumVal = 0;
    aN = zeros(k, 1);
    bN = zeros(k, 1);
    for iter = 1:k
        aN(iter) = (1/pi) .* sum(fOrig.*cos(iter*time).*dx);
        bN(iter) = (1/pi) .* sum(fOrig.*sin(iter*time).*dx);
        
        sumVal = sumVal + aN(iter).*cos(iter.*time) + bN(iter).*sin(iter.*time);
    end
    
    fReconst = a0/2 + sumVal;

end

function [] = PlotFunctions(time, fOrig, nVec)
    Colors = ['y','m','c','r','g','b','k'];
    figure; 
    for n=1:7
    
        co=Colors(n);
    
        fReconsta = FindFourierCoeff(time, fOrig, nVec(n));
    
        subplot(2,4,1);plot(time, fOrig, 'k');title('Input Signal');
        subplot(2,4,n+1);plot(time, fReconsta,co);
        str_a = sprintf('Reconstructed Signal N = %d', nVec(n));
        title(str_a);
    end
end