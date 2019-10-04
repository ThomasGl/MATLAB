%% Rectangular Pulse Function Overlap presentation + Gibbs Phenomena


%Parameters declaration:
T1=1;T=4;N=[1,3,7,19,79];

%Function definition:
time = -T/2:0.01:T/2;    
    
x = zeros(size(time));
x(time<=-T1) = 0;
x( (time > -T1) & (time < T1) ) = 1;
x(time>=T1) = 0;

%The color vector below is critical for an efficient data visualization
Colors = ['y','m','c','r','g','b'];
    
figure;
plot(time,x,'k');title('Pulse Function Reconstruction');hold on;
legend('Original Signal')
%The following loop generates the legend for each function accordingly
for j=1:5
   fReconsta = FindFourierCoeff(time, x, N(j));
   txt=['Reconstructed Signal N=',num2str(N(j))];
   plot(time, fReconsta,Colors(j),'DisplayName',txt);
   hold on;         
end
hold off;
legend show;
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