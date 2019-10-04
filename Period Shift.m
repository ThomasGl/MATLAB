%% Periodic Shits to Sin-Shape Signal

n=linspace(-pi,pi,20);
w0=3*pi/4;%same as 6pi/8
x=sin(w0*n);
y=zeros(size(n));
   
%Ploting 3 periods
for p = 1:3
   w1=w0*p; 
   x=sin(w1*n);
   y=x;
   subplot(3,1,p);stem(y);title(sprintf('%d Period Shift', p));
end

