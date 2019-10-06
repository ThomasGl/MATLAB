%% Full Analysis on Decimation, Interpolation, Upsampling & Downsampling Process to DSP

%% Original Function Definition
sample_space=linspace(-25,25,51);x_value = zeros(size(sample_space));

for i= 1:length(sample_space)
x_value(i)= ((2/sample_space(i))*sin(sample_space(i)*pi/8))^2;
end

x_value(sample_space==0)=(2*pi/8)^2;

figure;
%% Fourier Transform of the Original Signal + Plottings
FTx=abs(fft(x_value));
omega=linspace(-2*pi,2*pi,length(x_value));
subplot(5,4,1);plot(omega,FTx);title('Original FT');
subplot(5,4,2);stem(sample_space,x_value);title('Original');

%% Downsample - Factor 3 - + Fourier Transform
sampled3=downsample(sample_space,3);
xd3=downsample(x_value,3);
FTxd3=abs(fft(xd3));
omegad3=linspace(-2*pi,2*pi,length(xd3));
subplot(5,4,3);plot(omegad3,FTxd3);title('Downsample 3 FT');
subplot(5,4,4);stem(sampled3,xd3,':r*');title('Downsample 3');

%% Decimation  - Factor 3 -  + Fourier Transform
samplede3=decimate(sample_space,3);
xde3=decimate(x_value,3);
FTxde3=abs(fft(xde3));
omegade3=linspace(-2*pi,2*pi,length(xde3));
subplot(5,4,5);plot(omegade3,FTxde3);title('Decimate 3 FT');
subplot(5,4,6);stem(samplede3,xde3,':g*');title('Decimate 3');

%% Downsample - Factor 5 - + Fourier Transform
sampled5=downsample(sample_space,5);
xd5=downsample(x_value,5);
FTxd5=abs(fft(xd5));
omegad5=linspace(-2*pi,2*pi,length(xd5));
subplot(5,4,7);plot(omegad5,FTxd5);title('Downsample 5 FT');
subplot(5,4,8);stem(sampled5,xd5,':m*');title('Downsample 5');

%% Decimation  - Factor 5 -  + Fourier Transform
samplede5=decimate(sample_space,5);
xde5=decimate(x_value,5);
FTxde5=abs(fft(xde5));
omegade5=linspace(-2*pi,2*pi,length(xde5));
subplot(5,4,9);plot(omegade5,FTxde5);title('Decimate 5 FT');
subplot(5,4,10);stem(samplede5,xde5,':k*');title('Decimate 5');

%% Upsample - Factor 3 - + Fourier Transform
sampleu3=upsample(sample_space,3);
xu3=upsample(x_value,3);
FTxu3=abs(fft(xu3));
omegau3=linspace(-2*pi,2*pi,length(xu3));
subplot(5,4,11);plot(omegau3,FTxu3);title('Up sample 3 FT');
subplot(5,4,12);stem(sampleu3,xu3,':k*');title('Up sample 3');

%% Interpolation  - Factor 3 -  + Fourier Transform
samplein3=interp(sample_space,3);
xin3=interp(x_value,3);
FTxin3=abs(fft(xin3));
omegain3=linspace(-2*pi,2*pi,length(xin3));
subplot(5,4,13);plot(omegain3,FTxin3);title('Interpol 3 FT');
subplot(5,4,14);stem(samplein3,xin3,':y*');title('Interpol 3');

%% Upsample - Factor 5 - + Fourier Transform
sampleu5=upsample(sample_space,5);
xu5=upsample(x_value,5);
FTxu5=abs(fft(xu5));
omegau5=linspace(-2*pi,2*pi,length(xu5));
subplot(5,4,15);plot(omegau5,FTxu5);title('Up sample 5 FT');
subplot(5,4,16);stem(sampleu5,xu5,':r*');title('Up sample 5');

%% Interpolation - Factor 5 - + Fourier Transform
samplein5=interp(sample_space,5);
xin5=interp(x_value,5);
FTxin5=abs(fft(xin5));
omegain5=linspace(-2*pi,2*pi,length(xin5));
subplot(5,4,17);plot(omegain5,FTxin5);title('Interpol 5 FT');
subplot(5,4,18);stem(samplein5,xin5,':g*');title('Interpol 5');

%% Resampling - Factor 2/3 - + Fourier Transform
sample23=resample(sample_space,2,3);
x23=resample(x_value,2,3);
FTx23=abs(fft(x23));
omegax23=linspace(-2*pi,2*pi,length(x23));
subplot(5,4,19);plot(omegax23,FTx23);title('Rational 2/3 FT');
subplot(5,4,20);stem(sample23,x23,':b*');title('Rational 2/3');

%% Sine Wave Manipulation

%sExtracting audio_file
mu=audioread('sin.wav');

%Down->Inter
mud4=downsample(mu,4);
muin4=interp(mud4,4);

%Inter->Down
mu4in=interp(mu,4);
mu4d=downsample(mu4in,4);

%Plotting + Fourier Transform
muf=abs(fft(mu));
muin4f=abs(fft(muin4));
mu4df=abs(fft(mu4d));
omega_orginal = linspace(-2*pi,2*pi,length(mu));
omega_2 = linspace(-2*pi,2*pi,length(muin4f));
omega_3 = linspace(-2*pi,2*pi,length(mu4df));
figure;
subplot(2,2,[1 2]);plot(omega_orginal,muf);title('Original');
subplot(2,2,3);plot(omega_2,muin4f);title('Down->Inter');
subplot(2,2,4);plot(omega_3,mu4df);title('Inter->Down');
