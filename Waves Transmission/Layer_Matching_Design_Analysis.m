%% Single Layer Matching
f=linspace(4E14,8E14,1E3);fo=562.5E12;th=linspace(-pi/3,pi/3,1E3);
c=3E8; A=120*pi;
%TE
Ze1=A./cos(th);
Ze2=A./sqrt(1.5-(sin(th)).^2);
Ze3=A./sqrt(2.25-(sin(th)).^2);

dt=c/(4*fo*sqrt(1.5));
b=zeros(1E3,1E3);
for i=1:1E3
    for k=1:1E3
        b((i),(k))=2.*pi.*f(i).*sqrt(minus(1.5,sin(th(k)).^2))./c;
    end
end

Zine=Ze2.*(Ze3+Ze2.*1j.*tan(b*dt))./(Ze2+Ze3.*1j.*tan(b*dt));
GammaIn=abs((Zine-Ze1)./(Zine+Ze1));

figure; subplot(1,2,1);contourf(th, f, GammaIn, 'ShowText', 'on')
colorbar; title('\Gamma_{IN}, TE Polarization')
xlabel('\theta [rad]');ylabel('\itf [Hz]');

%TM
Zm1=A.*cos(th);
Zm2=A.*sqrt(1.5-sin(th).^2)./1.5;
Zm3=A.*sqrt(2.25-sin(th).^2)./2.25;

Zinm=Zm2.*(Zm3+Zm2.*1j.*tan(b*dt))./(Zm2+Zm3.*1j.*tan(b*dt));
GammaIn=abs((Zinm-Zm1)./(Zinm+Zm1));

subplot(1,2,2);contourf(th, f, GammaIn, 'ShowText', 'on')
colorbar; title('\Gamma_{IN}, TM Polarization')
xlabel('\theta [rad]');ylabel('\itf [Hz]');

%% theta=0, frequency sweep
%TE Polarization
Ze10=A;
Ze20=A/sqrt(1.5);
Ze30=A/sqrt(2.25);

b0=2.*pi.*f.*sqrt(1.5)./c;

Zin0=Ze20.*(Ze30+Ze20.*1j.*tan(b0*dt))./(Ze20+Ze30.*1j.*tan(b0*dt));
GammaIn=abs((Zin0-Ze10)./(Zin0+Ze10));

figure; subplot(1,2,1);plot(f, GammaIn); ylabel('\Gamma_{IN}');
xlabel('\itf'); title('\Gamma_{IN}| \theta=0 | TE Polarization');

%TM Polarization

%Z's & beta are the same as TE
%Gamma is the same, we get the same plot

subplot(1,2,2); plot(f, GammaIn); ylabel('\Gamma_{IN}'); xlabel('\itf');
title('\Gamma_{IN}| \theta=0 | TM Polarization')

%% frequency=fo, theta sweep
%TE Polarization
b1=2*pi*fo.*sqrt(minus(1.5,sin(th).^2))./c;
 
Zin=Ze2.*(Ze3+Ze2.*1j.*tan(b1*dt))./(Ze2+Ze3.*1j.*tan(b1*dt));
GammaIn=abs((Zin-Ze1)./(Zin+Ze1));

figure; subplot(1,2,1);plot(th, GammaIn); ylabel('\Gamma_{IN}');
xlabel('\it\theta'); title('\Gamma_{IN}| \itf=f_0 | TE Polarization');

%TM Polarization
%beta is the same as TE
Zinm=Zm2.*(Zm3+Zm2.*1j.*tan(b1*dt))./(Zm2+Zm3.*1j.*tan(b1*dt));
GammaIn=abs((Zinm-Zm1)./(Zinm+Zm1));

subplot(1,2,2);
plot(th, GammaIn); ylabel('\Gamma_{IN}'); xlabel('\it\theta');
title('\Gamma_{IN}| \itf=f_0 | TM Polarization')

%% 2 layers matching
t=[1.257 1.773];
dt1=c/(4*fo*sqrt(t(1))); dt2=c/(4*fo*sqrt(t(2)));

%TE
Ze2=A./sqrt(t(1)-(sin(th)).^2);
Ze3=A./sqrt(t(2)-(sin(th)).^2);
Ze_gl=A./sqrt(2.25-(sin(th)).^2);

b1=zeros(1E3,1E3); b2=zeros(1E3,1E3);
for i=1:1E3
    for k=1:1E3
        b1((i),(k))=2.*pi.*f(i).*sqrt(minus(t(1),sin(th(k)).^2))./c;
        b2((i),(k))=2.*pi.*f(i).*sqrt(minus(t(2),sin(th(k)).^2))./c;
    end
end

geGl=(Ze_gl - Ze3)./(Ze_gl + Ze3);
Gamma = geGl.*exp(-2j*dt2.*b2);
ZL = Ze3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze2)./(ZL + Ze2);
Gamma = Gamma.*exp(-2j*dt1.*b1);
ZL = Ze2.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Ze1)./(ZL + Ze1));

figure; subplot(1,2,1); contourf(th, f, GammaIn, 'ShowText', 'on');
colorbar; title('\Gamma_{IN} | TE | N = 2'); xlabel('\theta [rad]');ylabel('f [Hz]');

%TM

Zm2=A.*sqrt(t(1)-sin(th).^2)./t(1);
Zm3=A.*sqrt(t(2)-sin(th).^2)./t(2);
Zm_gl=A.*sqrt(2.25-sin(th).^2)./2.25;

gmGl=(Zm_gl - Zm3)./(Zm_gl+Zm3);
Gamma = gmGl.*exp(-2j*dt2.*b2);
ZL=Zm3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL-Zm2)./(ZL+Zm2);
Gamma=Gamma.*exp(-2j*dt1.*b1);
ZL = Zm2.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL-Zm1)./(ZL+Zm1));

subplot(1,2,2); contourf(th, f, GammaIn, 'ShowText', 'on');colorbar;
title('\Gamma_{IN} | TM | N = 2'); xlabel('\theta [rad]');ylabel('\itf [Hz]');

%% theta=0, frequency sweep

%TE Polarization
Ze20=A/sqrt(t(1));
Ze30=A/sqrt(t(2));
Ze_gl0=A/sqrt(2.25);

b01=2.*pi.*f.*sqrt(t(1))./c; b02=2.*pi.*f.*sqrt(t(2))./c;

geGl0=(Ze_gl0 - Ze30)/(Ze_gl0 + Ze30);
Gamma = geGl0.*exp(-2j*dt2.*b02);
ZL = Ze30.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze20)./(ZL + Ze20);
Gamma = Gamma.*exp(-2j*dt1.*b01);
ZL = Ze20.*(1+Gamma)./(1-Gamma);

GammaIn=abs((ZL-Ze10)./(ZL+Ze10));

figure; subplot(1,2,1);plot(f, GammaIn); ylabel('\Gamma_{IN}');
xlabel('\itf'); title('\Gamma_{IN}| \theta=0 | TE Polarization | N = 2');

%TM Polarization

%Z's & beta are the same as TE
%Gamma is the same, we get the same plot

subplot(1,2,2); plot(f, GammaIn); ylabel('\Gamma_{IN}'); xlabel('\itf');
title('\Gamma_{IN}| \theta=0 | TM Polarization | N = 2')

%% frequency=fo, theta sweep
%TE Polarization
b01=2*pi*fo.*sqrt(minus(t(1),sin(th).^2))./c;
b02=2*pi*fo.*sqrt(minus(t(2),sin(th).^2))./c;
 
Gamma = geGl.*exp(-2j*dt2.*b02);
ZL = Ze3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze2)./(ZL + Ze2);
Gamma = Gamma.*exp(-2j*dt1.*b01);
ZL = Ze2.*(1+Gamma)./(1-Gamma);

GammaIn=abs((ZL-Ze1)./(ZL+Ze1));

figure; subplot(1,2,1);plot(th, GammaIn); ylabel('\Gamma_{IN}');
xlabel('\it\theta'); title('\Gamma_{IN}| \itf=f_0 | TE Polarization | N = 2')

%TM Polarization

Gamma = gmGl.*exp(-2j*dt2.*b02);
ZL = Zm3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Zm2)./(ZL + Zm2);
Gamma = Gamma.*exp(-2j*dt1.*b01);
ZL = Zm2.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Zm1)./(ZL + Zm1));

subplot(1,2,2); plot(th, GammaIn); ylabel('\Gamma_{IN}');
xlabel('\it\theta'); title('\Gamma_{IN}| \itf=f_0 | TM Polarization | N = 2')

%% 3 layers matching
t=[1.131 1.493 1.970]; dt1=c/(4*fo*sqrt(t(1))); dt2=c/(4*fo*sqrt(t(2))); dt3=c/(4*fo*sqrt(t(3)));
%TE
Ze2=A./sqrt(t(1)-(sin(th)).^2);
Ze3=A./sqrt(t(2)-(sin(th)).^2);
Ze4=A./sqrt(t(3)-(sin(th)).^2);

b3=zeros(1E3,1E3);
for i=1:1E3
    for k=1:1E3
        b1((i),(k))=2.*pi.*f(i).*sqrt(minus(t(1),sin(th(k)).^2))./c;
        b2((i),(k))=2.*pi.*f(i).*sqrt(minus(t(2),sin(th(k)).^2))./c;
        b3((i),(k))=2.*pi.*f(i).*sqrt(minus(t(3),sin(th(k)).^2))./c;
    end
end

geGl=(Ze_gl - Ze4)./(Ze_gl + Ze4);
Gamma = geGl.*exp(-2j*dt3.*b3);
ZL = Ze4.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze3)./(ZL + Ze3);
Gamma = Gamma.*exp(-2j*dt2.*b2);
ZL = Ze3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze2)./(ZL + Ze2);
Gamma = Gamma.*exp(-2j*dt1.*b1);
ZL = Ze2.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Ze1)./(ZL + Ze1));
figure; subplot(1,2,1); contourf(th, f, GammaIn, 'ShowText', 'on');
colorbar; title('\Gamma_{IN} | TE | N = 3'); xlabel('\theta [rad]');ylabel('f [Hz]');

%TM
Zm2=A.*sqrt(t(1)-sin(th).^2)./t(1);
Zm3=A.*sqrt(t(2)-sin(th).^2)./t(2);
Zm4=A.*sqrt(t(3)-(sin(th)).^2)./t(3);

gmGl=(Zm_gl - Zm4)./(Zm_gl+Zm4);
Gamma = gmGl.*exp(-2j*dt3.*b3);
ZL = Zm4.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Zm3)./(ZL + Zm3);
Gamma = Gamma.*exp(-2j*dt2.*b2);
ZL = Zm3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Zm2)./(ZL + Zm2);
Gamma = Gamma.*exp(-2j*dt1.*b1);
ZL = Zm2.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Zm1)./(ZL + Zm1));
subplot(1,2,2); contourf(th, f, GammaIn, 'ShowText', 'on');colorbar;
title('\Gamma_{IN} | TM | N = 3'); xlabel('\theta [rad]');ylabel('\itf [Hz]');

%% theta=0, frequency sweep

%TE Polarization
Ze20=A/sqrt(t(1));
Ze30=A/sqrt(t(2));
Ze40=A/sqrt(t(3));

b01=2.*pi.*f.*sqrt(t(1))./c; b02=2.*pi.*f.*sqrt(t(2))./c; b03=2.*pi.*f.*sqrt(t(3))./c;

geGl0=(Ze_gl0 - Ze40)./(Ze_gl0 + Ze40);
Gamma = geGl0.*exp(-2j*dt3.*b03);
ZL = Ze40.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze30)./(ZL + Ze30);
Gamma = Gamma.*exp(-2j*dt2.*b02);
ZL = Ze30.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze20)./(ZL + Ze20);
Gamma = Gamma.*exp(-2j*dt1.*b01);
ZL = Ze20.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Ze10)./(ZL + Ze10));

figure; subplot(1,2,1);plot(f, GammaIn); ylabel('\Gamma_{IN}');
xlabel('\itf'); title('\Gamma_{IN}| \theta=0 | TE Polarization | N = 3');

%TM Polarization

%Z's & beta are the same as TE
%Gamma is the same, we get the same plot

subplot(1,2,2); plot(f, GammaIn); ylabel('\Gamma_{IN}'); xlabel('\itf');
title('\Gamma_{IN}| \theta=0 | TM Polarization | N = 3')

%% frequency=fo, theta sweep
%TE Polarization
b01=2*pi*fo.*sqrt(minus(t(1),sin(th).^2))./c;
b02=2*pi*fo.*sqrt(minus(t(2),sin(th).^2))./c;
b03=2*pi*fo.*sqrt(minus(t(3),sin(th).^2))./c;

Gamma = geGl.*exp(-2j*dt3.*b03);
ZL = Ze4.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze3)./(ZL + Ze3);
Gamma = Gamma.*exp(-2j*dt2.*b02);
ZL = Ze3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze2)./(ZL + Ze2);
Gamma = Gamma.*exp(-2j*dt1.*b01);
ZL = Ze2.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Ze1)./(ZL + Ze1));
figure; subplot(1,2,1);plot(th, GammaIn); ylabel('\Gamma_{IN}');
xlabel('\it\theta'); title('\Gamma_{IN}| \itf=f_0 | TE Polarization | N = 3')

%TM Polarization

Gamma = gmGl.*exp(-2j*dt3.*b03);
ZL = Zm4.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Zm3)./(ZL + Zm3);
Gamma = Gamma.*exp(-2j*dt2.*b02);
ZL = Zm3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL-Zm2)./(ZL+Zm2);
Gamma = Gamma.*exp(-2j*dt1.*b01);
ZL = Zm2.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Zm1)./(ZL + Zm1));
subplot(1,2,2); plot(th, GammaIn); ylabel('\Gamma_{IN}');
xlabel('\it\theta'); title('\Gamma_{IN}| \itf=f_0 | TM Polarization | N = 3')

%% 4 layers matching
t=[1.0682 1.301 1.71 2.085];

dt1=c/(4*fo*sqrt(t(1))); dt2=c/(4*fo*sqrt(t(2)));
dt3=c/(4*fo*sqrt(t(3))); dt4=c/(4*fo*sqrt(t(4)));

%TE
Ze2=A./sqrt(t(1)-(sin(th)).^2);Ze3=A./sqrt(t(2)-(sin(th)).^2);
Ze4=A./sqrt(t(3)-(sin(th)).^2);Ze5=A./sqrt(t(4)-(sin(th)).^2);

b4=zeros(1E3,1E3);
for i=1:1E3
    for k=1:1E3
        b1((i),(k))=2.*pi.*f(i).*sqrt(minus(t(1),sin(th(k)).^2))./c;
        b2((i),(k))=2.*pi.*f(i).*sqrt(minus(t(2),sin(th(k)).^2))./c;
        b3((i),(k))=2.*pi.*f(i).*sqrt(minus(t(3),sin(th(k)).^2))./c;
        b4((i),(k))=2.*pi.*f(i).*sqrt(minus(t(4),sin(th(k)).^2))./c;
    end
end

geGl=(Ze_gl - Ze5)./(Ze_gl + Ze5);
Gamma = geGl.*exp(-2j*dt4.*b4);
ZL = Ze5.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze4)./(ZL + Ze4);	
Gamma = Gamma.*exp(-2j*dt3.*b3);
ZL = Ze4.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze3)./(ZL + Ze3);
Gamma = Gamma.*exp(-2j*dt2.*b2);
ZL = Ze3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze2)./(ZL + Ze2);
Gamma = Gamma.*exp(-2j*dt1.*b1);
ZL = Ze2.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Ze1)./(ZL + Ze1));
figure; subplot(1,2,1); contourf(th, f, GammaIn, 'ShowText', 'on');
colorbar; title('\Gamma_{IN} | TE | N = 4'); xlabel('\theta [rad]');ylabel('f [Hz]');

%TM

Zm2=(A.*sqrt(t(1)-(sin(th)).^2))./t(1);Zm3=(A.*sqrt(t(2)-(sin(th)).^2))./t(2);
Zm4=(A.*sqrt(t(3)-(sin(th)).^2))./t(3);Zm5=(A.*sqrt(t(4)-(sin(th)).^2))./t(4);

gmGl=(Zm_gl - Zm5)./(Zm_gl+Zm5);
Gamma = gmGl.*exp(-2j*dt4.*b4);
ZL = Zm5.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Zm4)./(ZL + Zm4);	
Gamma = Gamma.*exp(-2j*dt3.*b3);
ZL = Zm4.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Zm3)./(ZL + Zm3);
Gamma = Gamma.*exp(-2j*dt2.*b2);
ZL = Zm3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Zm2)./(ZL + Zm2);
Gamma = Gamma.*exp(-2j*dt1.*b1);
ZL = Zm2.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Zm1)./(ZL + Zm1));
subplot(1,2,2); contourf(th, f, GammaIn, 'ShowText', 'on');colorbar;
title('\Gamma_{IN} | TM | N = 4'); xlabel('\theta [rad]');ylabel('\itf [Hz]');

%% theta=0, frequency sweep
%TE Polarization
Ze20=A/sqrt(t(1));Ze30=A/sqrt(t(2));Ze40=A/sqrt(t(3));Ze50=A/sqrt(t(4));

b01=2.*pi.*f.*sqrt(t(1))./c; b02=2.*pi.*f.*sqrt(t(2))./c;b03=2.*pi.*f.*sqrt(t(3))./c; b04=2.*pi.*f.*sqrt(t(4))./c;

geGl0=(Ze_gl0 - Ze50)./(Ze_gl0 + Ze50);
Gamma = geGl0.*exp(-2j*dt4.*b04);
ZL = Ze50.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze40)./(ZL + Ze40);	
Gamma = Gamma.*exp(-2j*dt3.*b03);
ZL = Ze40.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze30)./(ZL + Ze30);
Gamma = Gamma.*exp(-2j*dt2.*b02);
ZL = Ze30.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze20)./(ZL + Ze20);
Gamma = Gamma.*exp(-2j*dt1.*b01);
ZL = Ze20.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Ze10)./(ZL + Ze10));
figure; subplot(1,2,1);plot(f, GammaIn); ylabel('\Gamma_{IN}');
xlabel('\itf'); title('\Gamma_{IN}| \theta=0 | TE Polarization | N = 4');

%TM Polarization

%Z's & beta are the same as TE
%Gamma is the same, we get the same plot

subplot(1,2,2); plot(f, GammaIn); ylabel('\Gamma_{IN}'); xlabel('\itf');
title('\Gamma_{IN}| \theta=0 | TM Polarization | N = 4')

%% frequency=fo, theta sweep
%TE Polarization
b01=2*pi*fo.*sqrt(minus(t(1),sin(th).^2))./c; b02=2*pi*fo.*sqrt(minus(t(2),sin(th).^2))./c;
b03=2*pi*fo.*sqrt(minus(t(3),sin(th).^2))./c; b04=2*pi*fo.*sqrt(minus(t(4),sin(th).^2))./c;

Gamma = geGl.*exp(-2j*dt4.*b04);  ZL = Ze5.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze4)./(ZL + Ze4);	
Gamma = Gamma.*exp(-2j*dt3.*b03);
ZL = Ze4.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze3)./(ZL + Ze3);
Gamma = Gamma.*exp(-2j*dt2.*b02);
ZL = Ze3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Ze2)./(ZL + Ze2);
Gamma = Gamma.*exp(-2j*dt1.*b01);
ZL = Ze2.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Ze1)./(ZL + Ze1));
figure; subplot(1,2,1);plot(th, GammaIn); ylabel('\Gamma_{IN}');xlabel('\it\theta');title('\Gamma_{IN}| \itf=f_0 | TE Polarization | N = 4')

%TM Polarization
Gamma = gmGl.*exp(-2j*dt4.*b04);  ZL = Zm5.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Zm4)./(ZL + Zm4);	
Gamma = Gamma.*exp(-2j*dt3.*b03);
ZL = Zm4.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Zm3)./(ZL + Zm3);
Gamma = Gamma.*exp(-2j*dt2.*b02);
ZL = Zm3.*(1 + Gamma)./(1 - Gamma);

Gamma = (ZL - Zm2)./(ZL + Zm2);
Gamma = Gamma.*exp(-2j*dt1.*b01);
ZL = Zm2.*(1 + Gamma)./(1 - Gamma);

GammaIn = abs((ZL - Zm1)./(ZL + Zm1));
subplot(1,2,2); plot(th, GammaIn); ylabel('\Gamma_{IN}');xlabel('\it\theta'); title('\Gamma_{IN}| \itf=f_0 | TM Polarization | N = 4')