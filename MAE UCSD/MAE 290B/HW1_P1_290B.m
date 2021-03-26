clc; clear all;
kh = linspace(0,2*pi, 1000);
dkh1 = sin(kh);
dkh2 = (2.*sin(kh)) - (sin(kh).*cos(kh));
plot(kh,dkh1);
hold on;
plot(kh, dkh2);
hold on;
plot(kh, kh);
xlabel('kh');
ylabel('k''h');
title('Plot of k''h vs kh')
legend('Central FDA', '2 Point Forward FDA', 'Actual kh');