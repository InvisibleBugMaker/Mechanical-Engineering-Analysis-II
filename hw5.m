clc
clear all
close all

%% Q2 ai
disp('-----------Q2 ai-----------')
figure(1)
hold on
sys = tf(1,[1 3 2])
step(sys)

t = [0:0.2:10];
x = 0.5+0.5.*exp(-2.*t)-exp(-t);
plot(t,x,'r*')
title('ai')
legend('Step Function','Hand Calculation')
hold off

%% Q2 aii
disp('-----------Q2 aii-----------')
figure(2)
hold on
sys = tf(1,[1 3 2])
impulse(sys)

t = [0:0.2:10];
x = exp(-t)-exp(-2.*t);
plot(t,x,'r*')
title('aii')
legend('Step Function','Hand Calculation')
hold off

%% Q2 aiii
disp('-----------Q2 aiii-----------')
figure(3)
hold on
a = [0 1; -2 -3];
b = [0; 1];
c = [1  0];
x0 = [1 ; 0];
sys = ss(a,b,c,[]);
initial(sys,x0)

t = [0:0.2:18];
x = 2*exp(-t)-exp(-2.*t);
plot(t,x,'r*')
title('aiii')
legend('Step Function','Hand Calculation')
hold off

%% Q2 aiv
disp('-----------Q2 aiv-----------')
figure(4)
hold on
a = [0 1; -2 -3];
b = [0; 1];
c = [1  0];
x0 = [0 ; 1];
sys = ss(a,b,c,[]);
initial(sys,x0)

t = [0:0.2:18];
x = exp(-t)-exp(-2.*t);
plot(t,x,'r*')
title('aiv')
legend('Step Function','Hand Calculation')
hold off

%% Q2 bi
disp('-----------Q2 bi-----------')
figure(5)
hold on
sys = tf(1,[1 2 2])
step(sys)

t = [0:0.2:10];
x = 0.5-0.5.*exp(-t).*(sin(t)+cos(t));
plot(t,x,'r*')
title('bi')
legend('Step Function','Hand Calculation')
hold off

%% Q2 bii
disp('-----------Q2 bii-----------')
figure(6)
hold on
sys = tf(1,[1 2 2])
impulse(sys)

t = [0:0.2:10];
x = exp(-t).*sin(t);
plot(t,x,'r*')
title('bii')
legend('Step Function','Hand Calculation')
hold off

%% Q2 biii
disp('-----------Q2 biii-----------')
figure(7)
hold on
a = [0 1; -2 -2];
b = [0; 1];
c = [1  0];
x0 = [1 ; 0];
sys = ss(a,b,c,[]);
initial(sys,x0)

t = [0:0.2:18];
x = exp(-t).*cos(t)+exp(-t).*sin(t);
plot(t,x,'r*')
title('biii')
legend('Step Function','Hand Calculation')
hold off

%% Q2 biv
disp('-----------Q2 biv-----------')
figure(8)
hold on
a = [0 1; -2 -2];
b = [0; 1];
c = [1  0];
x0 = [0 ; 1];
sys = ss(a,b,c,[]);
initial(sys,x0)

t = [0:0.2:18];
x = exp(-t).*sin(t);
plot(t,x,'r*')
title('biv')
legend('Step Function','Hand Calculation')
hold off