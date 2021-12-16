%Zhaoyi Jiang
%ME565 HW3
clc
clear all
close all
%Q1
disp('Q1 matlab part 1-----------------------------------------------')
x_bound = 300;
y_bound = 300;
x = linspace(-1.5,1.5,x_bound); y = linspace(-1.5,1.5,y_bound);
max_iteration = 30000;
matrix = zeros(x_bound,y_bound);
u = matrix;

matrix(:,:) = 2;
matrix(:,1:end/2) = 1;
u_1 = matrix;
u_2 = matrix;

for ii=1:1:x_bound
    for jj=1:1:y_bound
        theta = atan2(y(jj),x(ii));
        u(ii,jj) = cos(theta);
    end
end
u_3 = u;

temp_map = zeros(4,y_bound);heat_map = zeros(4,x_bound,y_bound);
temp_map2 = zeros(4,y_bound);heat_map2 = zeros(4,x_bound,y_bound);
temp_map3 = zeros(4,y_bound);heat_map3 = zeros(4,x_bound,y_bound);
count = 1;
stoptime = floor(linspace(1,max_iteration,4));

for itertions=1:max_iteration
    for i=1:1:x_bound
        for j=1:1:y_bound
            if ((x(i)^2+y(j)^2)>1.5)
                matrix(i,j) = u_1(i,j);
                u_2(i,j) = u_1(i,j);
                u(i,j) = u_3(i,j);
            elseif ((x(i)^2+y(j)^2)<=1.5)&&itertions==1
                matrix(i,j) = 0;
                u(i,j) = 0;
                if i>=y_bound/2
                    u_2(i,j) = -1;
                else
                    u_2(i,j) = 1;
                end
            end
        end
    end
    if itertions==stoptime(count)
        temp_map(count,:) = matrix(end/2,:);
        heat_map(count,:,:) = matrix(:,:);
        temp_map2(count,:) = u_2(end/2,:);
        heat_map2(count,:,:) = u_2(:,:);
        temp_map3(count,:) = u(:,end/2);
        heat_map3(count,:,:) = u(:,:);
        count = count + 1;
    end
    cal = del2(matrix);
    matrix(2:x_bound-1,2:y_bound-1) = matrix(2:x_bound-1,2:y_bound-1) + cal(2:x_bound-1,2:y_bound-1);
    cal2 = del2(u_2);
    u_2(2:x_bound-1,2:y_bound-1) = u_2(2:x_bound-1,2:y_bound-1) + cal2(2:x_bound-1,2:y_bound-1);
    cal3 = del2(u);
    u(2:x_bound-1,2:y_bound-1) = u(2:x_bound-1,2:y_bound-1) + cal3(2:x_bound-1,2:y_bound-1);
end

figure(1)
title('Case 1')
subplot(4,1,1); 
imagesc(squeeze(heat_map(1,:,:))); 
colorbar;
title('start'); 
axis square;
caxis([0 2]);
subplot(4,1,2); 
imagesc(squeeze(heat_map(2,:,:))); 
colorbar;
title('beginning'); 
axis square;
caxis([0 2]);
subplot(4,1,3); 
imagesc(squeeze(heat_map(3,:,:)));
colorbar;
title('middle'); 
axis square;
caxis([0 2]);
subplot(4,1,4); 
imagesc(squeeze(heat_map(4,:,:))); 
colorbar;
title('steady state');
axis square;
caxis([0 2]);

figure(2)
hold on;
plot(x,temp_map(1,:),'k-');
plot(x,temp_map(2,:),'r-.');
plot(x,temp_map(3,:),'g-');
plot(x,temp_map(4,:),'b-');
title('temp distribution case1')
legend('start','beginning','middle','ss')
hold off;

figure(3)
title('Case 1, [-1,1] IC')
subplot(4,1,1); 
imagesc(flipud(squeeze(heat_map2(1,:,:)))); 
colorbar;
title('start'); 
axis square;
caxis([-1 2]);
subplot(4,1,2);
imagesc(flipud(squeeze(heat_map2(2,:,:))));
colorbar;
title('beginning');
axis square;
caxis([-1 2]);
subplot(4,1,3); 
imagesc(flipud(squeeze(heat_map2(3,:,:))));
colorbar;
title('middle'); 
axis square;
caxis([-1 2]);
subplot(4,1,4); 
imagesc(flipud(squeeze(heat_map2(4,:,:)))); 
colorbar;
title('steady state'); 
axis square;
caxis([-1 2]);

figure(4)
hold on;
plot(x,temp_map2(1,:),'k-');
plot(x,temp_map2(2,:),'r-.');
plot(x,temp_map2(3,:),'g-');
plot(x,temp_map2(4,:),'b-');
legend('start','beginning','middle','ss')
title('temp distribution case1 [-1,1]')
hold off;

figure(5)
title('Case 2')
subplot(4,1,1); 
imagesc(flipud(squeeze(heat_map3(1,:,:)))); 
colorbar;
title('start'); 
axis square;
caxis([-1 1]);
subplot(4,1,2); 
imagesc(flipud(squeeze(heat_map3(2,:,:)))); 
colorbar;
title('beginning'); 
axis square;
caxis([-1 1]);
subplot(4,1,3); 
imagesc(flipud(squeeze(heat_map3(3,:,:)))); 
colorbar;
title('middle'); 
axis square;
caxis([-1 1]);
subplot(4,1,4); 
imagesc(flipud(squeeze(heat_map3(4,:,:)))); 
colorbar;
title('steady state'); 
axis square;
caxis([-1 1]);

figure(6)
hold on;
plot(x,temp_map3(1,:),'k-');
plot(x,temp_map3(2,:),'r-.');
plot(x,temp_map3(3,:),'g-');
plot(x,temp_map3(4,:),'b-');
legend('start','beginning','middle','ss')
title('temp distribution case2')
hold off;


%%
%Q2
disp('Q2 matlab part 1-----------------------------------------------')
n = [1:100];
a_n = 2./(pi.*n).*sin(n.*pi./2);
b_n = zeros(1,100);
figure
hold on
stem(1:100,a_n)
stem(1:100,b_n,'*')
title('Q2 part1 an and bn');
legend('an','bn');
xlabel('n')
ylabel('an and bn')

disp('Q2 matlab part 2-----------------------------------------------')
n = 10;
x1 = [-2:0.01:2];
f1 = zeros(1,numel(x1));
f1 = 1/2;
for i = 1:n
    a_n = 2*sin(i*pi/2)/(i*pi);
    f1=f1+a_n*cos(i*pi*x1/2);
end

n = 10;
x2 = [-2:0.2:2];
f2 = zeros(1,numel(x2));
f2 = 1/2;
for i = 1:n
    a_n = 2*sin(i*pi/2)/(i*pi);
    f2=f2+a_n*cos(i*pi*x2/2);
end

x3=[-2 -1 -1 1 1 2];
y3=[0 0 1 1 0 0];
figure
plot(x1,f1,x2,f2,x3,y3)
title('Q2 part2')
legend('x=0.01','x=0.2','Square wave')
xlabel('x')
ylabel('y')

% The spacing of x does not change the plor too much.
% Smaller x spacing only makes the plor more smooth,
% but not more accurate.
% The key factor to make the plot accurate is n.



