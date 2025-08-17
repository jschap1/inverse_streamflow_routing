% Calculating Kalman gain for the four-cell basin
%
% 8/31/2023 JRS

% one gage case

H = [0,0,1,0,1,1,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,1,1,0,1;
    0,0,0,0,0,0,0,0,0,0,1,0];

tc = ones(4,4);
TC = [tc, zeros(4,4), zeros(4,4); 
    zeros(4,4), tc, zeros(4,4);
    zeros(4,8), tc];

SC = eye(12);
SC(6,5) = 0.5;
SC(5,6) = 0.5;
% SC = eye(12) + diag(0.5*ones(11,1),-1) + diag(0.5*ones(11,1),1);
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
x = transpose([x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12]);
phi = x*transpose(x);

P = phi.*TC.*SC;
K = P*H'/(H*P*H');
K6 = K(6,:);

%% Plots

% Show when x6+ will be negative, as a function of K, y^-y, and x6/y^

yt = 1; % just choose a value (mm/day)
x6prior = linspace(0,100);
othercellsprior = 3;
x3prior = othercellsprior;
x5prior = othercellsprior;
x8prior = othercellsprior;
yhat = x3prior+x5prior+x6prior+x8prior;

K6_rho0 = (x6prior.^2)./(x3prior.^2+x5prior.^2+x6prior.^2+x8prior.^2);
x6post_rho0 = x6prior + K6_rho0.*(yt-yhat);

K6_rho05 = (x5prior.*x6prior+2*x6prior.^2)./(2*x3prior.^2+2*x5prior.^2+2*x6prior.^2+2*x8prior.^2 + 2*x5prior.*x6prior);
x6post_rho05 = x6prior + K6_rho05.*(yt-yhat);

K6_rho1 = (x6prior.^2 + + x5prior.*x6prior)./(x3prior.^2+x5prior.^2+x6prior.^2+x8prior.^2 + 2*x5prior.*x6prior);
x6post_rho1 = x6prior + K6_rho1.*(yt-yhat);

% Figure

sz = 25;
fs=16;
lw=3;

rat0 = (x6prior.^2)./(x3prior.^2+x5prior.^2+x6prior.^2+x8prior.^2);
figure, plot(rat0)

figure
plot(K6_rho0, x6post_rho0, 'linewidth', lw)
hold on
plot(K6_rho05, x6post_rho05, 'linewidth', lw)
plot(K6_rho1, x6post_rho1, 'linewidth', lw)
xlabel('Proportion of prior uncertainty')
title(['yt = ' num2str(yt)])
ylabel('x6post')
set(gca, 'fontsize', fs)

% [zmin, zind] = min(abs(yhat-yt));
% xline(x6prior(zind)./yhat(zind), 'linewidth', lw)





%%



figure(1)

subplot(1,3,3)
plot(yhat-yt, x6post_rho0, 'linewidth', lw)
hold on
plot(yhat-yt, x6post_rho05, 'linewidth', lw)
plot(yhat-yt, x6post_rho1, 'linewidth', lw)
ylabel('x6post')
xlabel('yhat - yt')
legend('K(\rho)=0','K(\rho)=0.5','K(\rho)=1')
title('x6prior = 10')
set(gca, 'fontsize', fs)

%% Two cell basin - simplest case

H = [1,1];
syms x1 x2
x = [x1;x2];
phi = x*transpose(x);
rho = eye(2); % no corr
rho = [1,1;1,1]; % full corr
P = phi.*rho;
K = P*H'/(H*P*H');

x1 = linspace(0,10);
x2 = linspace(0,10);

y = [0.1, 1, 3];
yhat = zeros(100,100,3);
K0 = zeros(100,100,3);
x1post_0 = zeros(100,100,3);
x1post_05 = zeros(100,100,3);
x1post_1 = zeros(100,100,3);
for k=1:3
    yt=y(k);
    for i=1:100
        for j=1:100
            yhat(i,j,k) = x1(i)+x2(j);
            K0 = x1(i)^2/(x1(i)^2+x2(j)^2);
            K05 = (x1(i)^2 + 0.5*x1(i)*x2(j))/(x1(i)^2+x1(i)*x2(j)+x2(j)^2);
            K1 = (x1(i)^2 + x1(i)*x2(j))/(x1(i)^2+2*x1(i)*x2(j)+x2(j)^2);
            x1post_0(i,j,k) = x1(i) + K0*(yt-yhat(i,j,k));
            x1post_05(i,j,k) = x1(i) + K05*(yt-yhat(i,j,k));
            x1post_1(i,j,k) = x1(i) + K1*(yt-yhat(i,j,k));
            
            K0_2 = x2(j)^2/(x1(i)^2+x2(j)^2);
            K05_2 = (x2(j)^2 + 0.5*x1(i)*x2(j))/(x1(i)^2+x1(i)*x2(j)+x2(j)^2);
            K1_2 = (x2(j)^2 + x1(i)*x2(j))/(x1(i)^2+2*x1(i)*x2(j)+x2(j)^2);
            
            x2post_0 = x2(j) + K0_2*(yt-yhat(i,j,k));
            x2post_05 = x2(j) + K05_2*(yt-yhat(i,j,k));
            x2post_1 = x2(j) + K1_2*(yt-yhat(i,j,k));            
        end
    end
end

%% Figure 

figure(2), clf

subplot(1,2,1)
imagesc(x1,x2,x1post_0(:,:,2)')
colorbar
title(['x_1^+, (y_t = ' num2str(y(2)) ')'])
xlabel('x_1^-')
ylabel('x_2^-')
set(gca, 'ydir', 'normal')
caxis([-1,3])
set(gca, 'fontsize', fs)

subplot(1,2,2)
imagesc(x1,x2,x1post_0(:,:,3)')
colorbar
title(['x_1^+, (y_t = ' num2str(y(3)) ')'])
xlabel('x_1^-')
ylabel('x_2^-')
set(gca, 'ydir', 'normal')
caxis([-1,3])
set(gca, 'fontsize', fs)

colormap(bluewhitered)

%%
figure(16)

subplot(2,3,1)
imagesc(x1,x2,x1post_0(:,:,1)')
colorbar
title(['x1post, y_t = ' num2str(y(1)) ', \rho=0'])
xlabel('x1prior')
ylabel('x2prior')
set(gca, 'ydir', 'normal')
caxis([-1,1])

subplot(2,3,2)
imagesc(x1,x2,x1post_0(:,:,2)')
colorbar
title(['x1post, y_t = ' num2str(y(2)) ', \rho=0'])
xlabel('x1prior')
ylabel('x2prior')
set(gca, 'ydir', 'normal')
caxis([-1,1])

subplot(2,3,3)
imagesc(x1,x2,x1post_0(:,:,3)')
colorbar
title(['x1post, y_t = ' num2str(y(3)) ', \rho=0'])
xlabel('x1prior')
ylabel('x2prior')
set(gca, 'ydir', 'normal')
caxis([-1,1])

subplot(2,3,4)
imagesc(x1,x2,x1post_1(:,:,1)')
colorbar
title(['x1post, y_t = ' num2str(y(1)) ', \rho=1'])
xlabel('x1prior')
ylabel('x2prior')
set(gca, 'ydir', 'normal')
caxis([-1,1])

colormap(cool)

subplot(2,3,5)
imagesc(x1,x2,x1post_1(:,:,2)')
colorbar
title(['x1post, y_t = ' num2str(y(2)) ', \rho=1'])
xlabel('x1prior')
ylabel('x2prior')
set(gca, 'ydir', 'normal')
caxis([-1,1])

subplot(2,3,6)
imagesc(x1,x2,x1post_1(:,:,3)')
colorbar
title(['x1post, y_t = ' num2str(y(3)) ', \rho=1'])
xlabel('x1prior')
ylabel('x2prior')
set(gca, 'ydir', 'normal')
caxis([-1,1])

%% What is bigger? K with corr or K without corr?

x1=1;
x2=2;


K0 = x1^2/(x1^2+x2^2)
K1 = (x1^2 + x1*x2)/(x1^2+x2^2+2*x1*x2)

%%

hold on
plot(yhat-yt, x6post)
plot(yhat-yt, x6post)
legend('x6/yhat=','','')

pointsize = 10;
scatter(x, y, pointsize, z);
