

%% given
clc
clear all
close all
% VSP problem
% depth discretizing 
depth=linspace(0,500,251);
p=length(depth);
% velocity model
v=(3.0+sqrt(depth)/sqrt(1000))*2000;
v(120:160)=v(120:160)-700;
v(50:80)=v(50:80)-600;
% slowness
s=1./v;
% kernel construction
G=triu(ones(p,p));
deltZ=depth(2)-depth(1);
G=G'*deltZ;
% frist arrival time (data)
t=G*s';
% plot
subplot(2,2,1)
plot(v,depth)
set(gca,'ydir','reverse')
xlabel('velocity/m/s')
ylabel('Depth/m')

%% SVD-method
mm=mean(t);
stndrd=mm*0.000002;
t_noisy=t+(stndrd*randn(length(t),1));

 [u1 s1 v1]=svd(G);% pp=rank(G);
 sss=s1(1:245,1:245);
 uuu=u1(:,1:245);
 vvv=v1(:,1:245);
 SS=inv(sss);
 UU=uuu';
 GG=(vvv*SS*UU);
 m_est.svd=GG*t_noisy;
 vv=1./m_est.svd;
 subplot(2,2,2)
 plot(vv,depth)
 set(gca,'ydir','reverse')
 xlabel('velocity')
 ylabel('Depth/m')
 title('SVD')

%% Tikhonov-method
%calculation of optimum alpha__L-curve
alpha=logspace(-3,1,15);

for i=1:length(alpha)
m_tknv(:,i)=(G'*G+alpha(i)*eye(size(G,2)))\(G'*t_noisy);
N1(i)=norm(m_tknv(:,i));
N2(i)=norm(t_noisy-(G*m_tknv(:,i)));

end
 

subplot(2,2,3)
loglog(N1,N2)
xlabel('||t-Gs||')
ylabel('||m||')
title('L-curve')



% opt_alphaa=alpha(6);  %optimum alpha--from L-curve
% S_tknv=((G'.*G)+(opt_alphaa.*eye)).*(G'.*t);
m_est=m_tknv(:,6);
subplot(2,2,4)
plot(1./m_est,depth)
hold on
plot(v,depth,'r')
set(gca,'ydir','reverse')
xlabel('velocity')
ylabel('Depth/m')
title('Tikhonov')
