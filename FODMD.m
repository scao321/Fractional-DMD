% [1] Igor Podlubny (2025). Mittag-Leffler function (https://www.mathworks.com/matlabcentral/fileexchange/8738-mittag-leffler-function), MATLAB Central File Exchange. Retrieved January 18, 2025.
% [2] R. Garrappa and M. Popolizio, Computing the matrix Mittag–Leffler function with applications to fractional calculus, Journal of Scientific Computing, 2018, 17(1), 129-153
% [3] Le Clainche, Soledad; Vega, José M.  (2020), "Book Higher Order Dynamic Mode Decomposition and its applications: MATLAB codes”, Mendeley Data, V1, doi: 10.17632/z8ks4f5vy5.1
clear
clc
close all
addpath("ml_matrix/") % [2] 
tf=10;
dt=0.001;
t=0:dt:tf;
xi=linspace(-10,10,100);
[Xgrid,T]=meshgrid(xi,t);
a=[2,1];
na=[0.5,0];
b=0;
nb=0;
u=0*ones(length(t),1);
y=cos(Xgrid').*mlf(0.5,1,(1+1i)*t.^(0.5))'+cos(10+2*Xgrid').*mlf(0.5,1,(-2+1i)*t.^(0.5))';%[1]
lod= 0.2; %length of the input data
y1=y(:,1:floor(lod*length(t)));
%%
lom=1;%length of the memory
memory=floor(lom*length(y1));
alpha=.5;
w=foweight(alpha-1,memory,length(xi));
dy1=creatfeature(y1(:,2:end)-y1(:,1:end-1),memory,1);
nablaY=w*dy1;
%% Dimension reduction [3]
x1=y1(:,1:end-1);
x2=nablaY;
[U,S,V]=svd(x1,'econ');
sigma=diag(S);
n=length(sigma);
NormS=norm(sigma,2);
r=0;
for k=1:n
    R(k)=norm(sigma(k:n),2)/NormS;
    if R(k)>1e-10
        r=r+1;
    end
end
U_r = U(:, 1:r); % truncate to rank-r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
Atilde = U_r' * x2 * V_r / S_r; % low-rank dynamics

[W_r, D] = eig(Atilde);
Phi=U_r*W_r;
lambda = diag(D); % discrete-time eigenvalues
%% FODMD mode amplitudes b
x11 = x1(:, 1);
b = pinv(Phi)*x11;

%% Reconstruction
 i=1;
 for t1=0:dt:tf % more smaller dt is, more accurate results are
 y2(:,i)=Phi*ml_matrix(D*t1^alpha/(dt^alpha),alpha,1)*b; % [2]
 i=i+1;
 end
%%
figure(1)
plot(linspace(-10,10,100),Phi(:,1)*b(1),LineWidth=2) 
hold on
plot(linspace(-10,10,100),Phi(:,2)*b(2),LineWidth=2)
plot(linspace(-10,10,100),cos(10+2*linspace(-10,10,100)),LineWidth=2,LineStyle="--") 
plot(linspace(-10,10,100),cos(linspace(-10,10,100)),LineWidth=2,LineStyle="--")
legend('FO Mode 1','FO Mode 2', 'Mode 1','Mode 2')
figure(2)
plot(t,real(y(1,:)),LineWidth=3)
hold on
plot(t,real(y2(1,:)),'--o',MarkerIndices=1:700:length(y2(1,:)),LineWidth=3,MarkerSize=8)
line([2 2], [-3,3], 'LineWidth', 2,'Color','black','LineStyle','-');
ylim([-3,3])
grid on
xlabel('time','FontSize',15)
ylabel('$f(t)$', 'Interpreter', 'latex','FontSize',15)
legend('Real','FO-DMD','',FontSize=13)
title('time evolution',FontSize=15)
pbaspect([3 1 1])
figure(3)
layout = tiledlayout(1,2);
nexttile
surfl(Xgrid',T',real(y))
shading interp; colormap(gray); 
xlabel('X')
ylabel('time')
nexttile
surfl(Xgrid',T',real(y2))
shading interp; colormap(gray); 
xlabel('X')
ylabel('time')

%%
norm(y2(1,:)-y(1,:))/norm(y(1,:),2)
%%
function xnew=creatfeature(x,L,numb)
[m,n]=size(x);
xnew=zeros(m,n);
for i=1:m
r=zeros(1,L);
r(1:numb)=x(i,1:numb);
xnew(L*(i-1)+1:L*i,:)=toeplitz(r,x(i,:));
end
end

function w=foweight(alpha,L,r)
w1=1;
w=zeros(r,L);
for i=2:L
w1(i)=w1(i-1)*(1-(alpha+1)/(i-1));
end
cells = repmat({w1}, 1, r);
w=blkdiag(cells{:});
end