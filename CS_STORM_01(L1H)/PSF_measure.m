function [A,cx,cy,sx,sy,sita,offset]=PSF_measure(X,Y,Z,h)
%用于快速拟合二维高斯分布的实验数据，方法是对x,y,ln(z)做二次多项式拟合。
%适用于那些无偏置，相对噪声小于阈值h的准二维高斯分布数据
%% 平滑图像,尽量消除造成拟合不稳定的凹陷.取那些大于阈值的点带入拟
%这里的滤波在计算单个分子的PSF函数时不能使用，因为会导致标准差拟合误差过大
% f=[1 2 1;2 4 2;1 2 1]/16; 
% Z=filter2(f,Z);    
% for k=1:1    
%     Z=filter2(f,Z);
% end

zmax=max(max(Z));        
xnew=X(find(Z>(zmax*h)));
ynew=Y(find(Z>(zmax*h)));
znew=Z(find(Z>(zmax*h)));
%% 将图像取对数，进行二次多项式拟合
zlognew=log(znew);
a = zeros(max(size(xnew)),6);
n=1;
for i1= 0:2
   for j1=0:2-i1
       a(:,n) = (xnew.^i1).*(ynew.^j1);
       n=n+1;
   end
end
p = (a\zlognew);%zlognew = a*c,c为系数矩阵，而p为c的逆矩阵
%% 由拟合出的多项式系数求出A,cx,cy,sx,sy, sita
c(1)=-p(6);c(2)=-p(3);c(3)=-p(5);c(4)=-p(4);c(5)=-p(2); c(6)=p(1);
sitap=0.5*acot((c(1)-c(2))/c(3));
sxp=sqrt(1/((c(1)-c(2))/cos(2*sitap)+c(1)+c(2)));
syp=sqrt(1/(-(c(1)-c(2))/cos(2*sitap)+c(1)+c(2)));

MA=([-2*(cos(sitap)^2/sxp^2+sin(sitap)^2/syp^2),sin(2*sitap)*(1/syp^2-1/sxp^2);...
    sin(2*sitap)*(1/syp^2-1/sxp^2),-2*(sin(sitap)^2/sxp^2+cos(sitap)^2/syp^2)]);
MB=2*[c(4),c(5)]';
cp=MA\MB;
cxp=cp(1);
cyp=cp(2);
Ap=exp(c(6)+(cos(sitap)*cxp+sin(sitap)*cyp).^2/2/sxp^2+...
  (-sin(sitap)*cxp+cos(sitap)*cyp).^2/2/syp^2);
%% 将线性拟合的结果作为初值带入fminnuc进行最小二乘拟合
tic
offsetp=0;
options = optimset('Display','off','TolFun',1e-4,'LargeScale','off'); %用fminunc拟合，将上面的拟合结果带入作为fminunc的初值
initpar=[Ap,cxp,cyp,sxp,syp,sitap,offsetp];
para=fminunc(@objfun,initpar,options,X,Y,Z);
A=para(1);cx=para(2)-0.5;cy=para(3)-0.5;sx=para(4);sy=para(5);sita=para(6);offset=para(7);%(x,y)-0.5是应为仿真图像点的坐标与真实点坐标差了0.5
time_fit2=toc

B=sum(sum(A));A=A./B;

function S=objfun(p,X,Y,Z)
%% 高斯分布有转动时的误差函数
%A=p(1)
%cx=p(2)
%cy=p(3)
%sx=p(4)
%sy=p(5)
%sita=p(6)
%offset=p(7)
S=sum(sum((p(1)*exp(-0.5/p(4)^2*(cos(p(6))*(X-p(2))+sin(p(6))*(Y-p(3))).^2-0.5/p(5)^2*(-sin(p(6))*(X-p(2))+cos(p(6))*(Y-p(3))).^2)+p(7)-Z).^2));

% zfitp=Ap*exp(-(cos(sitap)*(X-cxp)+sin(sitap)*(Y-cyp)).^2/2/sxp^2-(-sin(sitap)*(X-cxp)+cos(sitap)*(Y-cyp)).^2/2/syp^2)+offsetp;

%写在主函数里的部分
% Z = double(Camara_image);
%  [x y] = find( Z==max(max(Z)) );
% Z = Z( (x-1):(x+1),(y-1):(y+1)) ;
%  [X,Y] = meshgrid(x-1:1:x+1,y-1:1:y+1);
%  B=sum(sum(Z));Z=Z./B;
%  k = 1/Z(2,2);
%  Z = k*Z;
