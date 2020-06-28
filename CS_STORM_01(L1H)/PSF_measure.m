function [A,cx,cy,sx,sy,sita,offset]=PSF_measure(X,Y,Z,h)
%���ڿ�����϶�ά��˹�ֲ���ʵ�����ݣ������Ƕ�x,y,ln(z)�����ζ���ʽ��ϡ�
%��������Щ��ƫ�ã��������С����ֵh��׼��ά��˹�ֲ�����
%% ƽ��ͼ��,�������������ϲ��ȶ��İ���.ȡ��Щ������ֵ�ĵ������
%������˲��ڼ��㵥�����ӵ�PSF����ʱ����ʹ�ã���Ϊ�ᵼ�±�׼�����������
% f=[1 2 1;2 4 2;1 2 1]/16; 
% Z=filter2(f,Z);    
% for k=1:1    
%     Z=filter2(f,Z);
% end

zmax=max(max(Z));        
xnew=X(find(Z>(zmax*h)));
ynew=Y(find(Z>(zmax*h)));
znew=Z(find(Z>(zmax*h)));
%% ��ͼ��ȡ���������ж��ζ���ʽ���
zlognew=log(znew);
a = zeros(max(size(xnew)),6);
n=1;
for i1= 0:2
   for j1=0:2-i1
       a(:,n) = (xnew.^i1).*(ynew.^j1);
       n=n+1;
   end
end
p = (a\zlognew);%zlognew = a*c,cΪϵ�����󣬶�pΪc�������
%% ����ϳ��Ķ���ʽϵ�����A,cx,cy,sx,sy, sita
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
%% ��������ϵĽ����Ϊ��ֵ����fminnuc������С�������
tic
offsetp=0;
options = optimset('Display','off','TolFun',1e-4,'LargeScale','off'); %��fminunc��ϣ����������Ͻ��������Ϊfminunc�ĳ�ֵ
initpar=[Ap,cxp,cyp,sxp,syp,sitap,offsetp];
para=fminunc(@objfun,initpar,options,X,Y,Z);
A=para(1);cx=para(2)-0.5;cy=para(3)-0.5;sx=para(4);sy=para(5);sita=para(6);offset=para(7);%(x,y)-0.5��ӦΪ����ͼ������������ʵ���������0.5
time_fit2=toc

B=sum(sum(A));A=A./B;

function S=objfun(p,X,Y,Z)
%% ��˹�ֲ���ת��ʱ������
%A=p(1)
%cx=p(2)
%cy=p(3)
%sx=p(4)
%sy=p(5)
%sita=p(6)
%offset=p(7)
S=sum(sum((p(1)*exp(-0.5/p(4)^2*(cos(p(6))*(X-p(2))+sin(p(6))*(Y-p(3))).^2-0.5/p(5)^2*(-sin(p(6))*(X-p(2))+cos(p(6))*(Y-p(3))).^2)+p(7)-Z).^2));

% zfitp=Ap*exp(-(cos(sitap)*(X-cxp)+sin(sitap)*(Y-cyp)).^2/2/sxp^2-(-sin(sitap)*(X-cxp)+cos(sitap)*(Y-cyp)).^2/2/syp^2)+offsetp;

%д����������Ĳ���
% Z = double(Camara_image);
%  [x y] = find( Z==max(max(Z)) );
% Z = Z( (x-1):(x+1),(y-1):(y+1)) ;
%  [X,Y] = meshgrid(x-1:1:x+1,y-1:1:y+1);
%  B=sum(sum(Z));Z=Z./B;
%  k = 1/Z(2,2);
%  Z = k*Z;
