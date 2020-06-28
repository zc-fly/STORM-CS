Canshu = 4;%荧光分子仿真图像的参数
%% Parament setting 
row = 32;
col = 32;
unit_pixel = 160;%nm
eps = 1.5;
base_line = 99;%detemined by the PSF measure result: offset
%% Camara image simulate
[Camara_image,x,y] = simulate_camara_image(row,col,unit_pixel,Canshu);
figure(1);
colormap(gray);
imagesc(Camara_image);
hold on;
plot(x+0.5,y+0.5,'+');%这里x,y为荧光分子的真实位置，但标记在图上时需要加上半个像素对应网格

%% CS-STORM

image_recover = CS_STORM_MY(Camara_image,unit_pixel,base_line,eps);
image_recover(image_recover>600) = 1000;
div =8;
x = div*x;
y = div*y;
figure(2);
colormap(gray);
imagesc(image_recover);
hold on;
plot(x+0.5,y+0.5,'.','Color',[1 0 0]);
%% PSF measurement
% col = 32;
% row = 32;
% X = 1:col;
% Y = 1:row;
% [X ,Y]= meshgrid(X,Y);
% [A,cx,cy,sx,sy,sita,offset]=PSF_measure(X,Y,Camara_image,0.01)%PSF measurement

