%% Parament setting 
row = 32;
col = 32;
unit_pixel = 160;%nm
base_line = 99;%detemined by the PSF measure result: offset
%% Camara image simulate
[Camara_image,x,y] = simulate_camara_image(row,col,unit_pixel);
figure(2);
imagesc(Camara_image);
hold on;
plot(x+0.5,y+0.5,'.','Color',[1 0 0]);

%%  CS-STORM
image_recover = CS_STORM_MY(Camara_image,unit_pixel,base_line);
image_recover(image_recover>520) = 1000;

div =8;
x = x*div;
y = y*div;
figure(3);
colormap(gray);
imagesc(image_recover);
hold on;
plot(x+0.5,y+0.5,'.','Color',[1 0 0]);