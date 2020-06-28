function [Camara_image,x,y]= simulate_camara_image(row,col,unit_pixel,Canshu);

%camara image parament
row_nm = row*unit_pixel;
col_nm = col*unit_pixel;

%PSF parament
A =1;%peak value
N = 1000;%1000 Is the number of possible Photon per pixel
wave_lenth = 670;%nm
NA = 0.8;
Standard_Deviation = 0.21*wave_lenth/NA;

%noise parament
noiseGaussianMean    = 100;    % background level of camera 
noiseGaussianSigma   = 10;   
noisePoissonMean     = 0;      

%generate radom fluorecent molecules
x = 4:Canshu:col-4;
y = 4:Canshu:row-4;
[x,y] = meshgrid(x,y);
x_nm = unit_pixel*x;
y_nm = unit_pixel*y;
molecules_number = size(x);
molecules_number = molecules_number(1)*molecules_number(2);

% generate global grid for camara image
i = unit_pixel/2:unit_pixel:col_nm-unit_pixel/2;
j = unit_pixel/2:unit_pixel:row_nm-unit_pixel/2;
[i,j] = meshgrid(i,j);

%PSF function
%generate camara image
Camara_image = zeros(row,col);
while(molecules_number)
        t = molecules_number;
        I = N*A*exp((-(i-x_nm(t)).^2-(j-y_nm(t)).^2)/(2*(Standard_Deviation^2))); %PSF function
        Camara_image = Camara_image+I;
        molecules_number = molecules_number-1;
        Camara_image(Camara_image < 1) = 0;
end


%add Gaussian and Poisson noise to simulated image 
%Gaussian noise
G_noise = ones(size(Camara_image))*noiseGaussianMean + floor(noiseGaussianSigma*randn(size(Camara_image)) );
Camara_image = Camara_image + G_noise;
%Poisson noise
P_noise = imnoise(uint16( noisePoissonMean*ones(size(Camara_image)) ),'poisson');
Camara_image = uint16(Camara_image) + P_noise;
Camara_image(Camara_image<0) = 0;        %negative noisy  = 0


Camara_image = double(Camara_image);
end
