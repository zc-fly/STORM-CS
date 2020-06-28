function [Camara_image,x,y]= simulate_camara_image(row,col,unit_pixel);

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
noiseGaussianMean    = 100;    % d.c. background level of camera (counts)
noiseGaussianSigma   =10;     % std dev of noise is 'b' in Thompson-2002
noisePoissonMean     = 0;      % Best set to zero for simplicity

%generate radom fluorecent molecules
molecules_number = 60;
x = rand(molecules_number,1)*(col-8)+4;
y = rand(molecules_number,1)*(row-8)+4;
% x = 4:canshu:col-4;
% y = 4:canshu:row-4;
% [x,y] = meshgrid(x,y);
x_nm = unit_pixel*x;
y_nm = unit_pixel*y;
% molecules_number = size(x);
% molecules_number = molecules_number(1)*molecules_number(2);

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

%% Rsn
%Rsn=10
% Now add camera (Gaussian) and shot (Poisson) noise to simulated data 
%Gaussian noise
G_noise = ones(size(Camara_image))*noiseGaussianMean + floor(noiseGaussianSigma*randn(size(Camara_image)) );
Camara_image1 = Camara_image + G_noise;
%Poisson noise
P_noise = imnoise(uint16( noisePoissonMean*ones(size(Camara_image)) ),'poisson');
Camara_image1 = uint16(Camara_image1) + P_noise;
Camara_image1(Camara_image1<0) = 0;        % Let negative noisy intensities = 0
Camara_image1 = double(Camara_image1);

end
