function image_recover = CS_STORM_MY(Camara_image,unit_pixel,base_line,eps);

%% Parament setting & prime factorization
%CCD Parament
[row,col] = size(Camara_image);
ccd_base_line = base_line;%determine by the PSF measure result
photon_per_count = 1;%in the image simulate,this value has been setted as 1000,here we give 1.
% CVX Arithmetic Parament
eps = 1;%CVX optimization error
scanoverlap = 1;%this value is set dependes on your image(row,col),ensure that you can scan the whole image
boxsize = 7;%pixel
RF = 8;
%L1-Homotopy Arithmetic Parament
maxIteration = 50;
isNonnegative = false;
lambda = 1e-2;
tolerance = 0.2;
stoppingCriterion = -1;
x0 = 0;
%prime factorization
Q = 1;
sub = factor(RF);
N1 = length(sub);
for i = 1 : N1
    subRF(i) = sub(i)*Q;
    Q = subRF(i);
end   
mar = 1;
for i = 1 : N1
    margin(i) = mar*subRF(i)/2;
end
%pixel value theshold using in Cascading control
T = 1;

%% CSR-Prereduction
% Step one /generate global grid for image/*divide unit pixel into 8 grid*
x_MASK = unit_pixel/2:unit_pixel:boxsize*unit_pixel-unit_pixel/2;%Camara image grid
y_MASK = unit_pixel/2:unit_pixel:boxsize*unit_pixel-unit_pixel/2;
[x_MASK,y_MASK] = meshgrid(x_MASK,y_MASK );
x_MASK = x_MASK(:);
y_MASK = y_MASK(:);
x_divMASK = unit_pixel/(2*subRF(1)):unit_pixel/subRF(1):boxsize*unit_pixel-unit_pixel/(2*subRF(1));%image grid after divided
y_divMASK = unit_pixel/(2*subRF(1)):unit_pixel/subRF(1):boxsize*unit_pixel-unit_pixel/(2*subRF(1));
for m = 1:margin(1)
    x_divMASK  = [x_divMASK(1)-unit_pixel/subRF(1) x_divMASK];
    x_divMASK = [x_divMASK x_divMASK(end)+unit_pixel/subRF(1)];
    y_divMASK = [y_divMASK(1)-unit_pixel/subRF(1) y_divMASK];
    y_divMASK = [y_divMASK y_divMASK(end)+unit_pixel/subRF(1)];
end
[x_divMASK,y_divMASK] = meshgrid(x_divMASK,y_divMASK );
x_divMASK = x_divMASK(:);
y_divMASK = y_divMASK(:);

% Step two/ Mesurement matrix generate/*depending on the PSF of the molecules*
len0 = length(x_divMASK);
for i = 1:len0
img_kernel = exp((-((x_MASK-x_divMASK(i))/unit_pixel ).^2-((y_MASK-y_divMASK(i))/unit_pixel).^2)/(2*(1.1046^2)));%PSF Function
MASK(:,i) = img_kernel;
end
c = sum(MASK);
A = sparse(MASK);%（挤出矩阵中的零元素，形成稀疏矩阵）
figure(3)
imagesc(A);
% Step three/ Image reconstruction with CS/*using CVX datebase*
%extract box image
%CVX optimization
Camara_image = (Camara_image - ccd_base_line) * photon_per_count;
Camara_image = Camara_image(:,:);
xbox = 1;
ybox = 1;
%pre-set image_recover
image_recover = zeros(row*subRF(3),col*subRF(3));
while true
        boximage = Camara_image(xbox:xbox+boxsize-1,ybox:ybox+boxsize-1);%extract box image
        b = boximage(:); 

        cvx_quiet true;%CVX optimization
        n = len0;
        L1sum = norm(b,1);
        L2sum = sqrt(L1sum); 
        cvx_begin
           variable x(n)
           minimize(c*x)
           subject to
              x >= 0;
              norm(A*x-b,2) <= eps*L2sum;
         cvx_end
         x = reshape(x(1:len0), [boxsize*subRF(1)+2*margin(1) boxsize*subRF(1)+2*margin(1)]);

        %%  Cascading control
        %referrence:Faster super-resolution imaging of high density molecules via a cascading algorithm based on compressed sensing
        len1 = length(subRF);
        for k = 2 : len1
            [supx,supy] = find( x > T );
            tempRF = subRF(k)/subRF(k-1);
            IX1 = [ ];
            IY1 = [ ];
            x1_divMASK = 1 : (boxsize+mar)*subRF(k);
            y1_divMASK = 1 : (boxsize+mar)*subRF(k);    
            m1 = length(x1_divMASK);
            n1 = length(y1_divMASK);
            validSet = zeros(m1,n1);
            
            L = length(supx); 
            if L ==0;
                x = validSet;
                continue;
            end   
            for i = 1 : L 
                HX_begin = (supx(i) - 1)*tempRF + 1;
                HX_end = supx(i)*tempRF;
                HY_begin = (supy(i) - 1)*tempRF + 1;
                HY_end = supy(i)*tempRF;
                IX = HX_begin : HX_end;
                IY = HY_begin : HY_end;
                validSet(IX,IY) = 1;
                [IX,IY] = meshgrid(IX,IY);
                IX = IX(:);
                IY = IY(:);
                IX1 = [IX1;IX];
                IY1 = [IY1;IY];
            end 
            
            sita0 = find(validSet);
            IX1 = (IX1 - margin(k) - 1)*(unit_pixel/(2^k))+unit_pixel/(2^(k+1));%(unit_pixel) is the edge of the box
            IY1 = (IY1 - margin(k) - 1)*(unit_pixel/(2^k))+unit_pixel/(2^(k+1));
            validSet(x1_divMASK,y1_divMASK) = 0;
            len2 = length(sita0);
            
            m2 = length(x_MASK);%初始化测量矩阵
            n2 = len2;
            MASK_2 = zeros(m2,n2);
            for i = 1 : len2
                img_kernel_2 = exp((-((x_MASK-IX1(i))/unit_pixel ).^2-((y_MASK-IY1(i))/unit_pixel).^2)/(2*(1.1046^2)));
                MASK_2(:,i) = img_kernel_2;
            end
            A2 = sparse(MASK_2);

            [x_out, iterationCount] = SolveHomotopy(A2, b, ...
                                                                        'maxIteration', maxIteration,...
                                                                        'isNonnegative', isNonnegative, ...
                                                                        'stoppingCriterion', stoppingCriterion, ...
                                                                        'groundtruth', x0, ...
                                                                        'lambda', lambda, ...
                                                                        'tolerance', tolerance);
            for i = 1 : len2%put the x into validSet
                t = sita0(i);
                validSet(t) = x_out(i);
            end
            x = validSet;       
        end
         x(x>200) = 200;%if some value is too big,others may can not be seen
         %move the recoverd box image into recoverd image*attention:move the margin away
         image_recover((xbox-1)*subRF(k)+1:((xbox-1)*subRF(k)+1)+boxsize*subRF(k)-1,(ybox-1)*subRF(k)+1:((ybox-1)*subRF(k)+1)+boxsize*subRF(k)-1)...
                               = x(margin(k)+1:end-margin(k),margin(k)+1:end-margin(k)); 
 %调试
% figure(6);
% imagesc((image_recover)), colormap(gray), hold on;
% rectangle('Position',[(ybox-1)*subRF(k)+1,(xbox-1)*subRF(k)+1,boxsize*subRF(k),boxsize*subRF(k)],'linewidth',1,'EdgeColor','y'),pause(0.5);


          %slip to another box image                  
          xbox = xbox + boxsize-scanoverlap;
          if xbox+boxsize -1 > row
              xbox = 1;
              ybox = ybox + boxsize-scanoverlap;
          end
          if ybox+boxsize -1 > col
              break
          end  
          
end   
end