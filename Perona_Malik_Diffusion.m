
%Janavi Kumar
%Perona and Malik Diffusion Algorithm For Image Processing
I = imread('Super_Heros.jpg'); 
I = rgb2gray(I);
I = im2double(I);

figure(1);
imshow(I);

m = size(I,1);
n = size(I,2);

% Add some noise to the image and display it
for i=1:m
    for j=1:n
        Noisy(i,j)=I(i,j)+((rand(1))-.5)/5;
    end
end

figure(2);
imshow(Noisy);

dt = 1; dx = 1; dy = 1;
tfinal = 1.1;
I_recovered = ones(m,n);
A = sparse(m*n);

for i = 1:m
    for j = 1:n
        index = j+(i-1)*n;
        if i==1 || j==1 || i==m || j==n
            A(index,index) = 1;
        else
            %isolate pixels from noisy image
            uxp = (Noisy(i,j+1)-Noisy(i,j))/(dx);
            uxm = (Noisy(i,j)-Noisy(i,j-1))/(dx);
            uyp = (Noisy(i+1,j)-Noisy(i,j))/(dx);
            uym = (Noisy(i,j)-Noisy(i-1,j))/(dx);
            
            %set center, right, left, top and bottom using Perona and Malik
            C = 1 + g(uxm,0) + g(uxp,0) + g(0,uym) + g(0,uyp);
            R = -g(uxp,0)*dt/dx/dx;
            L = -g(uxm,0)*dt/dx/dx;
            T = -g(0,uyp)*dt/dy/dy;
            B = -g(0,uym)*dt/dy/dy;
            
            %set new pixels
            A(index,index) = C;
            A(index,index+1) = R;
            A(index,index-1) = L;
            A(index,index+n) = T;
            A(index,index-n) = B;
            
        end
    end
end

%solve the system for U at time tn+1
for t=0.0:dt:tfinal
    Un = zeros(m*n,1);
    for i=1:m
        for j=1:n
            Un(j+(i-1)*n) = Noisy(i,j);
        end
    end
    Unp1 = A\Un;
end

%create the smoothed out figure and display
for i = 1:m
    for j = 1:n
        I_recovered(i,j) = Unp1(j+((i-1)*n));
    end
end
figure(3);
imshow(I_recovered);

%function to run the Perona-Malik diffusion on given data 
function outval = g(a,b)
    v = 1; %set the coefficient corresponding to space/area that is inspected for diffusion
    K = 3; %set the constant corresponding to edge sensitivity
    s = abs(sqrt(a^2 + b^2));
    outval = v/(1+(s^2/K)); %calculate the diffusion coefficient (D in report)
end


