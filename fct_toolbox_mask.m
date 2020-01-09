function [mask_rescaled_v] = fct_toolbox_mask(ntones, duration_sound, fs, radius, coef_smoothing, shape)

% To  mask a 2D image (duration_sound*fs bins on the
% x-axis, and ntones bins on the y axis) with a (smoothed) disc/square/squircle
%%
% E. Ponsot, January 7th 2020
%%

mask_m=zeros(ntones,ntones);
y_center=floor(ntones/2);
x_center=floor(ntones/2);

if strcmp(shape,'disc') 
    % to do a disc
    for ii=1:size(mask_m,2)
        for jj=1:size(mask_m,1)
            if ((ii-x_center)^2+(jj-y_center)^2) < radius^2
                mask_m(jj,ii)=1;
            end
        end
    end

    
elseif strcmp(shape,'square')
     % do a square (2nd version)
    [x,y]=meshgrid(linspace(-1,1,ntones),linspace(-1,1,ntones));
    n=1000; % choose how much rounding: change n (must be even), you will get more or less rounding.
    r=(x.^n+y.^n).^(1/n);
    size_squicle_n=0.5;
    mask_m=(1-abs(double(r>size_squicle_n)));


elseif strcmp(shape,'squircle')
    % do a squircle
    [x,y]=meshgrid(linspace(-1,1,ntones),linspace(-1,1,ntones));
    n=4; % choose how much rounding: change n (must be even), you will get more or less rounding.
    r=(x.^n+y.^n).^(1/n);
    size_squicle_n=0.5;
    mask_m=(1-abs(double(r>size_squicle_n)));

elseif strcmp(shape,'shortsurround')
    for ii=1:size(mask_m,2)
        for jj=1:size(mask_m,1)
            if  abs(jj-y_center)<2*radius
                mask_m(jj,ii)=1;
            end
            
        end
    end   
end

if radius>0
% smoothed the edges of the shape you found
w_gauss_v = gausswin(floor(radius*coef_smoothing)); % choose the max values as the radius of the disc (1 generally)
mask_m=conv2(mask_m,w_gauss_v,'same');
mask_m=conv2(mask_m',w_gauss_v,'same');
end

[mask_rescaled_v] = imresize(mask_m, [ntones, duration_sound * fs]);

% normalize the amplitude of the mask
mask_rescaled_v=mask_rescaled_v./(max(max(abs(mask_rescaled_v))));


