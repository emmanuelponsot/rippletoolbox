function [mask_rescaled_v] = fct_toolbox_smooth(ntones, duration_sound, fs, size_edge, radius, coef_smoothing)

% To create a 2D mask (duration_sound*fs bins on the
% x-axis, and ntones bins on the y axis),
% parameters of this mask are: size_edge, radius, coef_smoothing
%%
% E. Ponsot, January 7th 2020
%%
mask_m=ones(ntones,ntones);

%% remove edges
for ii=1:size(mask_m,2)
    for jj=1:floor(size_edge)
        mask_m(ii,jj)=0;
    end
    for jj=size(mask_m,1)-floor(size_edge)+1:size(mask_m,1)
        mask_m(ii,jj)=0;
    end
end

for jj=1:size(mask_m,1)
    for ii=1:floor(size_edge)
        mask_m(ii,jj)=0;
    end
    for ii=size(mask_m,1)-floor(size_edge)+1:size(mask_m,1)
        mask_m(ii,jj)=0;
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

