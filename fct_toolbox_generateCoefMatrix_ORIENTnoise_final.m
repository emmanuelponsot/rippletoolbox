function [coef_noise_generated_m, scale_rate_values_m, level_noise_generated_m, phi_noise_generated_m, phi_start_target_v]= fct_toolbox_generateCoefMatrix_ORIENTnoise_final(fmin, fmax, fs, duration_sound, ntones, f0, d, spacing_target,phi_start_list_v, angle_degree_target,nb_angles, SD_noise_dB, nb_noises_n)

%%
% fct_toolbox_generateCoefMatrix_ORIENTnoise_final
%%
% This function creates a signal made of a Gabor/Ripple target embedded in
% an orientation noise made of spectrotemporal rotations of the target
% ripple
%
% ATTENTION !!!: here the ripple is defined as a visual gabor by its angle
% (in degree) and the spacing of the different stripes.
% the equivalence in terms of rate and scale values is given at the output
% of the function
%
% INPUTS:

% fmin, fmax: spanning the fmin to fmax region (in Hz)
% fs : sampling frequency in Hz
% duration_sound : in seconds
% ntones : number of tones (log-spaced across the [fmin-fmax] region) for
% the carrier
% f0 : the first tone freq.
% d : amplitude modulation depth of the ripple (between 0 and 1)
% spacing_target: actually controls the space between two flanks of the ripple 
% phi_start_list_v : vector of the possible starting phases of the ripples
% (chosen randomly among the list)
% angle_degree_target : angle of the target component in degrees
% nb_angles : number of components you want for the noise, there will be ripples equally spaced around 180 degrees 
% SD_noise_dB: SD of the normal distrib in dB
% nb_noises_n: number of noises you want to generate
%
% OUTPUTs:

% coef_noise_generated_m: matrix of amp coefs for the enveloppe of the noise
% scale_rate_values_m: scale and rate values of the superimposed ripples 
% level_noise_generated_m: matrix of the levels samples in the MPS space
% phi_noise_generated_m: matrix of the phases sampled in the MPS space
% phi_start_target_v: vector of the phases at the target location
%
%%
%  E. Ponsot / last seen on 07/01/2020
%%


%% Parameters

% create matrices
level_noise_generated_m=zeros(nb_noises_n, nb_angles);
scale_rate_values_m=zeros(nb_angles, 2);
phi_noise_generated_m=zeros(nb_noises_n, nb_angles);
phi_start_target_v=zeros(1,nb_noises_n);

% create the frequenciy vector and the amplitudes/phases of the different components of the carrier
t_v = linspace(0, duration_sound, duration_sound * fs);
f_v=zeros(1,ntones);
for ii=1:ntones
    f_v(ii) = fmin * (fmax / fmin).^((ii - 1)/(ntones - 1));
end
x_v=log(f_v/f0);



%% Loop

% to make a rotation of a gabor
coef_noise_generated_m=zeros(nb_noises_n, ntones, duration_sound * fs);

% to make a rotation of a gabor
halfSize_x = duration_sound * fs/2;
halfSize_y = ntones/2;
[xx,yy] = meshgrid(linspace(-halfSize_x,halfSize_x,2*halfSize_x), linspace(-halfSize_y,halfSize_y,2*halfSize_y));

for  nn=1:nb_noises_n
    
    coef_noise_m=zeros(ntones, duration_sound * fs);

    % loop on the different orientations to create the orientation noise  
    for jj=1:nb_angles
        
        % calculate the angle
        angle_degree_noise_nn=angle_degree_target+(jj-1)*(180/nb_angles);
        
        % make the appropriate convertion in scale and rate values
        [rate_temp,scale_temp] = fct_toolbox_conversion_AngleTORateScale(angle_degree_noise_nn,spacing_target,fmin,fmax,duration_sound);
        scale_rate_values_m(jj,1)= scale_temp ;
        scale_rate_values_m(jj,2)= rate_temp ;
        
        % draw the phases
        draw_n=randperm(length(phi_start_list_v));
        phi_start_noise_nn=phi_start_list_v(draw_n(1));
        phi_noise_generated_m(nn, jj)=phi_start_noise_nn;
        if jj==1
            phi_start_target_v(nn)=phi_start_noise_nn;
            if angle_degree_target < 0
                disp('ok target')
            else
                disp('ok non-target')
            end
        end
        
        % draw the levels from a normal distrib, truncated
        level_noise_nn=fct_normal_trunc(1,SD_noise_dB,3);
        level_noise_generated_m(nn, jj)=level_noise_nn;
        
        % construct the noise
        coef_noise_m(:,:) = coef_noise_m(:,:)  +  10^((level_noise_nn)/20).*(  1 + d*sin(2*pi*spacing_target*duration_sound/(ntones)*(xx*cosd(angle_degree_noise_nn)*ntones/(duration_sound * fs) + yy*sind(angle_degree_noise_nn))+phi_start_noise_nn));
    end
    coef_noise_generated_m(nn,:,:)=coef_noise_m(:,:);
end


end
