function [angle_degree,spacing_param] = fct_toolbox_conversion_RateScaleTOAngle(rate,scale,fmin,fmax,duration)
% function that can be used to calculate the parameter scaling_param and angle alpha
% from the rate (Hz) and scale (cycl/oct) 
% useful to merge the two ways of generating ripples signals
% derived using an analytical expression
% E. Ponsot 2019

spacing_param=sqrt(rate^2+(scale/duration*(log(fmax)-log(fmin))).^2 );
angle_degree=acosd(rate/spacing_param);

end

