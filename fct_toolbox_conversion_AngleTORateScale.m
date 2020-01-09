function [rate,scale] = fct_toolbox_conversion_AngleTORateScale(angle_degree,spacing_param,fmin,fmax,duration)
% function that can be used to calculate the rate (Hz) and scale (cycl/oct)
% from the parameter scaling_param and alpha
% useful to merge the two ways of generating ripples signals
% (derived using an analytical expression)
% E. Ponsot 2019

scale=spacing_param*duration*sind(angle_degree)/(log(fmax)-log(fmin));
rate=cosd(angle_degree)*spacing_param;

if rate<=0 && scale <=0
    rate=abs(rate);
    scale=abs(scale);
elseif rate>=0 && scale <=0
    rate=-rate;
    scale=abs(scale);
end
    
end