function [drawings_m] = fct_toolbox_normal_trunc(size_mat ,SD ,trunc_n)

% To  draw the values a matrix of points of size size_mat from a normal distribution with truncated at
% ± trunc_n*SD
%%
% E. Ponsot, January 7th 2020
%%

drawings_m=zeros(size_mat);

if min(size(drawings_m)) > 1  % if this is a matrix
    
    for ii=1:size(drawings_m,1)
        for jj=1:size(drawings_m,2)
            temp_value_n=randn(1);
            while abs(temp_value_n)>trunc_n
                temp_value_n=randn(1);
            end
            drawings_m(ii,jj)=temp_value_n*SD;
        end
    end
    
else % this is a vector
    
    for ii=1:size(drawings_m)
        temp_value_n=randn(1);
        while abs(temp_value_n)>trunc_n
            temp_value_n=randn(1);
        end
        drawings_m(ii)=temp_value_n*SD;
    end
    
end
