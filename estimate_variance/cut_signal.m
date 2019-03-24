function x = cut_signal(in,num)
%Description
%this function is cut signal on same pieces
%
%
%
%input
%
%
%output
%
%
    m = length(in)/num;
    x = zeros(m,num);
    bias = 1;
    g  = num;
    for i = 1:1:m
        k = 1;
        for j = bias:1:g

            x(i,k) = in(j);
        k = k + 1;
        end
      
        bias =  bias + num;
        g = g + num;
    end
end