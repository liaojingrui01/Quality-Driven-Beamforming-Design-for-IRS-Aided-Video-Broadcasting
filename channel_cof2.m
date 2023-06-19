function [ h ] = channel_cof2( m,n, dmax )
%%  只需输入一个所需的信道的个数即可，输出是信道L个信道实现
d = randi([dmax dmax], m, n);
%beta = 10.^(0)*d.^(-2.2);
beta = 1;
g = 1/sqrt(2)*(randn(m ,n)+sqrt(-1)*(randn( m, n)));   %%   L

h=sqrt(beta).*g;                          

end

