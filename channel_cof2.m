function [ h ] = channel_cof2( m,n, dmax )
%%  ֻ������һ��������ŵ��ĸ������ɣ�������ŵ�L���ŵ�ʵ��
d = randi([dmax dmax], m, n);
%beta = 10.^(0)*d.^(-2.2);
beta = 1;
g = 1/sqrt(2)*(randn(m ,n)+sqrt(-1)*(randn( m, n)));   %%   L

h=sqrt(beta).*g;                          

end

