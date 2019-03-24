clear
clc
% In general, you can generate N
% random numbers in the interval (a,b)
% with the formula r = a + (b-a).*rand(N,1).

n = 10;
len = 2^n;
tt = linspace(0,1,len); 
x = zeros(1,len);
SQRT_SNR = 10;
[x,s_n] = wnoise(1,n,SQRT_SNR);%Генерирование сигнала с аддитивным белым Гауссовым шумом

% x = (x * (SQRT_SNR/std(x)));
% wn = randn(1,len);
% s_n = x + wn;

figure
plot(x)
figure
plot(s_n)

num = 64;
x_arr = cut_signal(s_n,num);
dct_arr = zeros(length(x_arr(:,1)),num);

for i = 1:1:length(x_arr(:,1))
    
% dct_arr(i,:) = dct(x_arr(i,:));
dct_arr(i,:) = DCT_ref(x_arr(i,:));

end

dct_arr(:,1) = [];
dct_out_arr = zeros(1,16);

for i = 1:1:length(dct_arr(:,1))
    
dct_out_arr(i) = median(abs( (dct_arr(i,:)) ) );

end

M=median(dct_out_arr);
rezult1=M*1.483;%Результат
Medblock=(1.483).*dct_out_arr;
fprintf('Результат восьмого измерения = %f\n',rezult1);
figure
bar(Medblock)





