clc
clear all
%������ 'blocks'
ind = linspace(0,1,2^10); 
%����������� ������-��� 100
[x,noisyx] = wnoise(1,10,10);%������������� ������� � ���������� ����� ��������� �����
figure
plot(x)
figure
plot(noisyx)
% noisyx = noisyx./(max(noisyx));
%���������� �������� �������
num = 64;
x_arr = cut_signal(noisyx,num);
for i = 1:1:length(x_arr(:,1))
    
   dct_arr(i,:) = DCT_ref(x_arr(i,:));
    
end
% dct_arr(:,1) = [];  
DKPabs = [];
Med = [];
for i = 1:1:length(dct_arr(:,1))
    
   DKPabs(i,:) = abs(dct_arr(i,:));
   tmp = DKPabs(i,:); 
   tmp(1) = [];
   Med(i) = median(tmp);
   tmp = [];
   
end

M=median(Med);
rezult1=M*1.483;%���������
Medblock=1.483*Med;
fprintf('��������� �������� ��������� = %f\n',rezult1);

bar(Medblock)
