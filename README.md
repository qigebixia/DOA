# DOA
coherent signals doa estimation by smoothy space 
clear all;
close all;
clc;

source_number=3;%信元数
sensor_number=10;%阵元数
N_x=1024; %信号长度
snapshot_number=N_x;%快拍数
w=[pi/4 pi/4 pi/3].';%信号频率
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2)+(2*pi*3e8)/w(3))/2;%信号波长  
d=0.5*l;%阵元间距
snr=10;%信噪比
j=sqrt(-1);

source_doa=[5 -5 0];%两个信号的入射角度
A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(1)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(2)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(3)*pi/180)/l)].';%阵列流型
s=10.^(snr/20)*exp(j*w*[0:N_x-1]);%仿真信号

x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%加了高斯白噪声后的阵列接收信号
R=x*x'/snapshot_number;
%对数据协方差矩阵斜对角线上的元素进行平均,即Toeplitz化
dd=zeros(2*sensor_number-1);
for i=-(sensor_number-1):(sensor_number-1)
    c=sum(diag(R,i))/(sensor_number-abs(i));%每一对角线取平均
    dd(i+sensor_number)=c;
end    

for k=1:sensor_number
      R(k,k)=dd(sensor_number);
end

for k=1:(sensor_number-1)
      R(k+1,k)=dd(sensor_number-1);
      R(k,k+1)=dd(sensor_number+1);
end

for k=1:(sensor_number-2)
      R(k+2,k)=dd(sensor_number-2);
      R(k,k+2)=dd(sensor_number+2);
end

for k=1:(sensor_number-3)
      R(k+3,k)=dd(sensor_number-3);
      R(k,k+3)=dd(sensor_number+3);
end

for k=1:(sensor_number-4)
      R(k+4,k)=dd(sensor_number-4);
      R(k,k+4)=dd(sensor_number+4);
end

for k=1:(sensor_number-5)
      R(k+5,k)=dd(sensor_number-5);
      R(k,k+5)=dd(sensor_number+5);
end

for k=1:(sensor_number-6)
      R(k+6,k)=dd(sensor_number-6);
      R(k,k+6)=dd(sensor_number+6);
end

for k=1:(sensor_number-7)
      R(k+7,k)=dd(sensor_number-7);
      R(k,k+7)=dd(sensor_number+7);
end

for k=1:(sensor_number-8)
      R(k+8,k)=dd(sensor_number-8);
      R(k,k+8)=dd(sensor_number+8);
end

for k=1:(sensor_number-9)
      R(k+9,k)=dd(sensor_number-9);
      R(k,k+9)=dd(sensor_number+9);
end
d=s(1,:);
Rx=x*d'/snapshot_number;
w=pinv(R)*Rx;
a=exp(-j*(0:sensor_number-1)'*pi*sin(pi*source_doa(2)/180));
W1p=pinv(R)*a/(a'*pinv(R)*a);
searching_doa=-90:0.1:90;
 for i=1:length(searching_doa)
   a_theta=exp(-j*(0:sensor_number-1)'*pi*sin(pi*searching_doa(i)/180));
   G(i)=abs(W1p'*a_theta);
 end
searching_doa=-90:0.1:90;%线阵的搜索范围为-90~90度
 for i=1:length(searching_doa)
   a_theta=exp(-j*(0:sensor_number-1)'*pi*sin(pi*searching_doa(i)/180));
   P(i)=abs((a_theta)'*w);
 end
plot(searching_doa,10*log(G),'b');
xlabel('入射角/degree');
ylabel('空间谱/dB');
grid on;
