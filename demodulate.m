function []=demodulate()
close all;
clc;
%% 参数
Code = 1;
code = Code;
SF = 3;
number_bins = 2^SF; 
%%信号的参数设置
T=10e-6; %pulse duration10us
B=50e6; %chirp frequency modulation bandwidth 30MHz
K=B/T;%chirp slope
Fs=2*B;Ts=1/Fs;%sampling frequency and sample spacing
N=int32(T/Ts);
f0 = 0*10e6; %发射信号时的瞬时频率，也就是信号有效区间发射信号的中心频率
t=linspace(0,T,N);
n=0:N-1;

A0 = 10e-3;%发射信号的振幅
Phi0 = pi/2; %发射信号的随机初相
%% 产生线性调频信号
f_0=code*B/number_bins;
k=0;
St1 = zeros(1,N);
for i=1:1:N
    if f0+f_0+K*t(i) < f0+B 
        k=k+1;
        St1(i)=exp(1j*(2*pi*(f0+f_0)*t(i)+pi*K*t(i)^2)); %线性调频信号复数表达式
    else
        St1(i)=exp(1j*(2*pi*(f0)*t(i-k)+pi*K*t(i-k)^2));
    end
end
           
St1_1=exp(1j*(2*pi*(f0+B)*t-pi*K*t.^2)); %downchirp 

%% 添加噪声
Noise = 0*(wgn(1,N,1)+1j*wgn(1,N,1));
St1_0 = St1;
St1 = (St1+Noise);
St1_2 = [St1,Noise];
St1_2 = [Noise,St1_2]; %创建一段含有upchirp的样本

% figure
% t1=linspace(0,T,3*N);
% plot(t1*1e6,real(St1_2));


%% 引入真实数据

    filename='data/test.dat';
    fi_1 = fopen(filename,'rb');
    xx_1=read_complex_binary(filename,1e6);
    cpdata = xx_1(3.6e5:3.8e5-1);
    cpdata_0 = cpdata(1:5000);
    cpdata_0 = cpdata_0(2000-100:4000-101);
    plot(imag(cpdata_0))
    
    
    
    filename1 = 'data/downchirp_imag.dat';
    filename2 = 'data/downchirp_real.dat';
    I = load(filename1,'%f');
    Q = load(filename2,'%f');
    cpdata_1 = complex(Q,I);
    cpdata_1 = resample(cpdata_1,2000,200);
    figure;

    plot(real(cpdata_1))
    
    St1 = cpdata_0';
    St1_1 = cpdata_1';
%    St1 = St1.*St1_1;
    N = 2000;
    St1_2 = cpdata(700:3000)';
    
    
    
%% 解调
% [M,index] = get_shift_fft(St1,St1_1,N,number_bins);
% index1 = max_frequncy_gradient_idx(St1,N,number_bins);
% disp(['真实为:']);
% disp(code)
% disp(['第一种解调为:']);
% disp(index);
% disp(['第二种解调为:']);
% disp(index1);
 %St1 = St1.*St1_1;
 %upchirp_ifeq = instantaneous_frequency(St1_0,N);
 %[correlation,index_chirp] = detect_upchirp(St1_2,upchirp_ifeq,N);
 %[correlation1,index_chirp1] = detect_upchirp1(St1_2,upchirp_ifeq,N);
 %[correlation2,index_chirp2] = detect_upchirp2(St1_2,St1_1,N)
 


%% 画图
t=linspace(0,T,N);

figure
subplot(411)
plot(t*1e6,real(St1));
xlabel('Time（us）');
ylabel('Amplitude(Watts)')
title('upchirp*downchirp');
grid on;
axis tight;
 
subplot(412)
plot(t*1e6,real(St1_1));
xlabel('Time（us）');
ylabel('Amplitude(Watts)')
title('downchirp信号的实部');
grid on;
axis tight;

St1 = abs(St1.*St1_1);
subplot(413)
freq=linspace(-Fs/2,Fs/2,N);
plot(freq*1e-6,abs(fft(St1))); %先对St做傅里叶变换得到频谱，在取幅度值，然后将其移动到频谱中心
xlabel('Frequency（MHz）');
ylabel('Amplitude(Watts)')
title('线性调频信号的频谱');
grid on;
axis tight;

x=St1;
y=fft(x,N);    %对信号进行快速Fourier变换
mag=abs(y);     %求得Fourier变换后的振幅
f=n*Fs/N;    %频率序列
mag_1=mag(1:N/2)+mag(N/2+1:N);
subplot(414)
plot(mag_1); %绘出Nyquist频率之前随频率变化的振幅
xlabel('Frequency（MHz）');
ylabel('Amplitude(Watts)');
title('Samples');
grid on;

end

%% FRFT
function power = my_frft(data)
    f=100;
    T=1e-3;

    xx = data;
    point = T*f;
    pa=0:0.1:2;
    u=1:1:length(xx);
    Pp3=zeros(length(u),length(pa));
    for k=1:1:length(pa)
        p3=frft(xx,pa(k));
        %vx3=frft(vx,a);
        %vy3=frft(vy,a);
        Pp3(:,k)=p3.*conj(p3)/point;
    end
    power=max(max(max(Pp3)));
    [a b c]=find(Pp3==power);

    i=1;
    while pa(b)==1              %当p=1 阶为1的时候降低步进值
        pa=1-1/(10^i):0.1/10^i:1+1/(10^i);
        u=1:1:length(xx);
        Pp3=zeros(length(u),length(pa));
        for k=1:1:length(pa)
            p3=frft(xx,pa(k));
            %vx3=frft(vx,a);
            %vy3=frft(vy,a);
            Pp3(:,k)=p3.*conj(p3)/point;
        end
        power=max(max(max(Pp3)));
        [a b c]=find(Pp3==power);
        if i==3
            break;
        end
        i=i+1;
    end
     figure
     surf(pa,u,Pp3);
     
    
    

end
%% 得到环境中的FRFT峰值
function power = huanjing(filename)
    xx_1=read_data(filename);
    figure(7)
    plot(abs(xx_1))
    start=3*1e5;
    window_len=1000;
    data=xx_1(start:start+window_len-1); 
    plot(abs(data))
    power=my_frft(data);
end
%% 波形平滑一些
function data = my_smooth(cpdata,thr)
    m=median(abs(cpdata));
    for i=1:1:length(cpdata)
        cha=abs(cpdata(i))/m;
        if cha>thr
            cpdata(i)=cpdata(i)*(1-0.4*thr);
        end
    end
    data=cpdata;
end
%% 得到瞬时频率
function ifreq=instantaneous_frequency(samples,window)
    if window < 2
        disp('WARNING : samples size < 2 !');
    end
    ifreq = zeros(1,window);
    for i=1:1:window-1
        iphase_1 = angle(samples(i));
        iphase_2 = angle(samples(i+1));
        while ( (iphase_2 - iphase_1) > pi) 
            iphase_2 = iphase_2 - 2*pi;
        end
        while ( (iphase_2 - iphase_1) < -pi)
            iphase_2 = iphase_2 + 2*pi;
        end
        ifreq(i) = iphase_2 - iphase_1 ;
    end
    ifreq(window) = ifreq(window-1);
end
%% 得到瞬时相位
function iphase=instantaneous_phase(samples,window)
    iphase(1) =  angle(samples(0));
    for i = 2:1:window
        iphase(i) = angle(samples(i));
        while ( (iphase(i) - iphase(i-1)) > pi) 
            iphase(i) = iphase(i) - 2*pi;
        end
        while ( (iphase(i) - iphase(i-1)) < -pi)
            iphase(i) = iphase(i) + 2*pi;
        end
    end
end
%% 第一种解调方式，直接求频率最大处
function bin_idx = max_frequncy_gradient_idx(samples,samples_per_symbol,number_of_bins)
    samples_ifreq = instantaneous_frequency(samples,samples_per_symbol);
    decim_factor = samples_per_symbol / number_of_bins;
    samples_ifreq_avg = zeros(1,number_of_bins);
 
    for i = 1:1:number_of_bins
        samples_ifreq_avg(i) = sum(samples_ifreq((i-1)*decim_factor+1:i*decim_factor));
        %samples_ifreq_avg(i) = samples_ifreq_avg(i) / decim_factor;
    end
    max_gradient = 0.1;
    gradient = 0;
    max_index = 0;
    for i=1:1:number_of_bins-1
        gradient = samples_ifreq_avg(i) -samples_ifreq_avg(i+1);
        if gradient > max_gradient
            max_gradient = gradient;
            max_index = i;
        end
    end
    
    
    bin_idx = mod((number_of_bins - max_index), number_of_bins);
   
    %bin_idx = (number_of_bins - max_index);
end
%% 第二种解调方式，upchirp*downchirp后做FFT，求bin域
function [M,bin_idx] = get_shift_fft(samples,downchirp,samples_per_symbol,number_of_bins)
    mult_hf = zeros(1,samples_per_symbol);
    for i=1:1:samples_per_symbol
        mult_hf(i) = samples(i) * downchirp(i);
    end
    get_fft = fft(mult_hf,samples_per_symbol);
    mag = abs(get_fft);
    mag1 = mag(1:samples_per_symbol/2) + mag(samples_per_symbol/2+1:samples_per_symbol);
    size_bins = floor(samples_per_symbol / (2 * number_of_bins));
    k = 0;
    Index_peak = zeros(1,floor(samples_per_symbol/size_bins));
    for i=1:size_bins:samples_per_symbol/2-size_bins
        peak = sum(mag1(i:i+size_bins));
        k = k+1;
        Index_peak(k) = peak;
    end
    [M,I] = max(Index_peak);
    M = max(mag1);
    bin_idx = I-1;
end
%% upchirp 乘上downchirp 做fft
function [M,outlier] = mult_fft(samples,downchirp,window)
    mult_hf = zeros(1,window);
    for i=1:1:window
        mult_hf(i) = samples(i) * downchirp(i);
    end
    get_fft = fft(mult_hf,window);
    mag = abs(get_fft);
    mag1 = mag(1:window/2) + mag(window/2+1:window);
    [M,~] = max(mag1);
    avg = mean(mag1);
    outlier = M/avg;
end
%% 使用点积检测是否有upchirp
function [max_correlation,index] = detect_upchirp(samples,upchirp_ifreq,window)
    samples_ifreq = instantaneous_frequency(samples,window*2);
    max_correlation = 0;
    index = 1;
    for i=1:1:window
        A=samples_ifreq(i:i+window-1);
        corr = abs(corrcoef(samples_ifreq(i:i+window-1), upchirp_ifreq));
        max_corr = corr(2);
        if max_corr > max_correlation
            max_correlation = max_corr;
            index = i;
        end
    end
end
%% 滑动窗口检测chirp（标准差）
function [max_correlation,index] = detect_upchirp1(samples1,upchirp_ifreq,window)
    samples_ifreq = instantaneous_frequency(samples1,window*2);
    max_correlation = 0;
    index = 1;
    for k=1:1:window
        max_corr = cross_correlate_ifreq(samples_ifreq(k:k+window),upchirp_ifreq,window);    
        if max_corr > max_correlation
            max_correlation = max_corr;
            index = k;
        end
    end
end
%% 检测upchirp 乘上downchirp 做fft找集中情况
function [max_correlation,index] = detect_upchirp2(samples1,ideal_chirp,window)
    samples = samples1(1:window*2);
    max_correlation = 0;
    index = 1;
    for k=1:1:window
        [~,max_corr] = mult_fft(samples(k:k+window),ideal_chirp,window);    
        if max_corr > max_correlation
            max_correlation = max_corr;
            index = k;
        end
    end
end
%% 使用标准差求相识度
function result = cross_correlate_ifreq(samples_ifreq,ideal_chirp,to_idx)
    result = 0;
    average = mean(samples_ifreq(1:to_idx));
    chirp_avg = mean(ideal_chirp);
    sd = std(samples_ifreq) * std(ideal_chirp);
    for i=1:1:to_idx
        result = result + ( samples_ifreq(i) - average ) * (ideal_chirp(i) - chirp_avg) / sd;
    end
    result = result / to_idx;
end
            
        




