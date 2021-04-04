close all
clc
% read audio files and sample them 
[data, fs] = audioread('C:\Users\Ashish\Desktop\noisy_speech_x.wav');
[data1, fs] = audioread('C:\Users\Ashish\Desktop\desired_signal.wav');
data = resample(data,16000,fs);
fram = buffer(data,124,0);
data1 = resample(data1, 16000, fs);
fram1 = buffer(data1,124,0);

% wiener filter 
 
  for j=1:403
% %       
  yxx = xcorr(fram(:,j),fram(:,j)); 
  yxx=yxx(length(yxx)/2:length(yxx));
  ydx = xcorr(fram1(:,j),fram(:,j));
  ydx=ydx(length(ydx)/2:length(ydx));
  ydd = xcorr(fram1(:,j),fram1(:,j));
   ydd=ydd(length(ydd)/2:length(ydd));

  yxx_row = yxx.';
  ydx_row = ydx.';

    T = toeplitz(yxx_row);
    Tinv = inv(T);
   
    
    Hopt= Tinv*ydx; 
    yn = conv(Hopt,fram(:,j)); 
  yx(:,j)=yn;
  
    n=length(yn);
    f=linspace(-fs/2, fs/2,n);
   ynw=fft(yn,n)/n;
  end 
      
  
ash=yx'
ann=ash(:,1) ;
 for i = 3:247
 ann = cat(1,ann,ash(:,i));
 end
 
 play=ann(1:49972); % result from the Wiener Filter (pure Wiener filtering)
 mmse_wiener= immse(play,data1);
 
 figure(1);
   plot(f,abs(fftshift(ynw)));
   title('freq resp of pure wiener');
 
 
 % NLMS FOLLWED BY WIENER FILTER (tHIS IS FOR BEST ACOUSTIC REGENERATION)
 % FOR BEST MMSE RESULT PLEASE CHANGE THE VALUE 256 TO 5, AND 0.5 TO 0.002627
 
 nlms = dsp.LMSFilter(256, 'StepSize',0.5, 'Method', 'Normalized LMS');%stepsize taken is 2/lambdamax
  [yplay, e, h] = step(nlms, play,data1);
  n1=length(yplay);
    f1=linspace(-fs/2, fs/2,n1);
   ynw1=fft(yplay,n1)/n1;
   figure(2);
   plot(f1,abs(fftshift(ynw1)));
   title('freq. resp wiener followed by nlms')
   
   figure(3);
   plot(yplay);
   title('time domain wave form of output of Winer+ nlms')
   
  
  lambdamax= xcorr(data,data);
  lambdamax=lambdamax(length(lambdamax)/2:length(lambdamax));  % notice that the first element is 761.315. this is the largest eigen value.
  
  myresult=yplay*3 ;
   MMSE=  min(((yplay-data1).*(yplay-data1))/49972);
   msematlab = immse(yplay,data1);
   
   
   
   % PURE NLMS :
   
   nlms3 = dsp.LMSFilter(100, 'StepSize', 0.43, 'Method', 'Normalized LMS');
   [yn3, e3, h3] = step(nlms3, data,data1);
    
  npure=length(yn3);
    fpure=linspace(-fs/2, fs/2,npure);
   ynwpure=fft(yn3,n1)/n1;
   figure(4);
   plot(fpure,abs(fftshift(ynwpure)));
   title('freq resp of pure nlms');
   
   figure(5);
   plot(ann);
   title('time domain output of wiener');
   
   figure(6);
   plot(data);
   title('time domain : input');
   
   figure(7);
   plot(data1);
   title('time domain : desired');
   
   figure(8);
   ndata=length(data);
    fdata=linspace(-fs/2, fs/2,ndata);
   ynwdata=fft(data,ndata)/ndata;
   plot(fdata,abs(fftshift(ynwdata)));
   title('freq domain : input');
   
  figure(9);
   ndata1=length(data1);
    fdata1=linspace(-fs/2, fs/2,ndata1);
   ynwdata1=fft(data1,ndata1)/ndata1;
   plot(fdata1,abs(fftshift(ynwdata1)));
   title('freq domain : desired');
   
   figure(10);
   plot(yn3);
   title('Time Domain: Pure NLMS')
   
   
   