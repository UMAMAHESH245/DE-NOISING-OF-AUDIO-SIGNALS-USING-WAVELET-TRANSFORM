clc 
clear all
close all
%------->> Stage 1:
% Online recording via microphone.
prompt = {'Enter the number of second to be recorded','Enter the frequency sampling rate'}; 
dlg_title = 'Audio Recording'; 
num_lin = 2; 
def_t_f = {'5','5000'}; 
record_par = inputdlg(prompt ,dlg_title ,num_lin,def_t_f); 
record_t = str2num(record_par{1}); 
Fs = str2num(record_par{2}); 
record_speech = audiorecorder(Fs,16,1); 
msg0 = msgbox('Start speaking'); 
h = waitbar(0,'Start speaking') 
for step=1:1 
 recordblocking(record_speech , record_t); 
 waitbar(step) 
end
msg1 = msgbox('End speaking'); 
pause(2); 
close(h) 
%SNR to be maintained
snrval=10; 
noisy_speech=getaudiodata(record_speech); 
fig1=figure('name','noisy recorded speech','color','w'); 
%plotting recorded sppech
plot((1:length(noisy_speech))/Fs,noisy_speech); 
lev_wname = inputdlg({'Enter the number of decomposition levels','Enter the wavename : 
dbN'},'Number of level and wavename'); 
level = str2num(lev_wname{1}); 
wavename = lev_wname{2}; 
% Find the wavelet and scaling functions
iteration = 8; 
[phi,psi,xval_dbn] = wavefun(wavename,iteration); 
% Find the wavelet and scaling filters
[Lo_d , Lo_r] = wfilters(wavename , 'l'); 
Hi_r = qmf(Lo_r); 
Hi_d = wrev(Hi_r); 
sum_Hi_d = sum(Hi_d); 
sum_Lo_d = sum(Lo_d); 
Nextpow2_Lo_d = nextpow2(length(Lo_d)); 
Nextpow2_Hi_d = nextpow2(length(Lo_d)); 
Nextpow2_Lo_r = nextpow2(length(Lo_r)); 
Nextpow2_Hi_r = nextpow2(length(Lo_r)); 
fftLo_d = fft(Lo_d,2^Nextpow2_Lo_d); 
fftHi_d= fft(Hi_d,2^Nextpow2_Hi_d); 
fftLo_r = fft(Lo_r,2^Nextpow2_Lo_r); 
fftHi_r = fft(Hi_r,2^Nextpow2_Hi_r); 
% Decompose the noisy signal at a given level using the wavelet filters
[Cad,L] = wavedec(noisy_speech(:,1) ,level ,wavename); 
% Extract the approximation coefficients
% Extract the detail coefficients
Capp = appcoef(Cad , L , Lo_d , Hi_d, level);%Extract the approximation coefficients for 
level (level)
a2sq = sum(Capp.^2);%energy of the wavelet approximation coefficients
Energy_coef = zeros(1,level+1); 
Energy_coef(1) = a2sq; 
SDEV = zeros(1,level); 
SDEV_COEF = zeros(1,level); 
for i= level : -1 : 1 
 Cdet = detcoef(Cad,L,i);% Extract the detail coefficients for level (i)
 d2sq = sum(Cdet.^2); 
 Array_det_coef{level-i+1} = Cdet; 
 d2sq = sum(Array_det_coef{level-i+1}.^2);%energy of the wavelet coefficients for every 
scale
 Energy_coef(level-i+2)= d2sq; 
 SDEV(i) = wnoisest(Cad , L ,i);% standard deviation approximation for detail 
coefficients for every scale
 SDEV_COEF(i) = std(Cdet); 
end
engcoef = sum(Energy_coef);% total energy of wavelet coefficients(approximation & detail)
engsig = sum(noisy_speech(:,1).^2);% total energy of the noisy signal
engerr = abs(engcoef - engsig);% energy preservation 
percent_energy_app=(Energy_coef(1)/sum(Energy_coef))*100; 
% Reconstruct the approximation signal at coarser scale(final level)from wavelet 
decomposition
App_sig = wrcoef('a',Cad,L,wavename,level); 
ReconstructArray_sig{1} = App_sig; 
% Reconstruct the detail signals at all levels from the wavelet decomposition
SDEV_REC = zeros(1,level); 
for j = level : -1 : 1 
 Det_sig = wrcoef('d' , Cad , L , wavename , j); 
 ReconstructArray_sig{level-j+2} = Det_sig; 
 SDEV_REC(j) = std(Det_sig); 
end
% Plotting illustration
fig3 = figure('name' , 'phi_psi functions ,filters and fft of filters','Color','w'); 
subplot(5,2,1); axis tight ; 
plot(xval_dbn, phi);title('Scaling function phi'); 
subplot(5,2,2); axis tight ; 
plot(xval_dbn, psi);title('Wavelet function psi'); 
subplot(5,2,3); axis tight ; 
stem(Lo_d,'r');title('Decomposition low pass filter'); 
subplot(5,2,4); axis tight ; 
stem(Hi_d,'r');title('Decomposition high pass filter'); 
subplot(5,2,5); axis tight ; 
stem(Lo_r,'b');title('Reconstruction low pass filter'); 
subplot(5,2,6); axis tight ; 
stem(Hi_r,'b');title('Reconstruction high pass filter'); 
subplot(5,2,7); axis tight
Freq_Lo_d =(2*pi)/(2^Nextpow2_Lo_d):(2*pi)/(2^Nextpow2_Lo_d): pi ; 
fa_ld = abs(fftLo_d(1:(2^(Nextpow2_Lo_d)/2))); 
plot(Freq_Lo_d,fa_ld);title('FFT of analysis low pass filter'); 
xlim([(2*pi)/(2^Nextpow2_Lo_d) , pi+0.5]); 
subplot(5,2,8); axis tight
Freq_Lo_r = (2*pi)/(2^Nextpow2_Lo_r):(2*pi)/(2^Nextpow2_Lo_r) :pi ; 
fa_hr = abs(fftLo_r(1:(2^(Nextpow2_Lo_r)/2))); 
plot(Freq_Lo_r,fa_hr);title('FFT of synthesis low pass filter'); 
xlim([(2*pi)/(2^Nextpow2_Lo_r) , pi+0.5]); 
subplot(5,2,9); axis tight
Freq_Hi_d =(2*pi)/(2^Nextpow2_Hi_d):(2*pi)/(2^Nextpow2_Hi_d) :pi ; 
fa_ld = abs(fftHi_d(1:((2^Nextpow2_Hi_d)/2))); 
plot(Freq_Hi_d,fa_ld);title('FFT of analysis high pass filter'); 
xlim([(2*pi)/(2^Nextpow2_Hi_d) , pi+0.5]); 
subplot(5,2,10); axis tight 
Freq_Hi_r =(2*pi)/(2^Nextpow2_Hi_r):(2*pi)/(2^Nextpow2_Hi_r) : pi ; 
fa_hr = abs(fftHi_r(1:((2^Nextpow2_Hi_r)/2))); 
plot(Freq_Hi_r,fa_hr);title('FFT of synthesis high pass filter'); 
xlim([(2*pi)/(2^Nextpow2_Hi_r) , pi+0.5]); 
%----------
fig4 = figure('name' , 'decomposition coefficients','Color','w'); 
%----------
subplot(level+2,1,1); axis tight; 
plot([1:length(noisy_speech)]/Fs , noisy_speech,'k'); 
title('Noisy speech signal'); 
subplot(level+2,1,2); axis tight; 
plot(Capp,'b') 
title('clear speech signal'); 
ylabel(['Ca',num2str(level)],'Color','b'); 
for f = level : -1 : 1 
 row = level+2; 
 no_fig = level-f+3; 
 s = level-f+1; 
 lbl = num2str(f); 
 subplot(row,1,no_fig);axis tight; 
 plot(Array_det_coef{s},'r') 
 title(['Detail coefficients at level',lbl]); 
 ylabel(['Cd',lbl],'Color','r') 
end
subplot(row,1,row) 
%----------
fig5 = figure('name' , 'reconstructed coefficients','Color','w') ; 
%----------
subplot(level+2,1,1); axis tight; 
plot([1:length(noisy_speech)]/Fs , noisy_speech,'k'); 
title(['Noisy speech signal : a',num2str(level),' + d',num2str(level),' = 
d',num2str(level-1)]); 
subplot(level+2,1,2); axis tight; 
plot( ReconstructArray_sig{1} , 'b') 
title(['Approximation signal at level',num2str(level)]); 
ylabel(['a',num2str(level)],'Color','b') 
for r = level : -1 : 1 
 row = level+2; 
 no_fig = level-r+3; 
 s = level-r+1; 
 lbl = num2str(r); 
 subplot(row,1,no_fig); axis tight; 
 plot(ReconstructArray_sig{level-r+2} ,'r') 
 title(['Reconstructed coefficients at level',lbl]); 
 ylabel(['d',lbl],'Color','r') 
end
subplot(row,1,row); 
xlabel('Time'); 
msg5 = msgbox('Set the values of (threshold value , soft or hard thresholding function , 
KeepApp) manually'); 
pause(3) 
std_glb = median(abs(Array_det_coef{level}))/0.6745; 
thr = std_glb*sqrt(2*log(length(noisy_speech))); 
thr_globstr = num2str(thr); 
def_thr_s_1 = {thr_globstr,'s','1'}; 
thr_sorh_k = inputdlg({'Enter the value of threshold','Enter the type of thresholding 
function soft or hard s or h','Threshold the approximation? 1:no or 0:yes'},'Setting 
parameters',1,def_thr_s_1); 
thr = str2num(thr_sorh_k{1}); 
s_or_h = thr_sorh_k{2}; 
KeepApp = str2num(thr_sorh_k{3}); 
% Choose the type of thresholding
gbl_or_lvd = menu('Type of thresholding','Global thresholding','Level-dependent 
thresholding'); 
msg2 = msgbox('De_noising ...'); 
seg_step = Fs*0.01; 
overlap = Fs*0.01; 
seg_len = seg_step + overlap; 
sp_len = length(noisy_speech); 
Nseg = floor(sp_len/(seg_step))-1; 
window = hamming(seg_len); 
de = hanning(2*overlap - 1)'; 
dewindow = [de(1:overlap) , ones(1,seg_len -2*overlap),de(overlap:end)]'./window; 
switch gbl_or_lvd 
% gbl : global thresholding
case 1 
% Reconstruct the signal from thresholded coefficients
denoised_speech = zeros(sp_len, 1); 
for i = 1 : Nseg 
 sp_Seg(:,i) = noisy_speech((i-1)*seg_step+1 :i*seg_step+overlap); 
 noisy_speechW(:,i) = window.*sp_Seg(:,i); 
 [Cad_seg , L] = wavedec(noisy_speechW(:,i),level,wavename); 
 Cdet_seg = detcoef(Cad_seg , L , 1); 
 sigma_seg = median(abs(Cdet_seg))/0.6745; 
 thr_seg = sigma_seg*sqrt(2*log(length(noisy_speechW(:,i)))); 
 [denoised_seg ,Cad_thr_seg , L_thr_seg ,L2norm_recovery_seg , cmp_score_seg] 
=wdencmp('gbl' ,Cad_seg , L , wavename , level ,thr_seg ,s_or_h , 1); 
 denoised_seg (:,i) = denoised_seg; 
 noisy_speechDe(:,i) = denoised_seg (:,i).*dewindow; 
 denoised_speech((i-1)*seg_step+1 : i*seg_step+overlap) = noisy_speechDe(:,i) + 
denoised_speech((i-1)*seg_step+1: i*seg_step+overlap); 
end
%lvd:level-dependentthresholding
case 2 
msg3 = msgbox('Choosing the thresholds Based on threshold selection rules'); 
pause(3) 
%Choose the noise model
noise_mod_menu=menu('Noise model','Unscaled white noise'); 
noise_model={'one'}; 
SCAL=noise_model(noise_mod_menu); 
f=char(SCAL{1}); 
% Chooose the threshold selection rule
thrrule = menu('Threshold selection rule','sqtwolog' ,'minimaxi'); 
menu_thrrule={'sqtwolog','minimaxi'}; 
ThrSelectRule=menu_thrrule(thrrule); 
tptr =ThrSelectRule{1}; 
denoised_speech = zeros(sp_len,1); 
for i=1:Nseg 
 sp_Seg(:,i)=noisy_speech((i-1)*seg_step+1: i*seg_step+overlap); 
 noisy_speechW(:,i)= window.*sp_Seg(:,i); 
 denoised_seg=wden(noisy_speechW(:,i),tptr,s_or_h,f,level,wavename); 
 denoised_seg(:,i)=denoised_seg; 
 noisy_speechDe(:,i)=denoised_seg(:,i).*dewindow; 
 denoised_speech((i1)*seg_step+1:i*seg_step+overlap)=noisy_speechDe(:,i)+denoised_speech((i1)*seg_step+1:i*seg_step+overlap); 
 THR=zeros(1,level); 
end
end
%Plotting illustration
fig9=figure('name','thresholding function illustration','Color','w'); 
lin_fun = linspace(-0.5 , 0.5 , 100); 
lin_fun_s = wthresh(lin_fun , 's' , thr); 
lin_fun_h = wthresh(lin_fun , 'h' , thr); 
subplot(1,3,1); 
plot(lin_fun,lin_fun,'k'); 
title('Originalfunction'); 
xlabel('coefficients before thresholding'); 
subplot(1,3,2); 
plot(lin_fun,lin_fun_s,'b'); 
title('Soft thresholded function'); 
text(thr , -0.05 , [ strcat('(',num2str(thr),',') 
,strcat(num2str(0),')')],'Color','b','HorizontalAlignment','center') 
subplot(1,3,3); 
plot(lin_fun,lin_fun_h,'r'); 
title('Hard thresholded function'); 
text(thr , -0.05 , [ strcat('(',num2str(thr),',') , strcat(num2str(0),')')],'Color' 
,'r','HorizontalAlignment','center') 
for fig = 2 : 3 
 subplot(1,3,fig) 
 xlabel('coefficients after thresholding'); 
end
fig18 = figure('name' ,'Noisy and Denoised recorded speech signals','Color','w'); 
subplot(2,1,1); axis tight; 
plot([1:length(noisy_speech)]/Fs , noisy_speech,'k'); 
xlabel('Time(s)'); ylabel('Amplitude'); 
title('noisy recorded speech signals'); 
legend('noisy recorded speech'); 
subplot(2,2,2); axis tight; 
plot([1:length(noisy_speech)]/Fs , noisy_speech,'k'); 
xlabel('Time(s)'); ylabel('Amplitude'); 
title('Noisy and De-noised speech signal'); 
hold on; 
plot([1:length(denoised_speech)]/Fs,denoised_speech,'g'); 
legend('noisy speech' , 'denoised speech'); 
fig19 = figure('name' , 'powerdistribution','Color','w'); 
Npow2 = pow2(nextpow2(length(denoised_speech))); 
denoised_speech_pad = fft(denoised_speech , Npow2); 
noisy_speech_pad = fft(noisy_speech , Npow2); 
freq_range = (0:(Npow2 - 1))*(Fs/Npow2); 
power_denoised_speech=denoised_speech_pad.*conj(denoised_speech_pad)/Npow2; 
power_noisy_speech=noisy_speech_pad.*conj(noisy_speech_pad)/Npow2; 
subplot(1,2,1);
axis tight; 
plot(freq_range , power_denoised_speech); 
xlabel('frequency Hz'); 
ylabel('power'); 
title('power distribution of denoised speech signal'); 
subplot(1,2,2);
axis tight; 
plot(freq_range , power_noisy_speech); 
xlabel('frequency Hz'); 
ylabel('power'); 
title('power distribution of noisy recorded speech signal');
