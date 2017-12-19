%FT:
%syms t;
%y =(heaviside(t)-heaviside(t-5))*(cos(2 * pi * 5 * t) + 2 * sin(2 * pi * 15 * t)) ...
%+ (heaviside(t-5)-heaviside(t-10))*(cos(2 * pi * 20 * t)) ...
%+ (heaviside(t-10)-heaviside(t-15))*(cos(2 * pi * 30 * t) + 0.6 * sin(2 * pi * 45 * t)) ...
%+ (heaviside(t-15)-heaviside(t-20))*sin(2 * pi * 50 * t);
%fy=fourier(y)
%ezplot(fy*conj(fy),[-2*pi*60,2*pi*60,0,10])
%title('FT(x(t))*conj(FT(x(t)))')

%STFT:
t = 0:0.001:25;
F = -60:0.01:60;
x = zeros('like',t)
for i=1:length(t)
   u = t(i);
   x(i) =(heaviside(u)-heaviside(u-5))*(cos(2 * pi * 5 * u) + 2 * sin(2 * pi * 15 * u)) ...
   + (heaviside(u-5)-heaviside(u-10))*(cos(2 * pi * 20 * u)) ...
   + (heaviside(u-10)-heaviside(u-15))*(cos(2 * pi * 30 * u) + 0.6 * sin(2 * pi * 45 * u)) ...
   + (heaviside(u-15)-heaviside(u-20))*sin(2 * pi * 50 * u);
end
id = 1:1:512
window = zeros('like',id)
delta = 0.01
for i=1:length(id)
   tau = (id(i) - 256.0) * 0.001
   window(i) = exp(-(tau^2) / (2 * delta^2));
   wi = window(i);
end
[S,F,T,P] = spectrogram(x,window,256,F,1000);
surf(T,F,abs(S),'edgecolor','none'); axis tight;
%view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
title('delta = 0.01')


%Real FT
%syms t;
%y =(heaviside(t)-heaviside(t-5))*(cos(2 * pi * 5 * t) + 2 * sin(2 * pi * 15 * t)) ...
%+ (heaviside(t-5)-heaviside(t-10))*(cos(2 * pi * 20 * t)) ...
%+ (heaviside(t-10)-heaviside(t-15))*(cos(2 * pi * 30 * t) + 0.6 * sin(2 * pi * 45 * t)) ...
%+ (heaviside(t-15)-heaviside(t-20))*sin(2 * pi * 50 * t);

%y = (heaviside(t-0.1) - heaviside(t-1)) * y;
%fy=fourier(y)
%ezplot(fy*conj(fy),[-2*pi*60,2*pi*60,0,0.8])
%title('FT(x(t))*conj(FT(x(t))) 100ms to 1000ms')

%Gabor:

% id = 1:1:2500
% xx = zeros('like',id)
% for j = 1:2500
%     %syms t
%     m = 0;
%     t = j * (20 ) / 2500 - m;
%     s = 0.1;
%     xx(j) = ((heaviside(t + m)-heaviside(t-5 + m))*(cos(2 * pi * 5 * (t + m)) + 2 * sin(2 * pi * 15 * (t + m))) ...
%            + (heaviside(t-5 + m)-heaviside(t-10 + m))*(cos(2 * pi * 20 * (t + m))) ...
%            + (heaviside(t-10 + m)-heaviside(t-15 + m))*(cos(2 * pi * 30 * (t + m)) + 0.6 * sin(2 * pi * 45 * (t + m))) ...
%            + (heaviside(t-15 + m)-heaviside(t-20 + m))*sin(2 * pi * 50 * (t + m)));
% end

% % M = zeros('like',id);
% % F = zeros('like',id);
% M = 0:0.1:20;
% %F = 1/2500:1/2500:1
% F = -0.5:1/2500:0.5;
% F = F * 100;
% S = zeros(length(F),length(M));
% l = 0
% for i = 1:201
%     m = (i-1) * 0.1
%     %m = 0
%     x = zeros('like',xx);
%     for j = 1:2500
%         t = j * (20) / 2500 - m;
%         s = 2;
%         x(j) = xx(j) * exp(-0.5 * (t)^2 / s^2);
%     end
%     fy = fft(x,2500);
%     gy = (fy + conj(fy)) * 0.5;
%     for j = 1:2500
%         k = j + 1250;
%         if k > 2500
%             k = k - 2500;
%         end
%         S(k,i) = gy(j);
%     end
%     %plot(abs(gy));
% end
% [M,F] = meshgrid(M,F);
% surf(M,F,abs(S),'edgecolor','none'); axis tight;
% %view(0,90)
% xlabel('Time (Seconds)'); ylabel('HZ');
% %ezmesh(gabor)