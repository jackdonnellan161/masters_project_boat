dt = 1/200;
t = 0:dt:30;
accel_data_clean = cos(2*pi*t);
noisevec = 1*randn(1,length(t));
accel_data_noisy = accel_data_clean + noisevec;
figure
subplot(211)
plot(t,accel_data_clean)
subplot(212)
plot(t,accel_data_noisy)

% x = [p v]'
p0 = 0;
v0 = 0;
x_clean = zeros(2,length(t));
x_clean(:,1) = [p0;v0];
x_noisy = zeros(2,length(t));
x_noisy(:,1) = [p0;v0];
A = [1 dt;0 1];
B = [0.5*dt^2;dt];
for ii = 2:length(t)
    x_noisy(:,ii) = A*x_noisy(:,ii-1)+B*accel_data_noisy(ii-1);
    x_clean(:,ii) = A*x_clean(:,ii-1)+B*accel_data_clean(ii-1);
end
figure
subplot(211)
plot(t,x_clean(1,:),t,x_clean(2,:))
subplot(212)
plot(t,x_noisy(1,:),t,x_noisy(2,:))