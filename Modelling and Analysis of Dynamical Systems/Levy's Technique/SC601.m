%% Generating Datapoints
s = tf("s");
G = (5 + 7*s + s^2)/( s^7 + 12*s^6 + 16*s^5 + 24*s^4 + 20*s^3 + 15*s^2 + 23*s + 1);
[R,I,F] = nyquist(G);
R = squeeze(R);
I = squeeze(I);
data = [F,R,I];
data = data(41:120,:);  % Significant Frequency Range
m = size(data,1);
omega = data(:,1);
R = data(:,2);
I = data(:,3);

%% Calculating the elements of matrix

L0 = sum(omega.^(0));
L2 = sum(omega.^(2));
L4 = sum(omega.^(4));

S0 = sum((omega.^0).*(R));
S2 = sum((omega.^2).*(R));
S4 = sum((omega.^4).*(R));
S6 = sum((omega.^6).*(R));
S8 = sum((omega.^8).*(R));

T1 = sum((omega.^1).*(I));
T3 = sum((omega.^3).*(I));
T5 = sum((omega.^5).*(I));
T7 = sum((omega.^7).*(I));
T9 = sum((omega.^9).*(I));

U2 = sum( (omega.^2) .* ((R.^2) + (I.^2)));
U4 = sum( (omega.^4) .* ((R.^2) + (I.^2)));
U6 = sum( (omega.^6) .* ((R.^2) + (I.^2)));
U8 = sum( (omega.^8) .* ((R.^2) + (I.^2)));
U10 = sum( (omega.^10) .* ((R.^2) + (I.^2)));
U12 = sum( (omega.^12) .* ((R.^2) + (I.^2)));
U14 = sum( (omega.^14) .* ((R.^2) + (I.^2)));

%% Solving the matrix equation

A = [L0 , 0, -L2 , T1 , S2 , -T3 , -S4 , T5 , S6 , -T7;
      0 , L2, 0 , -S2 , T3 , S4 , -T5 , -S6 , T7, S8;
      L2, 0, -L4, T3, S4, -T5, -S6, T1, S8, -T9;
      T1 , -S2, -T3 , U2 , 0 , -U4 , 0 , U6 , 0 ,-U8;
      S2 , T3, -S4 , 0 , U4 , 0 , -U6 , 0 , U8 ,0;
      T3 , -S4, -T5 , U4 , 0 , -U6 ,0 , U8 , 0 , -U10;
      S4 , T5, -S6 , 0 , U6 , 0 , -U8 , 0 , U10, 0;
      T5 , -S6, -T7 , U6 , 0 , -U8 , 0 , U10 ,0, -U12;
      S6 , T7, -S8 , 0 , U8 , 0 , -U10 , 0 , U12, 0;
      T7 , -S8, -T9, U8, 0, -U10, 0, U12, 0, -U14 ];

B = [S0;
     T1;
     S2;
     0;
     U2;
     0;
     U4;
     0;
     U6;
     0];

%params = inv(A)*B;
params = inv(A'*A)*A'*B;   % contains b0,b1,b2,a1,a2,a3,a4,a5,a6,a7 in that order

%% Comparing obtained tf with original

G_obt = (4.9246 -5.3638*s -3.8713*s^2)/(-4.97*s^7 - 16.0965*s^6 - 0.0204*s^5 - 28.5742*s^4 + 22.037*s^3 - 34.7882*s^2 + 18.95*s + 1);

% Nyquist plot are almost identical for Original TF and Obtained TF
nyquist(G,G_obt)
legend('G_Original','G_Obtained')
% Now generating datapoints for frequency range 0.01 to 4 rad/sec

w = 0; % Frequency variable
err_vector = zeros(80,1);
i = 1;
for w = 0.01:0.05:4
    G_jw = ((5 - 1*w^2) + (w*7i)) / ( (1 - 15*w^2 + 24*w^4 -12*w^6) + w*(23 - 20*w^2 + 16*w^4 - 1*w^6)*i);  % Original TF
    G_obt_jw = ((4.9246 - (-3.8713)*w^2) + (w*(-5.3638)*1i)) / ( (1 - (-34.7882)*w^2 + (-28.5742)*w^4 -(-16.0965)*w^6) + w*(18.95 - (22.037)*w^2 + (-0.0204)*w^4 - (-4.97)*w^6)*i);  % Obtained TF
    err_vector(i,1) = (G_jw - G_obt_jw);
    i = i+1;
end
error = sum(abs(err_vector))/80;
fprintf("The error is ")
disp(error)
fprintf("which is a very small value. This is also confirmed from Nyquist Plot shown.")