syms s;

% Денчик
% sigma = 5;
% Tmax = 0.026;
% Emax = 1.3;
% Kg = 62500;

% Шиша
% sigma = 7;
% Tmax = 0.028;
% Emax = 1.2;
% Kg = 62500;
% Ws = 160 / (s * (0.00035 * s^2 + 0.0405 * s + 1));
sigma = input("Введите сигму: ");
Tmax = input("Введите Tпмакс: ");
Emax = input("Введите Eмакс: ");
disp("Введите передаточную функцию разомкнутой системы");
disp("Пример 160 / (s * (0.00035 * s^2 + 0.0405 * s + 1))")
Ws = input("Введите передаточную функцию разомкнутой системы: ");


% Wky = (0.035 * s + 1)(0.02777 * s + 1)(0.01695 * s + 1) / +...
% (31.361 * s + 1)(0.0008615 * s + 1)(0.0004498 * s + 1);
% W = tf([1.35565, 167.53, 6559.98, 82287.8], [1, 3384.01, 2.58074*10^4, 82287.8]);
% display(W);
% bode(W,1e-2:0.1:1e4); % ЛАФЧХ
% grid on

% Серба
% sigma = 5;
% Tmax = 0.030;
% Emax = 2.0;
% Kg = 62500;
% Ws = 1083000 / (s^3 + 144 * s^2 + 3096 * s);


% Рома
% sigma = 7;
% Tmax = 0.028;
% Emax = 1.9;
% Kg = 62500;
% Ws = 73.1 / (3.6e-4 * s ^ 3 + 0.049 * s ^ 2 + s);

% Журов
% sigma = 7;
% Tmax = 0.038;
% Emax = 1.4;
% Kg = 62500;
% Ws = 125 / ((0.026 * s + 1) * (0.0055 * s + 1) * s);

[num, den] = numden(Ws);
display(num); display(den);
cf = sym2poly(den);
num = num / cf(3);
den = expand(den / cf(3));


B = containers.Map('KeyType', 'double', 'ValueType', 'double');
B(5) = 6.5; B(10) = 6.7; B(20) = 6.9; B(25) = 8.8; B(30) = 11.3; 
B(35) = 14.1; B(40) = 16.9;

k = cell2mat(keys(B));
b = 0;
for i = 1:length(B)
    key = k(i);
    if (key > sigma)
        break;
    end
    b = B(key);
end


w_mid = b / Tmax;
w_begin = 0.01;
w_low = 0.16 * w_mid;
w_high = 6.5 * w_mid;
w_end = 10000;


disp("Из таблици вычисляем B ");
disp(b);
disp("w_mid = B / Tпмакс");
disp("Средняя частота ");
disp(w_mid);
disp("Нижняя частота ");
disp(w_low);
disp("Высока частота ");
disp(w_high);

Ktr = 1 / (Emax * 62500 * 1e-6);
disp("Ктр = ");
disp(Ktr);


w_mid_y = 0;
w_low_y = w_mid_y + 20 * abs(log10(w_low) - log10(w_mid));
w_begin_y = w_low_y + 40 * abs(log10(w_begin) - log10(w_low));
w_high_y = w_mid_y - 20 * abs(log10(w_mid) - log10(w_high));
w_end_y = w_high_y - 40 * abs(log10(w_high) - log10(w_end));

x_lower = [w_begin, w_low, w_mid, w_high, w_end];
y_lower = [w_begin_y, w_low_y, w_mid_y, w_high_y, w_end_y];


K = round(num, 4);


den = factor(den, 'FactorMode', 'real');
dend1 = round(coeffs(den(3)), 4); dend2 = round(coeffs(den(4)), 4);

T1 = min(dend1(1), dend1(2)) / max(dend1(1), dend1(2));
T2 = min(dend2(1), dend2(2)) / max(dend2(1), dend2(2));

wc1 = round(min(1 / T1, 1 / T2));
wc2 = round(max(1 / T2, 1 / T1));
wc_k = 1;

disp("Твой коефицент усиления ");
disp(K);
disp("Графики должны пересекаться. При том зелёный график должен быть ниже");
disp("Если графики не пересекаются перезапусти прогу и измени коефициент усиления");
disp("При первом запуске нажми n");

choice = input("Увеличить коэфициент усиления? [y/n]: ", 's');
    if (ischar(choice) && lower(choice) == 'y')
        disp("Совет: изменяй коефициент в 5, 10, 100 раз");
        K = input('НОВЫЙ коэффициент усиления K: ');
    end
disp("20log(K) = ");
l = 20 * log10(K);
disp(l);
% k_kala = 1e5 / Ktr;
% display(k_kala);

wc_k_y = 20 * log10(K);
w_begin_y = wc_k_y + 20 * abs(log10(wc_k) - log10(w_begin));
wc1_y = wc_k_y - 20 * abs(log10(wc_k) - log10(wc1));
wc2_y = wc1_y - 40 * abs(log10(wc1) - log10(wc2));
w_end_y = wc2_y - 60 * abs(log10(wc2) - log10(w_end));

x_higher = [w_begin, wc_k, wc1, wc2, w_end];
y_higher = [w_begin_y, wc_k_y, wc1_y, wc2_y, w_end_y];

disp("Точки перегиба графиков");
disp("Точка 2 ");
disp(wc1);
disp("Точка 3 ");
disp(w_low);
disp("Точка 4 ");
disp(wc2);
disp("Точка 5 ");
disp(w_high);

semilogx(x_lower, y_lower, 'g');
grid on
hold on
semilogx(x_higher, y_higher);
