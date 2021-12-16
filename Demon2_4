% продолжение номера 13
% на вход сюда передаточная функция из графика
% заебало

syms v;
syms z;
syms s;
syms x;

% Журов

Wky_s = (s * 66.667 + 1) * (s * 47.62 + 1) * (s * 10 + 1) * (s * 2 + 1) / +...
    ((s * 3703703.703 + 1) * (s + 1) * (s * 0.597 + 1) * (s * 0.5571 + 1));
Ws = 562.5 * 8750 / ((0.026 * s + 1) * (0.0055 * s + 1) * s);
T = 0.001309;
% v_mid = 0.083;
% Wky = (0.033 * s + 1)*(0.026 * s + 1)*(0.0055 * s + 1) /+...
%     ((21.8 * s + 1)*(0.00042 * s + 1)*(0.00042 * s + 1));

%Чел
% Wky_s = (1 / 0.023 * s + 1) * (1 / 0.023 * s + 1) * (1 / 0.045 * s + 1) / +...
%     ((1 / 0.0000007 * s + 1) * (s + 1) * (1 / 1.698 * s + 1) * (1 / 1.767 * s + 1));
% Ws = 528.125 * 9375 / (s * (0.011 * s + 1) * (0.022 * s + 1));
% T = 0.002618;




% Демон
% Wky_s = (1 / 0.033 * s) * (2.16 * s + 1) / ((1 / 0.005 * s + 1) * (s + 1));
% Ws = 40 / (s * (0.1 * s + 1));
% T = 0.1;

[n, d] = numden(Ws);
Ws_tf = tf(sym2poly(n), sym2poly(d));
display(Ws_tf);

Wky = simplify(subs(Wky_s, s, (z - 1) / (z + 1)));
display(Wky);
[num, den] = numden(Wky);
Wky_z = tf(sym2poly(num), sym2poly(den));
display(zpk(Wky_z));


% display(d2c(Wky_z));

disp("Реализацив в виде Z-формы");
Z_ky = filt(sym2poly(num), sym2poly(den), T);
display(Z_ky);

% пункт 14
% Даааа ебать того рот

Wky_fromz = d2c(Z_ky, 'tustin');
display(Wky_fromz);
W = Ws_tf * Wky_fromz;
display(W);
W = W / (1 + W);
Wd = c2d(W, T);
step(Wd);



% разностное уровнение рот ебал 
% может правильно а может и нет
% по всей видимотси не тем нум ден
num = sym2poly(num);
den = sym2poly(den);
% display(num); display(den);
s = "";
for i = 1:max(size(num))
    s = s + num2str(num(i)) + " * e[(k - " + num2str(i - 1) + ")T] ";
    if i ~= max(size(num))
        s = s + "+ ";
    end
end
s = s + " = ";
for i = 1:max(size(den))
    s = s + num2str(den(i)) + " * e[(k - " + num2str(i - 1) + ")T] ";
    if i ~= max(size(den))
        s = s + "+ ";
    end
end
disp(s);
