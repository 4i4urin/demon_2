% продолжение номера 13
% на вход сюда передаточная функция из графика
% заебало

syms v;
syms z;
syms s;
syms x;

% Журов

w = zeros(1, 8);
disp(['Введите точки абсциссы точек перегиба полученного графика слева' ...
    ' направо (их должно быть 8):']);
for i = 1:8
    w(i) = 1 / input("Точка " + i + ": ");
end

Wky_s = (w(2) * s + 1) * (w(3) * s + 1) * (w(4) * s  + 1) * (w(5) * s + 1) / +...
    ((w(1) * s + 1) * (w(6) * s + 1) * (w(7) * s + 1) * (w(8) * s + 1));

T = input("Введиете дискретность графика T = ");

Ws = input("Введите передаточную функцию разомкнутой системы: ");
[n, d] = numden(Ws);
Ws_tf = tf(sym2poly(n), sym2poly(d));
display(Ws_tf);

Wky = simplify(subs(Wky_s, s, (z - 1) / (z + 1)));
display(Wky);
[num, den] = numden(Wky);
Wky_z = tf(sym2poly(num), sym2poly(den));
disp("После подстановки w = (1 - z) / (z + 1)");
disp("В отчёте s заменить на Z");
display(zpk(Wky_z));


% display(d2c(Wky_z));

disp("Реализацив в виде Z-формы");
Z_ky = filt(sym2poly(num), sym2poly(den), T);
display(Z_ky);

disp("По Wку(z) получаем разностное уравнение");
num = sym2poly(num);
den = sym2poly(den);
num = num / den(1);
den = den / den(1);
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

% пункт 14
% Даааа ебать того рот

Wky_fromz = d2c(Z_ky, 'tustin');
display(Wky_fromz);
W = Ws_tf * Wky_fromz;
display(W);
W = W / (1 + W);
Wd = c2d(W, T);
step(Wd);
grid on;
hold on;
S = stepinfo(Wd);
disp(S);



