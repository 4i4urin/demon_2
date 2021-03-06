% 16 пункт 
% 15 прошёл как по маслу
% попробуем 16 бахнуть
% и 17 за одно
% коэффициенты подгонки зарешали

syms s;
syms x;

% Журов
% Ws = 562.5 * 8750 / ((0.026 * s + 1) * (0.0055 * s + 1) * s);
% Wky_s = (s * 66.667 + 1) * (s * 47.62 + 1) * (s * 10 + 1) * (s * 2 + 1) / +...
%     ((s * 3703703.703 + 1) * (s + 1) * (s * 0.597 + 1) * (s * 0.5571 + 1));
% из 15 пунта
Ws = 221 / ((0.026 * s + 1) * (0.0055 * s + 1) * s);

% Чел
% Ws = 528.125 * 9375 / (s * (0.011 * s + 1) * (0.022 * s + 1));

[num, den] = numden(Ws);
Ws_tf = tf(sym2poly(num), sym2poly(den));
[A, B, C, D] = ssdata(Ws_tf);

disp("A = ");
disp(A);
disp("B = ");
disp(B);
disp("C = ");
disp(C);
disp("D = ");
disp(D);

% изменяя коэфициент подгонки меняешь графки
% коэффициент подгонки 200 
P = [-1, -0.5 + 0.866i, -0.5 - 0.866i] * 200;
K = place(A, B, P);
disp("А теперь вычислим матрицу коэффициентов обратных связей");
disp("K = ");
disp(K);
disp("Вычислим матрицу обратной связи наблюдателя");
% коэффициент подгонки 200 
pn = [-1, -0.5 + 0.866i, -0.5 - 0.866i] * 200;
L = place(A', C', pn)';
display(L);
disp("вычислим матрицы динамического регулятора для расчета передаточной функции");
[Ar, Br, Cr, Dr] = reg(A, B, C, D, K, L);
disp("Ar = ");
disp(Ar);
disp("Br = ");
disp(Br);
disp("Cr = ");
disp(Cr);
disp("Dr = ");
disp(Dr);
[numr, denr] = ss2tf(Ar, Br, Cr, Dr);
Wreg = tf(numr, denr);
display(Wreg);
Wz = feedback(Ws_tf, Wreg);
display(Wz);
step(Wz);
grid on;
S = stepinfo(Wz);
disp(S);
