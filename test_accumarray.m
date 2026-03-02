ic = [1; 1; 2; 2; 3; 3];
y = [0; 0; 0; 1; 1; 1];
n_unique = 3;
pt_y = double(accumarray(ic, double(y > 0), [n_unique, 1], @any));
disp(pt_y);
