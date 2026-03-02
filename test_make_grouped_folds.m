addpath('utils');
id_list_cell = {'PT1'; 'PT1'; 'PT2'; 'PT3'; 'PT3'; 'PT4'};
y = [0; 1; 0; 1; 1; 0];
n_folds = 2;

fold_id = make_grouped_folds(id_list_cell, y, n_folds);
disp(fold_id);
