addpath('../../utils');
id_list_cell = {'A'; 'A'; 'B'; 'C'; 'C'; 'D'};
y = [1; 0; 0; 1; 1; 0];
n_folds = 2;
try
    fold_id = make_grouped_folds(id_list_cell, y, n_folds);
    disp(fold_id);
catch e
    disp(e.message);
end
