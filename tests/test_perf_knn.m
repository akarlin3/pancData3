n_tr = 500; p = 50; k = 5;
X_tr = randn(n_tr, p);
missing_idx = rand(n_tr, p) < 0.1;
X_tr(missing_idx) = NaN;
pat_id_tr = randi(100, n_tr, 1);
tic;
[X_tr_imp] = knn_impute_train_test(X_tr, [], k, pat_id_tr, []);
t = toc;
fprintf('Current logic: %f\n', t);
