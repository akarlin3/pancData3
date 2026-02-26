function tests = test_landmark_cindex
    tests = functiontests(localfunctions);
end

function test_landmark_left_truncation(testCase)
    %% test_landmark_left_truncation validates the data science logic for 
    % Immortal Time Bias prevention and left-truncation (shifting) 
    % necessary for a Landmark Survival Model (e.g. C-index calculation)

    % 1. Create Mock Time-Dependent Data for 4 patients
    % Patient 1: Event before landmark (at t=2). Expected: Dropped completely (Immortal Time Bias)
    % Patient 2: Censored before landmark (at t=2.5). Expected: Dropped completely
    % Patient 3: Event after landmark (at t=5). Expected: Kept, times shifted
    % Patient 4: Censored after landmark (at t=6). Expected: Kept, times shifted

    landmark_time = 3;

    % Columns: PatID, t_start, t_stop, event, Covariate (X)
    raw_data = [
        1,  0,  2,   1, 10; % Pt 1 fails at 2
        2,  0, 2.5,  0, 20; % Pt 2 censored at 2.5
        3,  0,  2,   0, 31; % Pt 3 interval 1 (ends before landmark)
        3,  2,  5,   1, 32; % Pt 3 interval 2 (fails at 5)
        4,  0,  4,   0, 41; % Pt 4 interval 1 (crosses landmark)
        4,  4,  6,   0, 42  % Pt 4 interval 2 (censored at 6)
    ];

    pat_id  = raw_data(:, 1);
    t_start = raw_data(:, 2);
    t_stop  = raw_data(:, 3);
    event   = raw_data(:, 4);
    X       = raw_data(:, 5);

    %% 2. Implement Landmark Logic
    
    % Step A: Identify valid patients (must be followed up past the landmark_time)
    unique_pats = unique(pat_id);
    valid_pats = [];
    for i = 1:numel(unique_pats)
        pid = unique_pats(i);
        p_idx = (pat_id == pid);
        % The final overall survival time for the patient
        max_time = max(t_stop(p_idx));
        
        % Patients who experience the event or are censored <= landmark_time 
        % must be excluded to prevent Immortal Time Bias.
        if max_time > landmark_time
            valid_pats = [valid_pats; pid];
        end
    end
    
    % Step B: Filter data to valid risk set
    keep_mask = ismember(pat_id, valid_pats);
    LM_pat_id  = pat_id(keep_mask);
    LM_t_start = t_start(keep_mask);
    LM_t_stop  = t_stop(keep_mask);
    LM_event   = event(keep_mask);
    LM_X       = X(keep_mask);
    
    % Step C: Shift time zero to landmark_time (Left-Truncation)
    LM_t_start_shifted = LM_t_start - landmark_time;
    LM_t_stop_shifted  = LM_t_stop  - landmark_time;
    
    % Step D: Drop intervals that fall entirely before or exactly at the landmark
    % (since they ended before our new baseline t=0)
    post_lm_mask = LM_t_stop_shifted > 0;
    
    LM_pat_id = LM_pat_id(post_lm_mask);
    % Cap start time at 0, since overlapping intervals effectively begin at the new baseline 0
    LM_t_start_final = max(0, LM_t_start_shifted(post_lm_mask)); 
    LM_t_stop_final  = LM_t_stop_shifted(post_lm_mask);
    LM_event_final   = LM_event(post_lm_mask);
    LM_X_final       = LM_X(post_lm_mask);

    %% 3. Strict Assertions
    
    % Patient exclusions
    testCase.verifyFalse(ismember(1, LM_pat_id), ...
        'Patient 1 should be dropped due to immortal time bias (event at t=2 <= LM=3).');
    testCase.verifyFalse(ismember(2, LM_pat_id), ...
        'Patient 2 should be dropped due to immortal time bias (censored at t=2.5 <= LM=3).');
    
    % Patient inclusions
    testCase.verifyTrue(ismember(3, LM_pat_id), ...
        'Patient 3 should be kept (event at t=5 > LM=3).');
    testCase.verifyTrue(ismember(4, LM_pat_id), ...
        'Patient 4 should be kept (censored at t=6 > LM=3).');
    
    %% Validate Left-Truncation (Time Shifting) for Patient 3
    idx3 = find(LM_pat_id == 3);
    testCase.verifyEqual(numel(idx3), 1, 'Patient 3 should have 1 interval effectively > 0 post-landmark.');
    
    % The original interval for Pt 3 was [2, 5]. Shifted by 3 -> [-1, 2]. 
    % Left-truncated max(0, -1) -> [0, 2]
    testCase.verifyEqual(LM_t_start_final(idx3), 0, 'Patient 3 start time left-truncated strictly to 0.');
    testCase.verifyEqual(LM_t_stop_final(idx3), 2, 'Patient 3 stop time properly shifted by LM (5-3 = 2).');
    testCase.verifyEqual(LM_event_final(idx3), 1, 'Patient 3 event indicator preserved.');
    testCase.verifyEqual(LM_X_final(idx3), 32, 'Patient 3 covariate value at risk preserved.');

    %% Validate Left-Truncation (Time Shifting) for Patient 4
    idx4 = find(LM_pat_id == 4);
    testCase.verifyEqual(numel(idx4), 2, 'Patient 4 should retain 2 intervals overlapping or post-dating the landmark.');
    
    % Interval 1 for Pt 4: Originally [0, 4]. Shifted by 3 -> [-3, 1].
    testCase.verifyEqual(LM_t_start_final(idx4(1)), 0, ...
        'Patient 4 interval 1 start left-truncated max(0, -3) to 0.');
    testCase.verifyEqual(LM_t_stop_final(idx4(1)), 1, ...
        'Patient 4 interval 1 stop shifted (4-3 = 1).');
        
    % Interval 2 for Pt 4: Originally [4, 6]. Shifted by 3 -> [1, 3].
    testCase.verifyEqual(LM_t_start_final(idx4(2)), 1, ...
        'Patient 4 interval 2 start properly shifted (4-3 = 1).');
    testCase.verifyEqual(LM_t_stop_final(idx4(2)), 3, ...
        'Patient 4 interval 2 stop properly shifted (6-3 = 3).');
end
