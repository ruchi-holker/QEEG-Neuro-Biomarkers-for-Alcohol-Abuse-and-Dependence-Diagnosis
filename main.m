% Date: 11/08/2021

%%
close all
clear
clc



%% RCSP with fetaures 
%% RCSP
fprintf('CSP processing started.\n')
W = cell(length(beta), nbands);
f_train = cell(2, length(beta));
for p = 1:length(beta)
    sel_ch = 3;  % Select columns 
    start = 1;
    feat_vec = zeros((25*6+12)*nbands, size(EEG, 2));
    for bank = 1: nbands
        X = reshape(cell2mat(EEG(bank, :)), [size(data, 1), size(data, 2), size(EEG, 2)]);
        W{p, bank} = RegCsp(X, labels, beta(p), gamma(p)); % No regularization
        vec2 = zeros(25*6, size(EEG, 2));
        vec3 = zeros(12, size(EEG, 2));
        for i = 1:size(X, 3)
            Prj_X = X(:,:,i) * W{p, bank}([1:sel_ch, 64-sel_ch+1:64], :).'; % Projection
            data_st(i).eeg_data=Prj_X.';  data_st(i).Fs = fs; data_st(i).ch_labels =[];
        end
        for i = 1:size(X, 3)
            %%% METHOD-2 %%%% (25 features per band)
%             vec2(:, i) = struct2array(generate_all_features(data_st(i), [], FEATURE_SET));
            [~, a] = generate_all_features(data_st(i), [], FEATURE_SET, true);
            vec2(:, i) = cell2mat(a);
            %%% METHOD-3 %%%% (12 features per band)
            vec3(:, i) = method_3(data_st(i).eeg_data);
        end
        vec2(isnan(vec2))= 0;  % If NaN then replace by 0
        stop = start + (25*6+12) - 1;
        feat_vec(start:stop, :) = [vec2; vec3];
        start = stop + 1;
    end
    csvwrite(['tr_park_comb_', num2str(beta(p)),'_', num2str(gamma(p)),'.csv'], [feat_vec.', labels.'])
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RCSP
fprintf('CSP processing started.\n')
f_test = cell(2, length(beta));
for p = 1:length(beta)
    sel_ch = 3;  % Select columns 
    feat_vec = zeros((25*6+12)*nbands, size(EEG, 2));
    start = 1;
    for bank = 1: nbands
        X = reshape(cell2mat(EEG(bank, :)), [size(data, 1), size(data, 2), size(EEG, 2)]);

        vec2 = zeros(25*6, size(EEG, 2));
        vec3 = zeros(12, size(EEG, 2));
        for i = 1:size(X, 3)
            Prj_X = X(:,:,i) * W{p, bank}([1:sel_ch, 64-sel_ch+1:64], :).'; % Projection
            data_st(i).eeg_data=Prj_X.';  data_st(i).Fs = fs; data_st(i).ch_labels =[];
        end
        parfor i = 1:size(X, 3)
            %%% METHOD-2 %%%% (25 features per band)
%             vec2(:, i) = struct2array(generate_all_features(data_st(i), [], FEATURE_SET));
            [~, a] = generate_all_features(data_st(i), [], FEATURE_SET, true);
            vec2(:, i) = cell2mat(a);
            %%% METHOD-3 %%%%
            vec3(:, i) = method_3(data_st(i).eeg_data);
        end
        vec2(isnan(vec2))= 0;  % If NaN then replace by 0
        stop = start + (25*6+12) - 1;
        feat_vec(start:stop, :) = [vec2; vec3];
        start = stop + 1;
    end
    csvwrite(['ts_park_comb_', num2str(beta(p)),'_', num2str(gamma(p)),'.csv'], [feat_vec.', labels.'])
end 
