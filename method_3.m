function feat_vec = method_3(EEG)
%{
% Approximate Entropy
m = 3;
ApEn = zeros(size(EEG, 1), 1);
for ch = 1:size(EEG, 1)
    eeg_band_data = real(ifft(EEG(ch, :)));

    Xm = buffer(eeg_band_data, 3, 2);
    Xm = Xm(:, 1:end-m+1);

    r = 0.2*std(Xm, [], 1);
    Cm = zeros(1, size(Xm, 2));
    for i= 1:size(Xm, 2)
        d = max(bsxfun(@minus, Xm, Xm(:, 1)), [], 1);
        Cm(i) = sum((r(1) - d)>0)/(size(Xm, 2));
    end
    phi_m = mean(log(Cm));

    Xm1 = buffer(eeg_band_data, 3+1, 2);
    Xm1 = Xm1(:, 1:end-m+1);

    r1 = 0.2*std(Xm1, [], 1);
    Cm1 = zeros(1, size(Xm1, 2));
    for i= 1:size(Xm1, 2)
        d = max(bsxfun(@minus, Xm1, Xm1(:, 1)), [], 1);
        Cm1(i) = sum((r1(1) - d)>0)/(size(Xm1, 2));
    end
    phi_m1 = mean(log(Cm1));

    ApEn(ch) = phi_m - phi_m1;
end
%}
% Permutation entropy
m = 3;
Pi = perms(1:m);  % Symbols
PeEn = zeros(size(EEG, 1), 1);
for ch = 1:size(EEG, 1)
    eeg_band_data = real(ifft(EEG(ch, :)));

    Xm = buffer(eeg_band_data, 3, 2);

    Xm = Xm(:, 1:end-m+1);

    [~, I] = sort(Xm, 1);
    p = zeros(1, size(Pi, 1));
    for j=1:size(Pi, 1)
        p(j) = sum(ismember(I.', Pi(1, :), 'rows'))/factorial(3)/size(Xm, 2);
    end
    PeEn(ch) = -p*log(p).';
end    

%% Feature vector
feat_vec = [PeEn(:)];

