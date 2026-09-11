function data = read_matrix_reliability(filename, target_snr, N)
fid = fopen(filename, 'r');
if fid == -1, error('Cannot open file: %s', filename); end
cleanup = onCleanup(@() fclose(fid));
best_distance = Inf;
data = [];
while true
    tag = fscanf(fid, '%s', 1);
    if isempty(tag), break; end
    if ~strcmp(tag, 'SNR'), error('Expected SNR in %s', filename); end
    snr = fscanf(fid, '%f', 1);
    sequence = fscanf(fid, '%d', N);
    if isempty(snr) || numel(sequence) ~= N
        error('Incomplete SNR block in %s', filename);
    end
    if abs(snr - target_snr) < best_distance
        best_distance = abs(snr - target_snr);
        data.relab_seq_prob = (sequence');
        data.snr = snr;
    end
end
if isempty(data), error('No SNR blocks in %s', filename); end
end