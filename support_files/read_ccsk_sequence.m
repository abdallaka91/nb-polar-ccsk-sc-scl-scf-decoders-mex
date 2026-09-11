function sequence = read_ccsk_sequence(filename, gf_size)
% Each text-file row contains: GF size, then the binary CCSK sequence.
fid = fopen(filename, 'r');
if fid == -1, error('Cannot open file: %s', filename); end
cleanup = onCleanup(@() fclose(fid));
rows = textscan(fid, '%f %s');
index = find(rows{1} == gf_size, 1);
if isempty(index), error('No CCSK sequence for GF%d.', gf_size); end
bits = rows{2}{index};
assert(numel(bits) == gf_size && all(bits == '0' | bits == '1'), ...
       'Invalid CCSK sequence for GF%d.', gf_size);
sequence = bits(:) - '0';
end