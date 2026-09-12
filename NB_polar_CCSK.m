% -----------------------------------------------------------------------------
% Name of file  : NB_polar_CCSK.m
% Description   : Monte Carlo simulation of Non-Binary Polar Codes (NB-PC)
%                 using CCSK modulation and frequency-domain MEX decoding.
%                 This script simulates encoding, channel transmission, and
%                 decoding over GF(q) for a given SNR and code length.
%
% Authors       : Abdallah Abdallah (<abdallah.abdallah@univ-ubs.fr>)
%                 Bertrand Le Gal (<bertrand.le-gal@univ-rennes.fr>)
%                 Camille Moniere (<camille.moniere@univ-ubs.fr>)
%                 Emmanuel Boutillon
% Organisation  : Lab-STICC, UMR 6284, Universite Bretagne Sud
% Date          : September 12, 2026
% Licence       : CeCILL-B, see LICENSE file:
%                 https://cecill.info/licences/Licence_CeCILL-B_V1-en.html
%
% Note          : This work has been funded by the French ANR project MIOT
%                 under grant Projet-ANR-24-CE93-0017
%                 Web site: https://project.inria.fr/miot/
% -----------------------------------------------------------------------------
% History:
%   Creation    : 13/10/2025 - Initial MATLAB simulator version
%   Update      : 12/09/2026 - Added MEX/DLL decoder interface
% -----------------------------------------------------------------------------
% FUNCTIONAL DESCRIPTION:
%   This script performs Monte Carlo simulations of non-binary polar codes
%   using Cyclic Code Shift Keying (CCSK) modulation. The decoder is provided
%   by the MEX/DLL files in support_files/.
%
%   The simulation includes:
%       - Reliability sequence loading from support_files/matrices/
%       - Polar encoding and CCSK modulation
%       - AWGN channel transmission
%       - Frequency-domain LLR computation
%       - MEX/DLL non-binary polar decoding
%       - Frame Error Rate (FER) estimation
%
% -----------------------------------------------------------------------------
% INPUT PARAMETERS (set by user in the Configuration section):
%   code_length        : N, original polar code length.
%   transmitted_length : Ns, length after shortening.
%   information_length : K, number of information symbols.
%   gf_size            : q, GF(q) field order.
%   snr_db             : Channel SNR in dB.
%   decoder_type       : MEX/DLL decoder selection string.
%   max_frames         : Number of simulated frames.
%   frames_per_call    : Number of frames sent to the MEX at once.
%   num_threads        : Number of decoder instances/threads used by the MEX.
%
% OUTPUT:
%   Console output with progressive FER estimation at the configured SNR.
% -----------------------------------------------------------------------------
clear
tic
%% Configuration
% MEX accepts N = 2, 4, 8, ..., 65536 (powers of two).
% MEX accepts q = 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096.
% This simulator also needs reliability data for (N,q) and a CCSK sequence for q.
code_length = 128;         % N: original code length
transmitted_length = 128;  % Ns: length after shortening
information_length = 90;   % K: number of information symbols
gf_size = 256;              % q: number of GF symbols
snr_db = -7.5;
max_frames =4e4;
frames_per_call = 200;
num_threads = 4;          % Number of decoder threads inside the MEX
random_seed = 0;
progress_interval = 100;   % Update the display every this many frames

decoder_type = 'dec4';
% Decoder choices supported by the MEX/DLL pair attached with this simulator.
% SC decoders:
%   dec1 or naive       : naive SC decoder, no fast cancellation/pruning.
%   dec3 or sc_spec_f32 : specialized SC decoder.
%   dec4 or pruned      : pruned SC decoder.
%
% SCL decoders. L is the list size, normally 2, 4, 8, 16, or 32.
% S is the split value, usually equal to L. M is the metric, allowed 1..9.
% D is debug, 0 or 1.
%   scl                 : alias for scl_naive_f32{4}
%   f_scl               : alias for f_scl_naive_f32{4}
%   sscl                : alias for scl_spec_f32{4}
%   f_sscl              : alias for f_scl_spec_f32{4}
%   pscl                : alias for scl_pruned_f32{4}
%   f_pscl              : alias for f_scl_pruned_f32{4}
%   scl_naive_f32{L},          f_scl_naive_f32{L}
%   scl_spec_f32{L},           f_scl_spec_f32{L}
%   scl_pruned_f32{L},         f_scl_pruned_f32{L}
%   scl_pruned_zc_f32{L-S-M-D}, f_scl_pruned_zc_f32{L-S-M-D}
%   scl_pruned_zc_f32_t{L-S-M}, f_scl_pruned_zc_f32_t{L-S-M}
%   scl_zc_f32{L-S-M-D},        f_scl_zc_f32{L-S-M-D}
%   scl_zc_f32_t{L-S-M},        f_scl_zc_f32_t{L-S-M}
%   scl_zc_f32_c{L-S-M},        f_scl_zc_f32_c{L-S-M}
%   scl_wave{L},                f_scl_wave{L}
%
% SCF decoders. T is the number of flip attempts/iterations.
%   scf                 : alias for scf_naive_f32{4}
%   f_scf               : alias for f_scf_naive_f32{4}
%   rscf                : alias for scf_related_f32{4}
%   f_rscf              : alias for f_scf_related_f32{4}
%   scf_naive_f32{T},          f_scf_naive_f32{T}
%   scf_related_f32{T},        f_scf_related_f32{T}
%


%% Setup
script_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(script_folder, 'support_files'));
rng(random_seed);

% Check that reliability data exist for the requested N and GF before reading.
matrix_folder = fullfile(script_folder, 'support_files', 'matrices', sprintf('GF%d', gf_size));
matrix_filename = sprintf('GF%dN%d.txt', gf_size, code_length);
matrix_path = fullfile(matrix_folder, matrix_filename);
if ~isfile(matrix_path)
    error('No reliability data available for N=%d and GF=%d. Missing file: %s', ...
          code_length, gf_size, matrix_path);
end
reliability_data = read_matrix_reliability(matrix_path, snr_db, code_length);
reliability_order = reliability_data.relab_seq_prob(:) + 1;

% Select the CCSK reference sequence for the chosen GF size
reference_sequence = read_ccsk_sequence( ...
    fullfile(script_folder, 'support_files', 'ccsk_sequences.txt'), gf_size);

% Precompute modulation and channel constants
shifted_sequences = zeros(gf_size, gf_size);
for symbol = 0:gf_size-1
    shifted_sequences(:, symbol+1) = circshift(reference_sequence, -symbol);
end
modulation_table = 1 - 2*shifted_sequences; % Bit 0 -> +1, bit 1 -> -1
reference_fft = fft(reference_sequence(:));
noise_std = sqrt(1 / 10^(snr_db/10));

% Allocate once; automatically choose N-Ns shortened positions and the frozen inputs.
% The MEX expects a most-to-least reliability permutation using indices 1:N.
[initialized, decoder_config] = nbdecode_mex('init', decoder_type, ...
    code_length, gf_size, information_length, transmitted_length, reliability_order, num_threads);
assert(initialized);
% Use these input positions for information symbols; all other inputs stay zero.
information_positions = decoder_config.information_positions;
% Send only these output positions, in the returned order.
transmitted_positions = decoder_config.transmitted_positions;

%% Simulation
frame_errors = 0;
frame_error_rate = 0;
progress_message = sprintf("SNR_dB = %.3f dB, FER = %d/%d = %.8f\n", ...
                           snr_db, 0, 0, 0);
fprintf(progress_message);

for first_frame = 1:frames_per_call:max_frames
    batch_frames = min(frames_per_call, max_frames - first_frame + 1);
    processed_frames = first_frame + batch_frames - 1;

    information_symbols = randi([0 gf_size-1], information_length, batch_frames);
    input_symbols = zeros(code_length, batch_frames); % Frozen symbols stay zero
    input_symbols(information_positions, :) = information_symbols;
    encoded_symbols = encode(input_symbols);

    transmitted_signal = zeros(gf_size, transmitted_length, batch_frames);
    for frame = 1:batch_frames
        transmitted_signal(:, :, frame) = modulation_table( ...
            :, encoded_symbols(transmitted_positions, frame) + 1);
    end

    noise = noise_std * randn(size(transmitted_signal));
    received_signal = transmitted_signal + noise;
    llrs = LLR_CCSK_FFT(received_signal, reference_fft, noise_std);
    probabilities = exp(-llrs);
    probabilities = probabilities ./ sum(probabilities, 1);

    % Input: single probabilities, q-by-Ns-by-batch_frames (not LLRs).
    % Reuses the decoder and fills shortened positions internally.
    % Output: N-by-batch_frames input symbols, including frozen positions.
    decoded_symbols = double(nbdecode_mex('decode', single(probabilities)));

    frame_errors = frame_errors + sum(any(decoded_symbols ~= input_symbols, 1));
    frame_error_rate = frame_errors / processed_frames;

    if floor(processed_frames / progress_interval) > floor((first_frame - 1) / progress_interval) ...
            || processed_frames == max_frames
        fprintf(repmat('\b', 1, length(char(progress_message))));
        progress_message = sprintf("SNR_dB = %.3f dB, FER = %d/%d = %.8f\n", ...
            snr_db, frame_errors, processed_frames, frame_error_rate);
        fprintf(progress_message);
    end
end

% Release the persistent decoder after the simulation.
nbdecode_mex('clear');
toc