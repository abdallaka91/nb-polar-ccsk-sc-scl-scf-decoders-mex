function L = LLR_CCSK_FFT(y, reference_fft, sigm)
    % y: q × Ns × number_of_frames
    % reference_fft: reference CCSK sequence FFT of length q
    % sigm: scalar noise standard deviation
    L = real(ifft(reference_fft .* conj(fft(y, [], 1)), [], 1));
    L = L * (2 / sigm^2);
end