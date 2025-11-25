function [spectra, fb] = Online_psd(data, fs, fmax, nff)
% online psd

bsize = size(data, 2);
nsc = floor(bsize/4.5);
nov = floor(nsc/2);
nff = max(nff, 2^nextpow2(nsc));

for m = 1:size(data,1)
    [spectra(:, m), fb] = pwelch(data(m, :), hamming(nsc), nov, nff, fs);
end

fuse = logical((fb > 0) .* (fb < fmax));
spectra = spectra(fuse, :);
fb = fb(fuse);