function [doppFilCoeff] = get_dopp_filter(fdop, fschan)

fsdop = 1e3;
delf = 1;
numHalfLobes = 20;
f = [-fsdop:delf:fsdop];
doppSpect = 1/(pi*fdop)*1./(1-(f/fdop).^2);
doppSpect(find(f>=fdop | f<=-fdop)) = 0;
doppSpect = sqrt(doppSpect);
doppFilCoeff = real(fftshift(ifft(fftshift(doppSpect))));

% capture the main lobe portion
mainLobeSampleNum = round(fsdop/delf);
[~, nmin]  = min(abs(round(-numHalfLobes*fsdop/fdop)-f));
[~, nmax] = min(abs(round(numHalfLobes*fsdop/fdop)-f));
doppFilCoeff   = doppFilCoeff(nmin:nmax);

[doppFilCoeff] = resample(doppFilCoeff, fschan, fsdop);
doppFilCoeff(end-round(fschan/fsdop)+1:end) = []; % to retain symmetry of filter coefficients

% detect 1st/2nd zero crossing
signVec = sign(doppFilCoeff);
signSignVec = signVec(2:end) - signVec(1:end-1);
zeroCrossingVec = find(signSignVec ~= 0);
firstZeroCrossingIndex = zeroCrossingVec(1);
doppFilCoeff(1:firstZeroCrossingIndex) = [];
doppFilCoeff(end-firstZeroCrossingIndex+1:end) = [];
doppFilCoeff = doppFilCoeff/sqrt(sum(doppFilCoeff.^2));

