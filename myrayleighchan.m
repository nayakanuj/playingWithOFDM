function [multipathChanOp, sigState, chanStateUpdated, idealChannel, doppFilResampled, awgnDoppLastSamp] = myrayleighchan(txSig, chStruct, pathGains, pathDelays, fs, fsChan, sigState, chanState, awgnDoppLastSamp)

npts = size(txSig, 2);
doppFilLen = length(chStruct.doppFilCoeff);
numPaths = length(pathDelays);
pathDelaysSamp = round(pathDelays*fs);
maxPathDelaySamp = max(pathDelaysSamp);
nptsChan = round(npts*fsChan/fs);

if nptsChan == 0
    nptsChan = 1;
end

% filter white gaussian noise for each path here - doppler
% generate white gaussian noise
if isequal(awgnDoppLastSamp, []) || ~isequal(size(awgnDoppLastSamp, 1), numPaths)
    awgnDoppLastSamp = (randn(numPaths, 1) + 1j*randn(numPaths, 1))*(1/sqrt(2));
end

wnoiseMat = [awgnDoppLastSamp (1/sqrt(2)*(randn(numPaths, nptsChan)+1j*randn(numPaths, nptsChan)))];

awgnDoppLastSamp = wnoiseMat(:, end);

% filter it using doppler filter
if isequal(chanState, []) || ~isequal(size(chanState, 1), numPaths) || ~isequal(size(chanState, 2), doppFilLen-1)
    if doppFilLen >= npts
        awgnDummy = [(1/sqrt(2))*(randn(numPaths, doppFilLen)+1j*randn(numPaths, doppFilLen)) wnoiseMat];
        % initialization
        [~, chanState] = filter(chStruct.doppFilCoeff.', 1, awgnDummy.');
        chanState = chanState.';
    else
        disp('doppler is too high'); % consider revision
    end    
end

for pathInd = 1:numPaths
    % consider the first sample of the gaussian noise (input to the doppler filter) of the next block for interpolation
    [realDoppFilOp, ~] = filter(chStruct.doppFilCoeff.', 1, real(wnoiseMat(pathInd, :).'), real(chanState(pathInd, :)));
    [imagDoppFilOp, ~] = filter(chStruct.doppFilCoeff.', 1, imag(wnoiseMat(pathInd, :).'), imag(chanState(pathInd, :)));
    
    % determine the doppler filter state discarding the first sample of the next block (of gaussian noise)
    [~, realChanStateUpdated] = filter(chStruct.doppFilCoeff.', 1, real(wnoiseMat(pathInd, 1:end-1).'), real(chanState(pathInd, :)));
    [~, imagChanStateUpdated] = filter(chStruct.doppFilCoeff.', 1, imag(wnoiseMat(pathInd, 1:end-1).'), imag(chanState(pathInd, :)));
    
    doppFilOp(:, pathInd) = (realDoppFilOp + 1j*imagDoppFilOp);
    chanStateUpdated(:, pathInd) = realChanStateUpdated + 1j*imagChanStateUpdated;
    doppFilOpUpsampledLen = round(fs/fsChan)*(nptsChan+1);
    doppFilResampled(:, pathInd) = interp1(linspace(1, npts+1, nptsChan+1), doppFilOp(:, pathInd), [1:npts+1]);
end

chanStateUpdated = chanStateUpdated.';

% apply gain to each path
doppFilResampled(npts+1:end, :) = []; % this introduces a minor error of 0 or 1 sample
rayleighChanProfile = (doppFilResampled.*(ones(npts, 1) * (10.^(pathGains/20)))).';
if isequal(sigState, []) || ~isequal(size(sigState, 1), numPaths) || ~isequal(size(sigState, 2), max(pathDelaysSamp))
    sigState = zeros(numPaths, maxPathDelaySamp);
end

txSigForChannel = zeros(numPaths, npts+maxPathDelaySamp);
% apply delay to the data
for pathInd = 1:numPaths
    txSigForChannel(pathInd, pathDelaysSamp(pathInd)+[1:npts]) = txSig;
    txSigForChannel(pathInd, 1:pathDelaysSamp(pathInd)) = sigState(pathInd, 1:pathDelaysSamp(pathInd));
end

% retain channel state
sigState = txSigForChannel(:, npts+1:end);

% multiply with the data
multipathChanOpAllPaths = txSigForChannel(:, 1:npts).*(rayleighChanProfile);

% add the multipath components
multipathChanOp = sum(multipathChanOpAllPaths, 1);

% ideal channel
idealChannel = zeros(1, npts);
for indPath = 1:numPaths
    idealChannel(pathDelaysSamp(indPath)+1) = idealChannel(pathDelaysSamp(indPath)+1) + mean(rayleighChanProfile(indPath,:).');
end

 
