function [varargout] = hNRFrequencyOffset(operationType,varargin)
%hNRFrequencyOffset Frequency offset estimation/correction
%   [FOFFSET] = hNRFrequencyOffset('COARSEFO',...) Coarse Frequency offset
%   estimation
%
%   FOFFSET = hNRFrequencyOffset('COARSEFO',CARRIER,WAVEFORM,SAMPLERATE)
%   estimates the frequency offset (FO) FOFFSET for the specified carrier
%   configuration CARRIER, waveform WAVEFORM and sample rate SAMPLERATE, by
%   means of correlation of the cyclic prefix.
%   
%   CARRIER      - Carrier configuration object, <a
%                  href="matlab:help('nrCarrierConfig')"
%                  >nrCarrierConfig</a>
%   WAVEFORM     - T-by-R matrix where T is the number of time domain
%                  samples and R is the number of receive antennas
%   SAMPLERATE   - Sample rate of the waveform
%   FOFFSET      - Output measured frequency offset in Hertz (Hz)
%                  This function can only accurately estimate frequency
%                  offsets of up to half a subcarrier spacing
%
%   FOFFSET =
%   hNRFrequencyOffset('COARSEFO',CARRIER,WAVEFORM,SAMPLERATE,TOFFSET)
%   estimates the frequency offset FOFFSET taking into account the sample
%   rate SAMPLERATE and timing offset TOFFSET.
%   The input timing offset TOFFSET is used to select the location within
%   the cyclic prefix for phase angle estimation
%   TOFFSET      - Positive scalar integer
%   Note that if TOFFSET is absent,the quality of the internal timing
%   estimate is subject to the length and signal quality of the input
%   WAVEFORM and therefore may result in inaccurate frequency offset
%   measurements
%
%   [OUT,AVERAGEFO] = hNRFrequencyOffset('FINEFO',...) 
%   Fine Frequency offset estimation
%
%   [OUT,AVERAGEFO] =
%   hNRFrequencyOffset('FINEFO',CARRIER,BWPCFG,WAVEFORM,REFGRID,...
%             IDEALGRID,CDMLENGTHS,SAMPLERATE,ACTIVESLOTS,CARRIERFREQUENCY)
%   estimates the average frequency offset AVERAGEFO by means of estimation
%   of the phase shift across channel estimates. A slot based grid of
%   known IQs REFGRID is used to calculate the channel coefficients. Note
%   that the quality of the fine FO estimate is subject to the resource
%   block  allocation length and signal quality of the input waveform
%   WAVEFORM and therefore may result in inaccurate frequency offset
%   estimates
%   CARRIER          - Carrier configuration object, <a
%                      href="matlab:help('nrCarrierConfig')"
%                      >nrCarrierConfig</a>
%   BWPCFG is a configuration object of type, <a
%                      href="matlab:help('nrWavegenBWPConfig')"
%                      >nrWavegenBWPConfig</a>
%   WAVEFORM         - T-by-R matrix where T is the number of time domain
%                      samples and R is the number of receive antennas
%   REFGRID          - Slot based grid containing known IQs of
%                      size K-by-L-by-P. K is the number of subcarriers,
%                      given by CARRIER.NSizeGrid * 12. L is the number of
%                      symbols spanning the duration of the waveform. P is
%                      the number of reference signal ports.
%   IDEALGRID        - Slot based grid containing known data IQs of
%                      size K-by-L-by-P.
%   CDMLENGTHS       - A 2-element row vector [FD TD] specifying the
%                      length of frequency domain code division
%                      multiplexing (FD-CDM) and time domain code division
%                      multiplexing (TD-CDM) despreading to perform. A
%                      value of 1 for an element indicates no CDM and a
%                      value greater than 1 indicates the length of the CDM
%   SAMPLERATE       - Sample rate of the waveform
%   ACTIVESLOTS      - A vector corresponding to the allocated slots
%                      present in the waveform WAVEFORM
%   CARRIERFREQUENCY - Carrier frequency in Hz
%   OUT              - Output waveform with the estimated fine frequency
%                      offset correction
%   AVERAGEFO        - Output FO averaged across the allocated slots
%                      present in the waveform WAVEFORM
%
%   [OUT,AVERAGEFO] = hNRFrequencyOffset('FINEFO',...,REFINDARRAY,MRB)
%   estimates the frequency offset FOFFSET for a transform precoded
%   waveform with reference signal indices REFINDARRAY and resource block
%   allocation length MRB
%   REFINDARRAY  - Array of structures, each containing 1-based linear
%                  indices addressing a K-by-L-by-P resource array. P is
%                  the number of reference signal ports and is inferred
%                  from the range of values in REFINDARRAY.
%   MRB          - Number of resource blocks allocated to the transform
%                  precoded waveform
%
%   FOFFSET =
%   hNRFrequencyOffset('INTEGERFO',CARRIER,WAVEFORM,REFGRID,SAMPLERATE,CARRIERFREQUENCY)
%   estimates the integer frequency offset FOFFSET in the time-domain
%   waveform WAVEFORM.
%   CARRIER          - Carrier configuration object, <a
%                      href="matlab:help('nrCarrierConfig')"
%                      >nrCarrierConfig</a>
%   WAVEFORM         - T-by-R grid where T is the number of time domain
%                      samples and R is the number of receive antennas
%   REFGRID          - Slot based grid containing known IQs of size
%                      K-by-N-by-P. N is the number of symbols spanning
%                      WAVEFORM.
%   SAMPLERATE       - sample rate of the waveform
%   CARRIERFREQUENCY - Carrier frequency in Hz
%   FOFFSET          - Estimated output frequency offset (Hz). It is an
%                      integer multiple of the subcarrier spacing
%                      CARRIER.SubcarrierSpacing
%
%   OUT = hNRFrequencyOffset('FOCORRECT',WAVEFORM,SAMPLERATE,FOFFSET)
%   corrects for a specified frequency offset FOFFSET in the time-domain
%   waveform by means of simple FM modulation
%   WAVEFORM     - T-by-R matrix.
%   SAMPLERATE   - Sample rate of the input waveform
%   FOFFSET      - Frequency offset correction to be applied in Hz
%   OUT          - Output waveform with the applied frequency offset
%                  correction

% Copyright 2022-2023 The MathWorks, Inc.

    switch operationType
        case 'coarseFO'
            carrier = varargin{1};
            waveform = varargin{2};
            sampleRate = varargin{3};
            if nargin == 5
                tOffset = varargin{4};
            end
            ofdmInfo = nrOFDMInfo(carrier);
            nLayers = size(waveform,2);
            overSamplingFactor = sampleRate/sum(ofdmInfo.SymbolLengths)/1000;
            waveformZeroPadded =  [waveform; zeros(round(ofdmInfo.Nfft*overSamplingFactor),nLayers)];
            if nargin == 4
                varargout{1} = frequencyOffsetCoarse(carrier,waveformZeroPadded,sampleRate);
            else
                varargout{1} = frequencyOffsetCoarse(carrier,waveformZeroPadded,sampleRate,tOffset);
            end
        case 'fineFO'
            carrier = varargin{1};
            bwpCfg = varargin{2};
            waveform = varargin{3};
            refGrid = varargin{4};
            idealGrid = varargin{5};
            cdmLengths = varargin{6};
            sampleRate = varargin{7};
            activeSlots = varargin{8};
            carrierFrequency = varargin{9};
            transformPrecoding = false;

            % Transform precoded waveform
            if nargin >10
                ind = varargin{10};
                mrb = varargin{11};
                transformPrecoding = true;
            end

            % Keep track of location where frequency offset compensated
            % 'waveform' needs to be inserted back using currentSlotCursor
            % and currentSlotLength
            currentSlotCursor = 0;
            currentSlotLength = 0;

            rxGrid = nrOFDMDemodulate(carrier,waveform,'SampleRate',sampleRate,'CarrierFrequency',carrierFrequency);

            % Work only on the relevant BWP in the waveform to simplify indexing
            bwpStart = bwpCfg.NStartBWP-carrier.NStartGrid;
            rxGridBwp = rxGrid(12*bwpStart + 1:12*(bwpStart+bwpCfg.NSizeBWP),:,:);
            if ~isempty(idealGrid)
                idealGridBwp = idealGrid(12*bwpStart + 1:12*(bwpStart+bwpCfg.NSizeBWP),:,:);
            end

            % Calculate waveform dimensions
            ofdmInfo = nrOFDMInfo(carrier);
            L = ofdmInfo.SymbolsPerSlot;            
            nSlots = floor(size(rxGrid,2)/L);
            slotwiseFineFOEnabled = mod(sampleRate,ofdmInfo.SampleRate) == 0;

            % Store waveform corrected for the estimated fine FO
            waveformFOCorrected = [];

            % fPrevSlot is used to keep track of previous estimated FO
            % p ensures phase continuity between frequency corrected outputs
            fPrevSlot = 0.0;
            p = 1;

            % Store the accumulated fine FO over nSlots
            fTotal = 0.0;
            for slotIdx = 0:nSlots-1
                currentSlot = carrier.NSlot;
                symIdx = slotIdx*L + 1:(slotIdx+1)*L;
                if ~isempty(find(currentSlot == activeSlots,1))
                    if length(find(sum(refGrid(:,symIdx,:)))) > 1
                        foffsetEstFine = frequencyOffsetFine(carrier,rxGridBwp(:,symIdx,:),refGrid(:,symIdx,:),cdmLengths);
                    else
                        if transformPrecoding
                            resourceIdx = find(currentSlot == activeSlots,1);
                            foffsetEstFine = frequencyOffsetFine(carrier,rxGridBwp(:,symIdx,:),idealGridBwp(:,symIdx,:),cdmLengths,...
                                                                                               ind(resourceIdx).ChannelIndices,mrb);
                        else
                            if ~isempty(idealGrid)
                                foffsetEstFine = frequencyOffsetFine(carrier,rxGridBwp(:,symIdx,:),idealGridBwp(:,symIdx,:),cdmLengths);
                            else
                                foffsetEstFine = 0;
                            end
                        end
                    end

                    % Correct this slot with the estimated fine FO
                    if slotwiseFineFOEnabled
                        rxSlot = nrOFDMModulate(carrier,rxGrid(:,symIdx,:),'SampleRate',sampleRate,'CarrierFrequency',carrierFrequency);
                        currentSlotLength = size(rxSlot,1);

                        % Select relevant portion of 'waveform' for FO
                        % correction
                        if foffsetEstFine
                            rxSlotFreqCorrected = frequencyOffsetCorrect(waveform(currentSlotCursor+1: currentSlotCursor+currentSlotLength,:),sampleRate,foffsetEstFine);
                        else
                            rxSlotFreqCorrected = waveform(currentSlotCursor+1: currentSlotCursor+currentSlotLength,:);
                        end

                        % Ensure phase continuity between frequency corrected
                        % outputs
                        t = size(rxSlotFreqCorrected,1)/sampleRate;
                        p = p*exp(-1i*2*pi*fPrevSlot*t);
                        rxSlotFreqCorrected = rxSlotFreqCorrected.*p;
                    end
                    fPrevSlot = foffsetEstFine;
                    fTotal = fTotal+fPrevSlot;
                    currentSlotCursor = currentSlotCursor + currentSlotLength;
                else
                    if slotwiseFineFOEnabled
                      rxSlot = nrOFDMModulate(carrier,rxGrid(:,symIdx,:),'SampleRate',sampleRate,'CarrierFrequency',carrierFrequency);
                      currentSlotLength = size(rxSlot,1);
                      rxSlotFreqCorrected = waveform(currentSlotCursor+1: currentSlotCursor+currentSlotLength,:);
                      currentSlotCursor = currentSlotCursor + currentSlotLength;
                    end
                end
                carrier.NSlot = carrier.NSlot+1;
                if slotwiseFineFOEnabled
                    waveformFOCorrected = [waveformFOCorrected; rxSlotFreqCorrected]; %#ok<*AGROW>
                end
            end
            if ~slotwiseFineFOEnabled
                waveformFOCorrected = frequencyOffsetCorrect(waveform,sampleRate,fTotal/length(activeSlots));
            end

            % Store corrected waveform and averaged fine FO
             varargout{1} = waveformFOCorrected;
             varargout{2} = fTotal/length(activeSlots);

        case 'integerFO'
            carrier = varargin{1};
            waveform = varargin{2};
            refGrid = varargin{3};
            sampleRate = varargin{4};
            carrierFrequency = varargin{5};

            % Perform time synchronization before integer FO estimation
            tmpMag = [];
            tmpLoc = [];

            % Sweep neighboring two subcarriers around the carrier frequency
            % and store a vector of correlation magnitudes and timing offset locations
            % as output by nrTimingEstimate
            for fIdx = -2:2
                tmpWaveform = frequencyOffsetCorrect(waveform,sampleRate,carrier.SubcarrierSpacing*1e3*fIdx);
                [loc,mag] = nrTimingEstimate(carrier,tmpWaveform,refGrid,'SampleRate',sampleRate,'CarrierFrequency',carrierFrequency);
                tmpMag = [ tmpMag max(mag,[],'all')];
                tmpLoc = [tmpLoc loc];
            end

            % Get the best correlation output and timing location
            [~,l] = max(tmpMag);
            waveform = waveform(1+tmpLoc(l):end,:);
            varargout{1} = frequencyIntegerCFO(carrier,waveform,refGrid,sampleRate,carrierFrequency);

        case 'FOCorrect'
            waveform = varargin{1};
            sampleRate = varargin{2};
            foffset = varargin{3};
            varargout{1} = frequencyOffsetCorrect(waveform,sampleRate,foffset);
        otherwise
            error('Unsupported operation.');
    end
end

function [foffset] = frequencyOffsetCoarse(carrier,waveform,varargin)
%    FOFFSET = frequencyOffsetCoarse('COARSE',CARRIER,WAVEFORM,SAMPLERATE)
%    estimates the frequency offset FOFFSET taking into account the
%    sample rate SAMPLERATE
%
%    FOFFSET =
%    frequencyOffsetCoarse('COARSE',CARRIER,WAVEFORM,SAMPLERATE,TOFFSET)
%    estimates the frequency offset FOFFSET taking into account the sample
%    rate SAMPLERATE and timing offset TOFFSET

    overSamplingFactor = 1;
    toffset = [];
    foffset = 0;
    ofdmInfo = nrOFDMInfo(carrier);

    if nargin >= 3
        sampleRate  = varargin{1};
        overSamplingFactor = sampleRate/1000/sum(ofdmInfo.SymbolLengths);
        if nargin == 4
            toffset = varargin{2};
        end
    end
    if overSamplingFactor ~= 1
        waveform = resample(waveform,sum(ofdmInfo.SymbolLengths)*1000,sampleRate);
    end

    % Get the number of samples per FFT, FFT duration, number of subframes,
    % number of slots, and the number of samples per slot.
    nFFT = ofdmInfo.Nfft;
    tFFT = 1/(carrier.SubcarrierSpacing*1e3);
    samplesPerSubframe = ofdmInfo.SampleRate/1000;
    samplesPerSlot = samplesPerSubframe/ofdmInfo.SlotsPerSubframe;
    L = ofdmInfo.SymbolsPerSlot;
    grid = nrOFDMDemodulate(carrier,waveform);%,'SampleRate',sampleRate);
    nSlots = floor(size(grid,2)/L);
    nSubframes = floor(nSlots/ofdmInfo.SlotsPerSubframe);

    % Ensure the waveform spans at least one slot plus FFT length
    if size(waveform,1) < (samplesPerSlot+nFFT)
        foffset = 0;
        warning('The input waveform must span at least 1 slot plus the length of the FFT');
        return;
    end

    % Derive the length of the cyclic prefixes (CP)s. The smaller CP length
    % is used for correlation with waveform. For sake of processing
    % convenience, the larger CP length is not used as it is present for a
    % relatively lesser number of symbols in a slot.
    cpLengths = ofdmInfo.CyclicPrefixLengths;
    cpLength = cpLengths(2);

    % Use peakCorr as a comparator for antenna correlation selection
    % Loop over each receive antenna
    peakCorr = -1;
    for i = 1:size(waveform,2)

        % Form two correlator inputs, the second delayed from
        % the first by nFFT.
        corrInp1 = waveform(1:end-nFFT,i);
        corrInp2 = waveform(1+nFFT:end,i);

        % Conjugate multiply the inputs and integrate over the CP length
        cpCorrRaw = corrInp1.*conj(corrInp2);
        cpcorr = conv(cpCorrRaw,ones(cpLength,1));
        cpcorr = cpcorr(cpLength:end);

        % Combine the correlation estimates, if needed, into a single subframe
        cpcorravg = cpcorr;
        if nSubframes > 0
            % Ensure cpcorravg spans a multiple of subframe length before
            % combining into a single subframe (zero pad if needed)
            sfZeroPadLen = length(cpcorr) - nSubframes*samplesPerSubframe;
            if sfZeroPadLen > 0 
               sfZeroPadLen = samplesPerSubframe-sfZeroPadLen;
            end
            cpcorravg = [cpcorr; zeros(sfZeroPadLen,1)];
            cpcorravg = sum(reshape(cpcorravg,samplesPerSubframe,length(cpcorravg)/samplesPerSubframe),2);
        end

        % Store absolute value of the averaged output for further
        % processing
        cpcorrmag = abs(cpcorravg);

        nSlotsToProcess = min(nSlots,ofdmInfo.SlotsPerSubframe);
        currentslotOffset = 0;
        cycShiftTot = [];
        sampleCount = 0;
        tmpSlotIdx = carrier.NSlot;

        % Select and process the antenna with the highest correlation peak
        % This is done to improve the estimation accuracy
        if (max(cpcorrmag) > peakCorr)

            % Update the peak value with this antenna
            peakCorr = max(cpcorrmag);

            % Loop over each slot in subframe
            for sIdx = 1:nSlotsToProcess

                % Extract the correlation samples for current slot
                slotIdxInSf = mod(tmpSlotIdx,ofdmInfo.SlotsPerSubframe);
                samplesPerSlot = sum(ofdmInfo.SymbolLengths(slotIdxInSf*ofdmInfo.SymbolsPerSlot + (1:ofdmInfo.SymbolsPerSlot)));
                cpcorrmagSlot = cpcorrmag(currentslotOffset(sIdx,1) + 1:currentslotOffset(sIdx,1) + samplesPerSlot);

                % For this slot, get the CP lengths
                symIdx = slotIdxInSf*ofdmInfo.SymbolsPerSlot + (1:ofdmInfo.SymbolsPerSlot);
                cpLengths = (ofdmInfo.CyclicPrefixLengths(symIdx));

                % Compute the number of samples through a slot where
                % each of the OFDM symbols start
                symbolStarts = cumsum(cpLengths+nFFT);
                symbolStarts0 = [0 symbolStarts(1:end-1)];

                % Extract the timing of the peak correlation relative to
                % the start of the nearest OFDM symbol in the slot. Choose
                % the closest candidate relative to any symbol start
                idx = 1:length(cpcorrmagSlot);
                tmpCycshift = idx(cpcorrmagSlot==max(cpcorrmagSlot));
                tmpCycshift = tmpCycshift(1);
                tail = samplesPerSlot-symbolStarts0(end);
                if tmpCycshift > symbolStarts0(end)
                    if (tmpCycshift-symbolStarts0(end)) > tail/2
                        tmpCycshift = mod(tmpCycshift+tail,samplesPerSlot)-tail;
                    end
                end
                candidates = -symbolStarts0+tmpCycshift;
                candidateidxs = 1:length(symbolStarts0);
                pos=candidateidxs(abs(candidates) == min(abs(candidates)));

                % Locate candidate postion for the chosen peak
                cycshift = min(candidates(pos));

                % If provided as an input, override cycshift with toffset
                if ~isempty(toffset)
                    cycshift = toffset;
                end

                % Store the cyclic shift across each slot in the
                % subframe
                cycShiftTot = [cycShiftTot cycshift];
                tmpSlotIdx = tmpSlotIdx + 1;
                sampleCount  = sampleCount + samplesPerSlot;
                currentslotOffset = [currentslotOffset; sampleCount];
            end

            % Use the median cyclic shift across each slot for better
            % accuracy
             cycShift = fix(median(cycShiftTot));

            % Form a vector of the locations of all OFDM symbols in the
            % original waveform.
            freqIdxs=[];
            prevSlotSampleCount = 0;

            for l=0:nSlots-1
                currentslotIdxInSf = mod(carrier.NSlot + l,ofdmInfo.SlotsPerSubframe);
                symIdx = currentslotIdxInSf*ofdmInfo.SymbolsPerSlot + (1:ofdmInfo.SymbolsPerSlot);
                cpLengths = (ofdmInfo.CyclicPrefixLengths(symIdx));
                symbolStarts = cumsum(cpLengths+nFFT);
                samplesInSlot = sum(ofdmInfo.SymbolLengths(currentslotIdxInSf*ofdmInfo.SymbolsPerSlot + (1:ofdmInfo.SymbolsPerSlot)));
                freqIdxs = [freqIdxs (prevSlotSampleCount + symbolStarts)];
                prevSlotSampleCount = prevSlotSampleCount + samplesInSlot;
            end
            freqIdxs(end) = [];

            estimates = [];
            cpLengths = ofdmInfo.CyclicPrefixLengths;

            % Form a vector of all samples of the correlator input,
            % centered around the middle of the cyclic prefix in each
            % symbol
            cpLength = cpLengths(2)/2;
            maxIdx = max(cycShift+freqIdxs+cpLength-1-fix(cpLength/4));
            if maxIdx > length(cpCorrRaw)
                estimates = 0;
            else
                for add=(fix(cpLength/2):fix(cpLength)-1) - fix(cpLength/4)
                    estimates = [estimates; cpCorrRaw(cycShift+freqIdxs+add)];
                end
            end

            % Average the estimates, take the angle and compute the
            % corresponding frequency offset.
            foffset = -angle(mean(estimates))/(2*pi*tFFT);
        end
    end
end

function [out] = frequencyOffsetCorrect(waveform,sampleRate,foffset)
%  OUT = frequencyOffsetCorrect(WAVEFORM,SAMPLERATE,FOFFSET)
%  Returns a waveform OUT corrected for a specified frequency offset
%  FOFFSET in the time-domain waveform WAVEFORM by means of simple FM modulation.
%  The sample rate SAMPLERATE is taken into account for FO correction.

    % Preallocate the output, taking into account the sample rate
    out=zeros(size(waveform));
    t=((0:size(waveform,1)-1)/sampleRate).';
    
    % Apply the frequency offset correction for each antenna
    for i=1:size(waveform,2)
        out(:,i) = waveform(:,i).*exp(-1i*2*pi*foffset*t);
    end
end

function [foffset] = frequencyOffsetFine(carrier,rxGrid,refGrid,cdmLengths,varargin)
%   OUT = frequencyOffsetFine(CARRIER,RXGRID,REFGRID,CDMLENGTHS)
%   estimates the average frequency offset FOFFSET by means of estimation
%   of the phase shift across channel estimates by considering these
%   inputs:
%   CARRIER      - Carrier configuration object
%   RXGRID       - Slot based grid containing demodulated IQs
%   REFGRID      - Slot based grid containing known IQs
%   CDMLENGTHS   - A 2-element row vector [FD TD] specifying the 
%                  length of FD-CDM and TD-CDM despreading to perform
%
%   OUT = frequencyOffsetFine(CARRIER,RXGRID,REFGRID,CDMLENGTHS,IND,MRB)
%   estimates the average frequency offset FOFFSET by considering these
%   inputs:
%   CARRIER      - Carrier configuration object
%   RXGRID       - Slot based grid containing demodulated IQs
%   REFGRID      - Slot based grid containing known IQs
%   CDMLENGTHS   - A 2-element row vector [FD TD] specifying the 
%                  length of FD-CDM and TD-CDM despreading to perform
%   IND          - Column vector of 1-based linear indices for a slot
%   MRB          - Number of resource blocks associated with the PUSCH
    transformPrecoding = 0;
    if nargin > 4
       transformPrecoding = 1;
       ind = varargin{1};
       mrb = varargin{2};
    end
    ofdmInfo = nrOFDMInfo(carrier);
    symbolsPerSlot = ofdmInfo.SymbolsPerSlot;
    totalNrSlots = floor(size(rxGrid,2)/symbolsPerSlot);
    L = ofdmInfo.SymbolsPerSlot;
    y = [];
    weights = 0;

    % Store number of receive antenna
    nrAnts = size(rxGrid,3);

    % for each slot
    for ii = 0:(totalNrSlots-1)
        % Pick the indices that are relevant for this slot
        symIdx = ii*L+1:(ii+1)*L;
        rxGridSlot = rxGrid(:,symIdx,:);
        refGridSlot = refGrid(:,symIdx,:);
        numSCs = size(rxGridSlot,1);
        if ~transformPrecoding
            ind = find(refGridSlot);
        end

        if transformPrecoding
            tmp = rxGridSlot(ind);
            tmp = nrTransformDeprecode(tmp,mrb);
            rxGridSlot(ind) = tmp;
            tmp = refGridSlot(ind);
            tmp = nrTransformDeprecode(tmp,mrb);
            refGridSlot(ind) = tmp;
        end

        knownGrid = zeros(numSCs,L,nrAnts);
        knownGrid(ind) = 1;

        % Remove first non-zero column because of differential product
        knownGridSummed = sum(knownGrid(:,:,1));
        firstColumn = find(knownGridSummed,1);

        % extract the current slot of the channel estimate
        % estimate the LS channel coefficients
        hestSlot = nrChannelEstimate(rxGridSlot,refGridSlot,'CDMLengths',cdmLengths,'AveragingWindow',[11 1]);
        hestSlot = sum(hestSlot,[3 4]);

        % Remove invalid values
        ind = isinf(hestSlot);
        hestSlot(ind) = 0;
        ind = isnan(hestSlot);
        hestSlot(ind) = 0;

        % Add differential product to remove drift within one OFDM
        % symbol
        fc = firstColumn;
        hestSlot(:,fc+1:L,:) = hestSlot(:,fc+1:L,:) .* conj(hestSlot(:,fc,:));

        y = [y; sum(hestSlot)];

        % Store number of symbols per column in 'weights'
        weights = weights + knownGridSummed;

        y(:,firstColumn) = 0;
    end
    if (numel(find(weights)) < 2)
        % At least two points are needed to infer frequency offset
        foffset = 0;
    else

        % Take the angle across multiple slots and antennas
        y = mean(angle(mean(y,1)),3);

        % Linear regression with each column weighted by number of symbols
        x = 1:numel(y);
        S = sum(weights);
        X2 = sum(weights.*x.^2)/S;
        XY = sum(weights.*x.*y)/S;
        Y = sum(weights.*y)/S;
        X = sum(weights.*x)/S;
        slope = (XY-X*Y)/(X2-X^2);  % slope of the best linear fit
        
        % compute OFDM symbol length 'N' and cyclic prefix length 'Ng'
        N = double(ofdmInfo.Nfft);
        Ng = double(ofdmInfo.CyclicPrefixLengths(2));

        % compute frequency offset estimate based on the phase shift
        % transition/slope across the symbols give by 'delta_f_tilde'
        delta_f_tilde = 1/(2*pi*(1+Ng/N))*slope;
        foffset = (delta_f_tilde*(carrier.SubcarrierSpacing*1e3));
    end
end

function [integerCfoHz] = frequencyIntegerCFO(carrier,waveform,refGrid,sampleRate,carrierFrequency)
%   OUT = frequencyIntegerCFO(CARRIER,WAVEFORM,REFGRID,SAMPLERATE,CARRIERFREQUENCY)
%   estimates the integer frequency offset FOFFSET by considering these
%   inputs:
%   CARRIER          - Carrier configuration object, <a
%                      href="matlab:help('nrCarrierConfig')"
%                      >nrCarrierConfig</a>
%   WAVEFORM         - Time domain samples
%   REFGRID          - Grid spanning one or more slots of known IQs with size
%                      K-by-N-by-P.
%   SAMPLERATE       - Sample rate of the input waveform
%   CARRIERFREQUENCY - Carrier frequency in Hz
%   OUT              - Output waveform with the applied frequency offset
%                      correction
%   The approach used follows Equation 49 of Ch. Nanda Kishore.; V.
%   Umapathi Reddy., "A Frame Synchronization and Frequency Offset
%   Estimation Algorithm for OFDM System and its Analysis," EURASIP Journal
%   on Wireless Communications and Networking, Volume 2006, Article ID
%   57018, Pages 1–16

    % Get OFDM related parameters
    ofdmInfo = nrOFDMInfo(carrier);
    L = ofdmInfo.SymbolsPerSlot;

    % Calculate the number of slots
    rxGrid = nrOFDMDemodulate(carrier,waveform,'SampleRate',sampleRate,'CarrierFrequency',carrierFrequency);
    nSlots = floor(size(rxGrid,2)/L);

    % Remodulate the waveform with default sample rate
    waveform = nrOFDMModulate(carrier,rxGrid);

    % OFDM modulate the reference grid
    refWaveform = nrOFDMModulate(carrier,refGrid);

    % Zero-pad to ensure waveform and rxWaveform are of same lengths
    refLen = size(refWaveform,1);
    rxLen = size(waveform,1);
    refWaveform = conj(refWaveform);
    if refLen < rxLen
        refWaveform = [refWaveform; zeros((rxLen - refLen), size(refWaveform,2))];
    end
    if rxLen < refLen
        refWaveform(rxLen+1:end,:) = [];
    end

    % Multiply received waveform with the complex conjugate of the time-domain waveform representing refGrid
    % The peak location of the above result is the integer FO
    [T,R] = size(waveform);
    P = size(refWaveform,2);
    mag = zeros([T R P],'like',waveform);

    % For each receive antenna
    for r = 1:R

        % For each port
        for p = 1:P

            % Multiply the given antenna of the received signal with the
            % given port of the reference signal
            x = waveform(:,r).*refWaveform(:,p);
            x = fft(x);
            mag(:,r,p) = abs(x);
        end
    end

    % Sum the magnitudes of the ports
    mag = sum(mag,3);
    [~,loc] = max(sum(mag,2));

    % Handle 'wrap-around' case, for negative integer FO values
    if loc > rxLen/2
        loc = -1*abs(loc-rxLen);
    end

    scIdx = fix(loc/(nSlots*ofdmInfo.SymbolsPerSlot));
    integerCfoHz = scIdx*carrier.SubcarrierSpacing*1e3;
end