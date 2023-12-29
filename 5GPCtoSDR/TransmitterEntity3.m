function [longwave, waveformInfo, trBlk, trBlkSizes] = TransmitterEntity3(harqEntity, txDataBits, nTxAnts, chan, channel, SNRdB, nRxAnts)
    
    % Wave varaibles
    [carrier, pdsch, ~, ~, codeRate, encodeDLSCH, decodeDLSCH] = genvar5g;

    % Generate PDSCH indices info, which is needed to calculate the transport
    % block size
    [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);

    % Calculate transport block sizes
    Xoh_PDSCH = 0;
    trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);
    
    estChannelGrid = getInitialChannelEstimate(channel,carrier);
    newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChannelGrid);
    
    %Data stream segmentation variable
    slstr = 1;
    slend = trBlkSizes;
    longwave = 0.00001+0.00001i;

    for nSlot = 0:10-1
        carrier.NSlot = nSlot;
    
        %data stream segmentation
        data = txDataBits(slstr:slend,1);
        slstr = slstr + trBlkSizes;
        slend = slend + trBlkSizes;

        % Get new transport blocks and flush decoder soft buffer, as required
        for cwIdx = 1:pdsch.NumCodewords
            if harqEntity.NewData(cwIdx)
                % Create and store a new transport block for transmission
                trBlk = data;
                setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);
                if harqEntity.SequenceTimeout(cwIdx)
                    resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
                end
            end
        end

        codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers,pdschInfo.G,harqEntity.RedundancyVersion);

        pdschSymbols = nrPDSCH(carrier,pdsch,codedTrBlock);

        precodingWeights = newPrecodingWeight;

        pdschSymbolsPrecoded = pdschSymbols*precodingWeights;

        dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
        dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);

        pdschGrid = nrResourceGrid(carrier,nTxAnts);
    
        [~,pdschAntIndices] = nrExtractResources(pdschIndices,pdschGrid);
        pdschGrid(pdschAntIndices) = pdschSymbolsPrecoded;

        % PDSCH DM-RS precoding and mapping
        for p = 1:size(dmrsSymbols,2)
            [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),pdschGrid);
            pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + dmrsSymbols(:,p)*precodingWeights(p,:);
        end

        %% OFDM Modulation

        [txWaveform,waveformInfo] = nrOFDMModulate(carrier,pdschGrid);
        
        %transport channel selections
        if strcmpi(chan,"OverTheAir")
            longwave = cat(1, longwave,txWaveform);
        elseif strcmpi(chan,"GaussianNoise")
            [txWaveform2,pathGains,sampleTimes] = channel(txWaveform);
            noise = generateAWGN(SNRdB,nRxAnts,waveformInfo.Nfft,size(txWaveform2));
            txWaveform2 = txWaveform2 + noise;
            longwave = cat(1, longwave,txWaveform2);
        else
            longwave = cat(1, longwave,txWaveform);
        end
    end 
    %removing empty start bit
    longwave(1) = [];
end

%% Local Functions
function noise = generateAWGN(SNRdB,nRxAnts,Nfft,sizeRxWaveform)
% Generate AWGN for a given value of SNR in dB (SNRDB), which is the
% receiver SNR per RE and antenna, assuming the channel does
% not affect the power of the signal. NRXANTS is the number of receive
% antennas. NFFT is the FFT size used in OFDM demodulation. SIZERXWAVEFORM
% is the size of the receive waveform used to calculate the size of the
% noise matrix.

    % Normalize noise power by the IFFT size used in OFDM modulation, as
    % the OFDM modulator applies this normalization to the transmitted
    % waveform. Also normalize by the number of receive antennas, as the
    % channel model applies this normalization to the received waveform by
    % default. The SNR is defined per RE for each receive antenna (TS
    % 38.101-4).
    SNR = 10^(SNRdB/10); % Calculate linear noise gain
    N0 = 1/sqrt(2.0*nRxAnts*double(Nfft)*SNR);
    noise = N0*complex(randn(sizeRxWaveform),randn(sizeRxWaveform));
end

function wtx = getPrecodingMatrix(PRBSet,NLayers,hestGrid)
% Calculate precoding matrix given an allocation and a channel estimate
    
    % Allocated subcarrier indices
    allocSc = (1:12)' + 12*PRBSet(:).';
    allocSc = allocSc(:);
    
    % Average channel estimate
    [~,~,R,P] = size(hestGrid);
    estAllocGrid = hestGrid(allocSc,:,:,:);
    Hest = permute(mean(reshape(estAllocGrid,[],R,P)),[2 3 1]);
    
    % SVD decomposition
    [~,~,V] = svd(Hest);
    
    wtx = V(:,1:NLayers).';
    wtx = wtx/sqrt(NLayers); % Normalize by NLayers
end

function estChannelGrid = getInitialChannelEstimate(channel,carrier)
% Obtain an initial channel estimate for calculating the precoding matrix.
% This function assumes a perfect channel estimate

    % Clone of the channel
    chClone = channel.clone();
    chClone.release();

    % No filtering needed to get channel path gains
    chClone.ChannelFiltering = false;    
    
    % Get channel path gains
    [pathGains,sampleTimes] = chClone();
    
    % Perfect timing synchronization
    pathFilters = getPathFilters(chClone);
    offset = nrPerfectTimingEstimate(pathGains,pathFilters);
    
    % Perfect channel estimate
    estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
end

function refPoints = getConstellationRefPoints(mod)
% Calculate the reference constellation points for a given modulation
% scheme.
    switch mod
        case "QPSK"
            nPts = 4;
        case "16QAM"
            nPts = 16;
        case "64QAM"
            nPts = 64;
        case "256QAM"
            nPts = 256;            
    end
    binaryValues = int2bit(0:nPts-1,log2(nPts));
    refPoints = nrSymbolModulate(binaryValues(:),mod);
end

function estChannelGrid = precodeChannelEstimate(estChannelGrid,W)
% Apply precoding matrix W to the last dimension of the channel estimate.

    % Linearize 4-D matrix and reshape after multiplication
    K = size(estChannelGrid,1);
    L = size(estChannelGrid,2);
    R = size(estChannelGrid,3);
    estChannelGrid = reshape(estChannelGrid,K*L*R,[]);
    estChannelGrid = estChannelGrid*W;
    estChannelGrid = reshape(estChannelGrid,K,L,R,[]);

end
    






