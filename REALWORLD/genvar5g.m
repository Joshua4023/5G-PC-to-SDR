function [carrier, pdsch, NHARQProcesses, rvSeq, codeRate, encodeDLSCH, decodeDLSCH] = genvar5g

%% Carrier config
carrier = nrCarrierConfig;
carrier.NCellID = 1;
carrier.SubcarrierSpacing = 15;
carrier.CyclicPrefix = 'normal';
carrier.NSizeGrid;
carrier.NSlot;
carrier.NFrame;
carrier.IntraCellGuardBands;
carrier
%Read-Only variables
%SymbolsPerSlot: Number of OFDM symbols per slot bases on CyclePrefix
%SlotsPerSubFrame: Number of Slots per 1ms frame, based of the subcarrierspacing
%SlotsPerFrame: Number of Slots per 10ms frame, based of the subcarrierspacing

%% Physical Downlink Shared Channel
pdsch = nrPDSCHConfig;
pdsch.NSizeBWP;
pdsch.NStartBWP;
pdsch.ReservedPRB;
pdsch.ReservedRE;
pdsch.Modulation = "16QAM";
pdsch.NumLayers = 1;
pdsch.MappingType;
pdsch.SymbolAllocation;
pdsch.PRBSet = 0:carrier.NSizeGrid-1;     % Full band allocation
pdsch.PRBSetType;
pdsch.VRBToPRBInterleaving;
pdsch.VRBBundleSize;
pdsch.NID;
pdsch.RNTI;

pdsch

%% PDSCH-DeModulation Refrence Signal
pdsch.DMRS.DMRSConfigurationType = 1;
pdsch.DMRS.DMRSReferencePoint;
pdsch.DMRS.DMRSTypeAPosition; 
pdsch.DMRS.DMRSAdditionalPosition;
pdsch.DMRS.DMRSLength = 2;
pdsch.DMRS.CustomSymbolSet; 
pdsch.DMRS.DMRSPortSet;
pdsch.DMRS.NIDNSCID; 
pdsch.DMRS.NSCID; 
pdsch.DMRS.NumCDMGroupsWithoutData;
pdsch.DMRS.DMRSDownlinkR16;
pdsch.DMRS.CDMGroups;
pdsch.DMRS.DeltaShifts;
pdsch.DMRS.FrequencyWeights;
pdsch.DMRS.TimeWeights;
pdsch.DMRS.DMRSSubcarrierLocations;
pdsch.DMRS.CDMLengths;

pdsch.DMRS

%% DL-SCH Configuration
NHARQProcesses = 1     % Number of parallel HARQ processes
rvSeq = 0              % 0 = disabled HARQ retransmission incase of an error

%% Coding rate
if pdsch.NumCodewords == 1 
    codeRate = 490/1024; %For up to 4 layers 
else
    codeRate = [490 490]./1024; %For above 4 layers
end

codeRate

%% Create DL-SCH encoder object
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = codeRate;

encodeDLSCH

%% Create DLSCH decoder object
decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = codeRate;
decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
decodeDLSCH.MaximumLDPCIterationCount = 6;

decodeDLSCH

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




