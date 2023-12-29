function [decbits, statusReport, pdschEq] = RecieverEntity(carrier, pdsch, decodeDLSCH, harqEntity, rxWaveform, offset, fileID2, newPrecodingWeight, trBlkSizes)

    precodingWeights = newPrecodingWeight;

    dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
    dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);

    [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);


    %% Timing Synchronization

    [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols);
    offset = hSkipWeakTimingOffset(offset,t,mag);
    
    rxWaveform = rxWaveform(1+offset:end,:);

    %% OFDM Demodulation

    rxGrid = nrOFDMDemodulate(carrier,rxWaveform);

    %% Channel Estimation

    % Perform practical channel estimation between layers and receive
    % antennas.
    [estChGridLayers,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsSymbols,'CDMLengths',pdsch.DMRS.CDMLengths);

    % Remove precoding from estChannelGrid before precoding
    % matrix calculation
    estChGridAnts = precodeChannelEstimate(estChGridLayers,conj(precodingWeights));

    % Get precoding matrix for next slot
    newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
    

    mesh(abs(estChGridLayers(:,:,1,1)));
    title('Channel Estimate');
    xlabel('OFDM Symbol');
    ylabel("Subcarrier");
    zlabel("Magnitude");
    
    %% Equalization

    [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChGridLayers);
    [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);
    
    %% PDSCH Decoding

    [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);

    % Scale LLRs by CSI
    
    csi = nrLayerDemap(csi);                                    % CSI layer demapping

    for cwIdx = 1:pdsch.NumCodewords
        Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); % Bits per symbol
        csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);                 % Expand by each bit per symbol
        dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % Scale    
    end

    %% DL-SCH Decoding

    
    decodeDLSCH.TransportBlockLength = trBlkSizes;

    [decbits,blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NumLayers, ...
    harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
    %decbits;

    fwrite(fileID2,decbits, 'logical');

    

    %% HARQ Process Update

    statusReport = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G); 

    

end



%% Local Functions
    
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
