function [eqSlotGrid,refSlotGrid,emissions,eqSym,refSym] = hDecodeSlots(carrier,rxGrids,hest,...
    cfgArray,pxschObj,varargin)
    %hDecodeSlots Create a slot based grid of equalized symbols and
    %   reference symbols by decoding an array of demodulated IQs
    %
    %   [EQSLOTGRID,REFSLOTGRID,~,EQSYM,REFSYM] =
    %   hDecodeSlots(CARRIER,RXGRIDS,HEST,CFGARRAY,PXSCHOBJ,...)
    %   Create a grid of equalized symbols EQSLOTGRID and ideal symbols.
    %   REFSLOTGRID, by decoding an array of demodulated IQs RXGRIDS
    %   EQSLOTGRID and REFSLOTGRID are 4-dimensional matrices of size
    %   K-by-N-by-L-by-E. K is the number of subcarriers, given by
    %   CARRIER.NSizeGrid * 12. N is the number of symbols spanning the
    %   duration of the waveform. L is the number of transmitted layers. E
    %   is the number of EVM window positions.
    %   EQSYM is an output cell array of IQ constellations for low and
    %   high EVM window locations . It is of dimensions E-by-1.
    %   Each array element is of dimensions numIqSamples-by-L
    %   REFSYM is an output cell array of reference IQ constellations
    %   for low and high EVM window locations.It is of dimensions E-by-1.
    %   Each array element is of dimensions numIqSamples-by-L.
    %   CARRIER is a configuration object of type, <a
    %                        href="matlab:help('nrCarrierConfig')"
    %                        >nrCarrierConfig</a>
    %   RXGRIDS is a cell array of demodulated IQs of length 1 or 2
    %   dependent on the mode of channel estimation, non-3GPP or 3GPP
    %   respectively. Each element of RXGRIDS must be an array of size
    %   K-by-N-by-R.
    %   HEST is a cell array of channel coefficients of length 1 or 2
    %   dependent on the mode of channel estimation, non-3GPP or 3GPP
    %   respectively. Each element of HEST must be of size K-by-N-by-R-by-P
    %   array. R is the number of receive antennas and P is the number of
    %   reference signal ports.
    %   CFGARRAY is an array of structures. Each structure contains the
    %   configuration and reference resource information. The contents of
    %   this structure are used for equalization, decoding, and building
    %   reference IQs. The array itself spans the length of the
    %   concatenated configs.
    %   PXSCHOBJ is a configuration object of type 'nrPDSCHConfig' or
    %   'nrPUSCHConfig'. The absence of these configuration types will
    %   instruct this function to decode PDCCH IQs.
    %   
    %   [...] =
    %   hDecodeSlots(...,SKIPCONTROLSYMBOLS,PDSCHENCODINGON,WAVEFORMPERIOD)
    %   Create a grid of equalized symbols and reference symbols by
    %   decoding an array of demodulated PDSCH IQs .
    %   SKIPCONTROLSYMBOLS when set to true, skips symbols occupied by
    %   PDCCH when measuring PDSCH EVM.
    %   PDSCHENCODINGON when set to true enables derate matching and LDPC
    %   decoding.
    %   WAVEFORMPERIOD is in duration of slots and is relevant only to
    %   PDSCH and PDCCH. It is used to index PDSCH or PDCCH configuration
    %   resources taking into account the waveform period.
    %
    %   [EQSLOTGRID,REFSLOTGRID,EMISSONS,EQSYM,REFSYM] = hDecodeSlots(...,
    %                      BWP,EMISSONS,RXGRIDINBANDEMISSIONS)
    %   Create a grid of equalized symbols and reference symbols by
    %   decoding an array of demodulated PUSCH IQs.
    %   BWP is a configuration object of type, <a
    %                        href="matlab:help('nrWavegenBWPConfig')"
    %                        >nrWavegenBWPConfig</a>
    %   EMISSIONS structure containing in-band emissions statistics.
    %   as defined in TS 38.101-1 (FR1) / TS 38.101-2 (FR2), Annex F.3.
    %   RXGRIDINBANDEMISSIONS is an array of demodulated IQs used for the
    %   purpose of measuring in-band emissions statistics. It is of
    %   dimensions K-by-N-by-R, where K is number of subcarriers of size
    %   carrier.NSizeGrid*12, N is the number of symbols spanning nSlots
    %   and R is the number of recieve antennas.

    % Copyright 2022-2023 The MathWorks, Inc.

    decodePdsch = false;
    decodePusch = false;
    decodePdcch = false;
    emissions = [];

    if iscell(pxschObj) && isa(pxschObj{1},'nrPDSCHConfig')
        decodePdsch = true;
        skipCtrlSymbols = varargin{1};
        pdschEncodingOn = varargin{2};
        waveformPeriod = varargin{3};
    elseif isa(pxschObj,'nrPUSCHConfig')
        decodePusch = true;
        bwp = varargin{1};
        emissions = varargin{2};
        rxGridInBandEmissions = varargin{3};
    else
        waveformPeriod = varargin{1};
        decodePdcch = true;
    end

    evm3GPP = length(rxGrids)>1;
    rxGridLow = rxGrids{1};
    HestLow = hest{1};

    if evm3GPP
        rxGridHigh = rxGrids{2};
        HestHigh = hest{2};
    end

    nEVMWindowLocations = 1;
    if evm3GPP
        nEVMWindowLocations = 2;
    end
    eqSym = cell(1+evm3GPP,1);                  % Equalized symbols for constellation plot, for each low/high EVM window location
    refSym = cell(1+evm3GPP,1);                 % Reference symbols for constellation plot, for each low/high EVM window location
    rxSym = cell(1+evm3GPP,1);                  % Store transform precoded IQs, for each low/high EVM window

    initialNSlot = carrier.NSlot;
    cfgLen = length(cfgArray);
    ofdmInfo = nrOFDMInfo(carrier);
    L = ofdmInfo.SymbolsPerSlot;
    nSlots = floor(size(rxGridLow,2)/L);

    if decodePdsch
        nLayers = cfgArray(1).PDSCH.NumLayers;
        activeSlots = [];
        pdschConfigLen = length(cfgArray);
        for pIdx = 1:pdschConfigLen
            if ~isempty(cfgArray(pIdx).Resources)
                activeSlots = [activeSlots cfgArray(pIdx).Resources.NSlot]; %#ok<*AGROW>
            end
        end
        activeSlots = unique(activeSlots);
    elseif decodePusch
        nLayers = cfgArray.PUSCH.NumLayers;
        transformPrecoding = pxschObj.TransformPrecoding;
        txScheme = pxschObj.TransmissionScheme;
        prbSet = cfgArray.PUSCH.PRBSet;
        bwpLen = bwp.NSizeBWP;
        if iscolumn(prbSet)
            prbSet = prbSet.';
        end
        nDeltaRB = bwpLen - size(prbSet, 2);
        activeSlots = [cfgArray.Resources.NSlot];
    else
        activeSlots = [cfgArray.Resources.NSlot];
        nLayers = 1;
    end

    refSlotGrid = zeros(size(HestLow,1),(initialNSlot+nSlots)*L,nLayers,nEVMWindowLocations);
    eqSlotGrid = refSlotGrid;

    slotRange = activeSlots(activeSlots < initialNSlot+nSlots);
    slotRange = slotRange(slotRange >= initialNSlot);
    numSCs = size(rxGridLow,1);

    nEVMWindowLocations = 1;
    if evm3GPP
        nEVMWindowLocations = 2;
    end

    if decodePdsch && pdschEncodingOn
        decodeDLSCH = nrDLSCHDecoder;
        decodeDLSCH.LDPCDecodingAlgorithm = 'Normalized min-sum';
        dlsch = nrDLSCH('MultipleHARQProcesses',false);
    end

    for slotIdx=slotRange
        carrier.NSlot = slotIdx;
        if decodePdsch
            % Skip this slot if waveformPeriod is not valid and slot is not allocated
            if ~isinf(waveformPeriod) && isempty(find(arrayfun(@(x) mod(x,waveformPeriod), slotRange) == mod(slotIdx,waveformPeriod),1))
                continue;
            end
            [channelIndices,channelSymbols,~,~,ptrsIndices,ptrsSymbols] = hSlotResources(cfgArray,mod(slotIdx,waveformPeriod));
        elseif decodePusch
            cfgArray.PUSCH.TransmissionScheme = 'nonCodeBook';
            [channelIndices,dmrsIndices,~,ptrsIndices] = ...
                nr5g.internal.wavegen.PXSCHResources(bwp, cfgArray.PUSCH, cfgArray.PUSCH.NID, slotIdx, [], []);

            % Extract the transport block size, rv and bit capacity (G) for the relevant slot
            resourceIdx = find(slotIdx == activeSlots);
            rv = cfgArray.Resources(resourceIdx).RV;
            tcr = cfgArray.PUSCH.TargetCodeRate;

            % Get the Transport block size from the resource
            tbs = [];
            if rv == 0
                tbs = cfgArray.Resources(resourceIdx).TransportBlockSize;
            elseif cfgArray.PUSCH.Coding
                rv0ResIdx = find([cfgArray.Resources(:).RV] == 0,1);
                tbs = cfgArray.Resources(rv0ResIdx).TransportBlockSize;
            end
            G = cfgArray.Resources(resourceIdx).G;

            if cfgArray.PUSCH.Coding
                puschIndCurrentSlot = setdiff(channelIndices,ptrsIndices);
                puschIndCurrentSlot = reshape(puschIndCurrentSlot,[],nLayers);
            else
                puschIndCurrentSlot = channelIndices;
            end
        else
            loc = find(mod(slotIdx,waveformPeriod) == activeSlots);
            if isempty(loc)
                continue;
            end
            channelIndices = cfgArray.Resources(loc).ChannelIndices;
            channelSymbols = cfgArray.Resources(loc).ChannelSymbols;
            ptrsIndices = [];
            ptrsSymbols = [];
        end

        % Do not include first two slot symbols for PDSCH EVM (used for
        % control as specified in TS 38.141-1 table 4.9.2.2-2 (NR-TMs) /
        % TS 38.101-1 table A.3.1-1 (FRCs))
        % This step is not performed when
        % * extrapolateHest is true, or
        % * label is absent, or
        % * coding is enabled
        % Ensure that the symbol allocation does not include the first two
        % OFDM symbols when the 'label' field is set. A non-empty 'label'
        % is assumed to contain a string associated with Release 15
        % NR-TMs/FRCs
        if decodePdsch && skipCtrlSymbols  && ~pdschEncodingOn
            idx = channelIndices(:,1) <= (2*numSCs);
            channelIndices(idx,:) = [];
            channelSymbols(idx,:) = [];
        end

        % Extract the relevant slot, channel estimates and PDSCH allocated REs
        currentSlotIdx = slotIdx - initialNSlot;
        currentRxGridLow = rxGridLow(:, currentSlotIdx*L+(1:L), :);
        currentHestLow = HestLow(:, currentSlotIdx*L+(1:L),:,:);
        [pxschRxLow,pxschHestLow] = nrExtractResources(channelIndices,currentRxGridLow,currentHestLow);
        noiseEst = 0;         % ZF based equalization, as per TS 38.104, Annex B.1 (FR1) / Annex C.1 (FR2)

        % Skip if no valid channel estimates present in current slot
        if ~sum(abs(pxschHestLow),'All')
            continue;
        end
        [eqGridLow,csiLow] = nrEqualizeMMSE(pxschRxLow,pxschHestLow,noiseEst);
        if evm3GPP
            currentRxGridHigh = rxGridHigh(:, currentSlotIdx*L+(1:L), :);
            currentHestHigh = HestHigh(:, currentSlotIdx*L+(1:L),:,:);
            [pxschRxHigh,pxschHestHigh] = nrExtractResources(channelIndices,currentRxGridHigh,currentHestHigh);
            [eqGridHigh,csiHigh] = nrEqualizeMMSE(pxschRxHigh,pxschHestHigh,noiseEst);
        end

        % Estimate and compensate phase error
        % Process only active configurations, one per slot. If more
        % than one configuration is present, process the first one as
        % the indices (pdschIndices, dmrsIndices) collectively include
        % all active ones in current slot
        if decodePdsch
            for idx = 1:length(pxschObj)
                period = cfgArray(idx).PDSCH.Period;
                if isempty(period)
                    period = 0;
                end

                % For repeated waveforms, pick slot index within the
                % waveform period
                tmpSlotIdx = mod(slotIdx,waveformPeriod);

                if any(mod(tmpSlotIdx,period) == cfgArray(idx).PDSCH.SlotAllocation) && cfgArray(idx).PDSCH.EnablePTRS
                    eqGridLow = [];
                    eqGridHigh = [];
                    eqGridLow = [eqGridLow; pdschCPE(carrier,pxschObj{idx},channelIndices,ptrsIndices,ptrsSymbols,currentRxGridLow,currentHestLow)];
                    if evm3GPP
                        eqGridHigh = [eqGridHigh; pdschCPE(carrier,pxschObj{idx},channelIndices,ptrsIndices,ptrsSymbols,currentRxGridHigh,currentHestHigh)];
                    end
                    break;
                end
            end
        end

        if decodePusch
            if cfgArray.PUSCH.EnablePTRS
                pxschObj.TransformPrecoding = transformPrecoding;
                eqGridLow = puschCPE(carrier,pxschObj,currentRxGridLow,currentHestLow);
                if evm3GPP
                    eqGridHigh = puschCPE(carrier,pxschObj,currentRxGridHigh,currentHestHigh);
                end
            end
        end

        % Accumulate equalized IQs
        eqSym{1} = [eqSym{1}; eqGridLow];
        if evm3GPP
            eqSym{2} = [eqSym{2}; eqGridHigh];
        end

        % Create vectors needed for decoding each pdschObj
        if decodePdsch
            trBlkSizes = zeros(cfgLen,1);
            G = zeros(cfgLen,1);
            rv = zeros(cfgLen,1);
            for pIdx = 1:cfgLen
                period = cfgArray(pIdx).PDSCH.Period;
                if isempty(period)
                    period = 0;
                end

                % For repeated waveforms, pick slot index within the
                % waveform period
                tmpSlotIdx = mod(slotIdx,waveformPeriod);
                if any(mod(tmpSlotIdx,period) == cfgArray(pIdx).PDSCH.SlotAllocation) && cfgArray(pIdx).PDSCH.Enable
                    currentSlotRange = [];
                    currentSlotRange = unique([currentSlotRange cfgArray(pIdx).Resources.NSlot]);

                    % Locate the resource index corresponding to the allocated slot
                    resourceIdx = find(tmpSlotIdx == currentSlotRange);
                    trBlkSizes(pIdx) = cfgArray(pIdx).Resources(resourceIdx).TransportBlockSize;
                    G(pIdx) = cfgArray(pIdx).Resources(resourceIdx).G;
                    if ~isempty(cfgArray(pIdx).Resources(resourceIdx).RV)
                        rv(pIdx) = cfgArray(pIdx).Resources(resourceIdx).RV;
                    end
                end
            end
        end

        % For low edge EVM and high edge EVM
        for e = 1:nEVMWindowLocations
            % Select the low or high edge equalizer output
            if (e == 1)
                eqGrid = eqGridLow;
                csi = csiLow;
            else
                eqGrid = eqGridHigh;
                csi = csiHigh;
            end
            if decodePdcch || (decodePusch && ~cfgArray.PUSCH.Coding)
                rxSymbols = eqGrid;
                if decodePdcch
                    % Decode the PDCCH soft output
                    pdcchLLRs = nrPDCCHDecode(eqGrid,carrier.NCellID,cfgArray.PDCCH.RNTI);

                    % Obtain reference symbols using hard slicing of the PDCCH LLRs
                    channelSymbols = nrPDCCH(double(pdcchLLRs<0),carrier.NCellID,cfgArray.PDCCH.RNTI);
                end
                if decodePusch
                    channelSymbols = cfgArray.Resources(resourceIdx).ChannelSymbols;
                    if strcmp(txScheme,'codebook')
                        W = nrPUSCHCodebook(nLayers,cfgArray.PUSCH.NumAntennaPorts,cfgArray.PUSCH.TPMI,cfgArray.PUSCH.TransformPrecoding);
                        channelSymbols = channelSymbols * pinv(W);
                    end
                    if pxschObj.TransformPrecoding
                        rxSymbols = nrTransformDeprecode(rxSymbols,length(pxschObj.PRBSet));
                        channelSymbols = nrTransformDeprecode(channelSymbols,length(pxschObj.PRBSet));
                    end
                end
            else
                if decodePdsch
                    rx = cell(1);
                    channelSymbols = [];
                    eqOffset = 0;
                    for pIdx = 1:cfgLen
                        % Decode rxSymbols. If CRC passes, re-encode them for
                        % obtaining reference IQs for EVM calculation. Extract
                        % the equalized symbols and csi for each pdschObj.
                        % Skip processing if not allocated in current slot.

                        period = cfgArray(pIdx).PDSCH.Period;
                        if isempty(period)
                            period = 0;
                        end

                        % Skip if slot is not allocated (taking into
                        % account waveform period)
                        if ~any(mod(tmpSlotIdx,period) == cfgArray(pIdx).PDSCH.SlotAllocation) || ~cfgArray(pIdx).PDSCH.Enable
                            continue;
                        end

                        if pdschEncodingOn
                            eqLen = length(cfgArray(pIdx).Resources(resourceIdx).ChannelIndices);
                        else
                            cIndices = cfgArray(pIdx).Resources(resourceIdx).ChannelIndices;
                            if skipCtrlSymbols
                                skipIdx = cIndices(:,1) <= (2*numSCs);
                                cIndices(skipIdx,:) = [];
                            end
                            eqLen = length(cIndices);
                        end
                        eq = eqGrid(1+eqOffset:eqOffset+eqLen,:);
                        currentCsi = csi(1+eqOffset:eqOffset+eqLen,:);
                        eqOffset = eqOffset+ eqLen;

                        % Remove power scaling before LLR computation
                        eq = eq.*db2mag(-1*cfgArray(pIdx).PDSCH.Power);
                        [dlschLLRs,rxSymbols] = nrPDSCHDecode(eq,pxschObj{pIdx}.Modulation,pxschObj{pIdx}.NID,pxschObj{pIdx}.RNTI,noiseEst);

                        % Adjust for channel power
                        rxSymbols{1} = rxSymbols{1}.*db2mag(cfgArray(pIdx).PDSCH.Power);

                        % Scale LLRs by CSI
                        numCWs = size(dlschLLRs,2);
                        currentCsi = nrLayerDemap(currentCsi);    % CSI layer demapping
                        for cwIdx = 1:numCWs
                            Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); % bits per symbol
                            currentCsi{cwIdx} = repmat(currentCsi{cwIdx}.',Qm,1);   % expand by each bit per symbol
                            dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* currentCsi{cwIdx}(:);   % scale
                        end

                        if pdschEncodingOn && trBlkSizes(pIdx)
                            decodeDLSCH.TransportBlockLength = trBlkSizes(pIdx);
                            decodeDLSCH.TargetCodeRate = cfgArray(pIdx).PDSCH.TargetCodeRate;
                            [demodBits,blkerr] = decodeDLSCH(dlschLLRs,pxschObj{pIdx}.Modulation,nLayers,rv(pIdx));
                            if any(blkerr)
                                warning('CRC failed on decoded data, using sliced received symbols, EVM may be inaccurate.');
                                % Attempt to gracefully handle bad CRC. Generate some
                                % bits for this slot by hard slicing LLRs, so that the
                                % simulation can continue
                                recodedBits = cellfun(@(x) x<0, dlschLLRs, 'UniformOutput', false);
                            else
                                trBlk = demodBits;
                                dlsch.TargetCodeRate = decodeDLSCH.TargetCodeRate;
                                setTransportBlock(dlsch,trBlk);
                                recodedBits = dlsch(pxschObj{pIdx}.Modulation,nLayers,G(pIdx),cfgArray(pIdx).PDSCH.RVSequence(1));
                            end
                            ref = nrPDSCH(recodedBits,pxschObj{pIdx}.Modulation,nLayers,pxschObj{pIdx}.NID,pxschObj{pIdx}.RNTI);
                            ref = ref*db2mag(cfgArray(pIdx).PDSCH.Power);
                        else
                            % Obtain reference symbols using hard slicing of the PDSCH LLRs
                            pdschSym = nrPDSCH(double(dlschLLRs{1}<0),pxschObj{pIdx}.Modulation,nLayers,pxschObj{pIdx}.NID,pxschObj{pIdx}.RNTI);
                            [~, ref] = nrPDSCHDecode(pdschSym,pxschObj{pIdx}.Modulation,pxschObj{pIdx}.NID,pxschObj{pIdx}.RNTI,noiseEst);
                            ref = ref{1};

                            % Adjust reference symbol power
                            ref = ref*db2mag(cfgArray(pIdx).PDSCH.Power);
                            if pdschEncodingOn && ~trBlkSizes(pIdx) || ~pdschEncodingOn
                                % Map the reference and rxSymbols to the
                                % relevant layers
                                ref = nrLayerMap(ref,nLayers);
                                rxSymbols = nrLayerMap(rxSymbols,nLayers);
                            end
                        end
                        if iscell(rxSymbols)
                            rxSymbols = cell2mat(rxSymbols);
                        end
                        rx{1} = [rx{1} ; rxSymbols];
                        channelSymbols = [channelSymbols;ref];
                    end
                elseif decodePusch
                    pxschObj.TransmissionScheme = 'nonCodebook';
                    pxschObj.TransformPrecoding = transformPrecoding;
                    [ulschLLRs,rxSymbols] = nrPUSCHDecode(carrier,pxschObj,eqGrid,noiseEst);

                    if pxschObj.TransformPrecoding
                        if pxschObj.EnablePTRS
                            % Remove csi values corresponding to PT-RS indices
                            tmpPtrsInd = [];
                            for ptrsIdx = 1:length(ptrsIndices)
                                tmpPtrsInd = [tmpPtrsInd find(ptrsIndices(ptrsIdx) == channelIndices)];
                            end
                            csi(tmpPtrsInd) = [];
                        else
                            MRB = length(pxschObj.PRBSet);
                            MSC = MRB * 12;
                            csi = nrTransformDeprecode(csi,MRB) / sqrt(MSC);
                            csi = repmat(csi((1:MSC:end).'),1,MSC).';
                            csi = reshape(csi,size(rxSymbols));
                        end
                    end

                    % Scale LLRs by CSI
                    csi = nrLayerDemap(csi);                      % CSI layer demapping
                    Qm = length(ulschLLRs) / length(rxSymbols);   % bits per symbol
                    csi = reshape(repmat(csi{1}.',Qm,1),[],1);    % expand by each bit per symbol
                    ulschLLRs = ulschLLRs .* real(csi);           % scale

                    decUL = nrULSCHDecoder('TargetCodeRate',tcr,'LDPCDecodingAlgorithm','Normalized min-sum',....
                        'MultipleHARQProcesses',false,'TransportBlockLength',tbs);
                    [demodBits,blkerr] = decUL(ulschLLRs,pxschObj.Modulation,nLayers,rv);
                    if any(blkerr)
                        warning('CRC failed on decoded data, using sliced received symbols, EVM may be inaccurate.');
                        % Attempt to gracefully handle bad CRC. Generate some
                        % bits for this slot by hard slicing LLRs, so that the
                        % simulation can continue
                        recodedBits = double(ulschLLRs<0);
                    else
                        encodeULSCH = nrULSCH;
                        encodeULSCH.TargetCodeRate = tcr;
                        trBlk = demodBits;
                        setTransportBlock(encodeULSCH,trBlk);
                        recodedBits = encodeULSCH(pxschObj.Modulation,nLayers,G,rv);
                    end

                    % Generate reference symbols using 'recodedBits'
                    pxschObj.TransformPrecoding = 0;
                    channelSymbols = nrPUSCH(carrier,pxschObj,recodedBits);
                end
            end
            if decodePdsch
                rxSymbols = cell2mat(rx);
            end

            if iscell(rxSymbols)
                rxSymbols = cell2mat(rxSymbols);
            end
            refSym{e} = [refSym{e}; channelSymbols];

            if decodePusch && transformPrecoding
                rxSym{e} = [rxSym{e}; rxSymbols];
            end

            % Map each reference-slot , equalized slot to its correct position in grid, layer, EVM-edge
            rbsPerSlot = numSCs*carrier.SymbolsPerSlot;
            tmpGrid = zeros(numSCs,carrier.SymbolsPerSlot);
            if (decodePdsch && pdschEncodingOn && any(trBlkSizes)) || (decodePusch && cfgArray.PUSCH.Coding)
                rxSymbols = reshape(rxSymbols,nLayers,length(rxSymbols)/nLayers).';
            end
            if decodePusch
                ind = zeros(size(puschIndCurrentSlot,1),nLayers);
            else
                ind = zeros(size(channelIndices,1),nLayers);
            end

            for layerIdx = 1:nLayers
                if layerIdx <= size(channelIndices,2)
                    if decodePusch
                        ind(:,layerIdx) = puschIndCurrentSlot(:,layerIdx) - rbsPerSlot*(layerIdx -1);
                    else
                        ind(:,layerIdx) = channelIndices(:,layerIdx) - rbsPerSlot*(layerIdx -1);
                    end
                    tmpGrid(ind(:,layerIdx)) = channelSymbols(:,layerIdx);
                    symIdx = (slotIdx)*L+1:(slotIdx+1)*L;
                    refSlotGrid(:,symIdx,layerIdx,e) = tmpGrid;
                    tmpGrid(:) = 0;
                    tmpGrid(ind(:,layerIdx)) = rxSymbols(:,layerIdx);
                    eqSlotGrid(:,symIdx,layerIdx,e) = tmpGrid;
                    tmpGrid(:) = 0;
                end
            end

            % In-band emissions measurement
            if decodePusch && ~isempty(emissions.DeltaRB)

                symIdx = (currentSlotIdx)*L+1:(currentSlotIdx+1)*L;
                symIdx(symIdx>size(rxGridInBandEmissions, 2)) = [];

                % YUnallocated contains REs outside the allocated region
                YUnallocated = rxGridInBandEmissions(:,symIdx,:);

                % YAllocated contains allocated REs
                YAllocated = zeros(size(YUnallocated));
                YAllocated([puschIndCurrentSlot; dmrsIndices]) = YUnallocated([puschIndCurrentSlot; dmrsIndices]);

                % Remove allocated RBs from YUnallocated
                YUnallocated([puschIndCurrentSlot; dmrsIndices]) = NaN;

                % Calculate the subcarriers indices spanning the PUSCH
                % allocated RBs, scInd. These subcarriers are removed
                % before in-band emission computation
                scInd = [];
                for rbIdx = 1:length(pxschObj.PRBSet)
                    rb = pxschObj.PRBSet(rbIdx);
                    scInd = [scInd rb*12+(1:12)];
                end

                % Compute absolute and relative emissions in
                % unallocated RBs for each slot
                scPower = sum(abs(YUnallocated(:,1:L)).^2,2,'omitnan');
                signalPower = sum(abs(YAllocated(:,1:L)).^2,2,'omitnan');
                scPower(scInd) = [];
                signalPower = signalPower(scInd);
                rbPower = sum(reshape(scPower, 12, nDeltaRB).', 2); %#ok<UDIM>
                signalRbPower = sum(signalPower)/(L*length(pxschObj.PRBSet));
                emissions.Absolute(:, slotIdx+1) = rbPower/(L);
                emissions.Relative(:, slotIdx+1) = emissions.Absolute(:, slotIdx+1)/signalRbPower;
            end
        end
    end
    if decodePusch && transformPrecoding
        eqSym = rxSym;
    end
end

function puschEq = puschCPE(carrier,pusch,rxGrid,estChannelGrid)
    % PUSCHEQ = puschCPE(CARRIER,PUSCH,RXGRID,ESTCHANNELGRID)
    % Slot-based Common phase error (CPE) estimation and correction for
    % PUSCH returns an equalized set of IQs, PUSCHEQ for which the CPE has
    % been estimated and compensated. PUSCHEQ is array spanning the length
    % of the PUSCH symbols x number of layers N.
    % PUSCH is an object of type 'nrPUSCHConfig'.
    % RXGRID is N-by-S-by-L grid extracted from the OFDM demodulated
    % waveform. S is the number of symbols in a slot.
    % ESTCHANNELGRID is N-by-S-by-P-by-R grid. P is the number of antenna
    % ports.

    % The channel estimation function must take in reference symbols
    % without applying codebook, to remove the overall channel (having
    % channel & precoder matrix). Generate PT-RS symbols for nonCodebook
    % transmission scheme before passing it to nrChannelEstimate function
    pusch.TransmissionScheme = 'nonCodebook';
    carrier.NSizeGrid = pusch.NSizeBWP;
    carrier.NStartGrid = pusch.NStartBWP;
    pusch.NSizeBWP = [];
    pusch.NStartBWP = [];

    % Generate reference PT-RS symbols, indices and puschIndicesInfo
    ptrsSymbols = nrPUSCHPTRS(carrier,pusch);
    [puschIndices, puschIndicesInfo, ptrsAllIndices] = nrPUSCHIndices(carrier,pusch);

    % Initialize temporary grid to store equalized symbols
    tempGrid = nrResourceGrid(carrier,pusch.NumLayers);

    % Extract PUSCH and PT-RS resource elements
    [puschRx, puschHest, ~, ~, ~, puschIndLayer] = nrExtractResources([puschIndices(:);ptrsAllIndices(:)],rxGrid,estChannelGrid,tempGrid);

    % Equalization
    % Set noiseEst to zero as it was not estimated along with the
    % corresponding estChannelGrid
    noiseEst = 0;
    puschEq = nrEqualizeMMSE(puschRx,puschHest,noiseEst);

    if pusch.TransformPrecoding
        % Use transform deprecoded symbols for estimating phase error
        puschRx = nrTransformDeprecode(puschEq,length(pusch.PRBSet));
        betaPTRS = nr5g.internal.pusch.ptrsPowerFactorDFTsOFDM(pusch.Modulation);
    else
        puschRx = puschEq;
        betaPTRS = 1;
    end
    tempGrid(puschIndLayer) = puschRx;

    % Estimate the residual channel at the PT-RS locations in
    % tempGrid.
    cpe = nrChannelEstimate(tempGrid,ptrsAllIndices,betaPTRS*ptrsSymbols,'CyclicPrefix',carrier.CyclicPrefix);

    % Sum estimates across subcarriers, receive antennas, and
    % layers. Then, get the CPE by taking the angle of the
    % resultant sum
    cpe = angle(sum(cpe,[1 3 4]));

    % Update the tempGrid with equalized symbols
    tempGrid(puschIndLayer) = puschEq;

    % Correct CPE in each OFDM symbol within the range of reference
    % PT-RS OFDM symbols
    if numel(puschIndicesInfo.PTRSSymbolSet) > 0
        symLoc = puschIndicesInfo.PTRSSymbolSet(1)+1:puschIndicesInfo.PTRSSymbolSet(end)+1;
        tempGrid(:,symLoc,:) = tempGrid(:,symLoc,:).*exp(-1i*cpe(symLoc));
    end

    % Extract PUSCH symbols to be passed in to the nrPUSCHDecode function
    [~,puschInd1] = nrExtractResources(puschIndices,tempGrid);
    puschEq = tempGrid(puschInd1);
end

function pdschEq = pdschCPE(carrier,pdsch,pdschIndices,ptrsIndices,ptrsSymbols,rxGrid,estChannelGrid)
%PDSCHEQ = pdschCPE(CARRIER,PDSCH,PDSCHINDICES,PTRSINDICES,PTRSSYMBOLS,RXGRID,ESTCHANNELGRID)
    % Slot-based Common phase error (CPE) estimation and correction for
    % PDSCH returns an equalized set of IQs PDSCHEQ, for which the CPE has
    % been estimated and compensated. PDSCHEQ is an array spanning the length
    % of the PDSCH symbols x number of layers N. 
    % PDSCH is an obj of type 'nrPDSCHConfig'.
    % PDSCHINDICES is an array of S-by-L, where S is the
    % length of the PDSCH allocation and L is the number of layers.
    % PTRSINDICES is an array of S-by-L, where S is the length of the DM-RS
    % allocation and L is the number of layers.
    % PTRSSYMBOLS is an array of S-by-L, where S is the length of the DM-RS
    % allocation and L is the number of layers.
    % RXGRID is a N-by-S-by-L grid extracted from the OFDM demodulated
    % waveform. S is the number of symbols in a slot.
    % ESTCHANNELGRID is N-by-S-by-P-by-R grid. P is the number of antenna
    % ports
    
    [~,pdschIndicesInfo] = nrPDSCHIndices(carrier,pdsch);

    % Set noiseEst to zero as it was not estimated along with the
    % corresponding estChannelGrid
    noiseEst = 0;

    % Get PDSCH resource elements from the received grid
    [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChannelGrid);

    % Equalization
    pdschEq = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

    % Initialize temporary grid to store equalized symbols
    carrier.NSizeGrid = pdsch.NSizeBWP;
    tempGrid = nrResourceGrid(carrier,pdsch.NumLayers);

    % Extract PT-RS symbols from received grid and estimated
    % channel grid
    [ptrsRx,ptrsHest,~,~,~,ptrsLayerIndices] = nrExtractResources(ptrsIndices,rxGrid,estChannelGrid,tempGrid);

    % Equalize PT-RS symbols and map them to tempGrid
    ptrsEq = nrEqualizeMMSE(ptrsRx,ptrsHest,noiseEst);
    tempGrid(ptrsLayerIndices) = ptrsEq;

    % Estimate the residual channel at the PT-RS locations in
    % tempGrid
    cpe = nrChannelEstimate(tempGrid,ptrsIndices,ptrsSymbols,'CyclicPrefix',carrier.CyclicPrefix);

    % Sum estimates across subcarriers, receive antennas, and
    % layers. Then, get the CPE by taking the angle of the
    % resultant sum
    cpe = angle(sum(cpe,[1 3 4]));

    % Map the equalized PDSCH symbols to tempGrid
    tempGrid(pdschIndices) = pdschEq;

    % Correct CPE in each OFDM symbol within the range of reference
    % PT-RS OFDM symbols
    if numel(pdschIndicesInfo.PTRSSymbolSet) > 0
        symLoc = pdschIndicesInfo.PTRSSymbolSet(1)+1:pdschIndicesInfo.PTRSSymbolSet(end)+1;
        tempGrid(:,symLoc,:) = tempGrid(:,symLoc,:).*exp(-1i*cpe(symLoc));
    end

    % Extract PDSCH symbols
    pdschEq = tempGrid(pdschIndices);
end