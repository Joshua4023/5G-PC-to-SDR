function [refGrid,idealGrid] = hReferenceGrid(carrier,bwpCfg,configArray,nSlots,phyChannelType,varargin)
    %hReferenceGrid Create a DM-RS based reference grid for a given number
    %   of slots. It also returns a grid of known data IQs.
    %   [REFGRID,IDEALGRID] = hReferenceGrid(CARRIER,BWPCFG,CONFIGARRAY,...
    %                                               NSLOTS,PHYCHANNELTYPE)
    %   creates a DM-RS based reference grid and known data based grid for a
    %   given number of slots and physical channel type.
    %   REFGRID contains the DM-RS symbols specified in CONFIGARRAY and
    %   is of dimensions K-by-N-by-L, where K is number of subcarriers of
    %   size carrier.NSizeGrid*12, N is the number of symbols spanning
    %   nSlots and L is the number of layers.
    %   IDEALGRID contains known data IQ symbols specified in CONFIGARRAY
    %   and is of dimensions K-by-N-by-L.
    %   CARRIER is a configuration object of type, <a
    %                        href="matlab:help('nrCarrierConfig')"
    %                        >nrCarrierConfig</a>
    %   BWPCFG is a configuration object of type, <a
    %                        href="matlab:help('nrWavegenBWPConfig')"
    %                        >nrWavegenBWPConfig</a>
    %   CONFIGARRAY is an array of structures. Each structure contains the
    %   configuration and reference resource information. The contents of
    %   this structure are used for equalization, decoding, and building
    %   reference IQs.
    %   NSLOTS is the number of slots for which the REFGRID is to be
    %   generated
    %   PHYCHANNELTYPE is a character string with three possible strings
    %   ('PDSCH','PDCCH',PUSCH')
    %   [REFGRID,IDEALGRID] = hReferenceGrid(...,WAVEFORMPERIOD) creates
    %   a reference grid of DM-RS and known data taking into account the
    %   waveform repetition duration WAVEFORMPERIOD. WAVEFORMPERIOD is
    %   in duration of slots and is relevant only to PDSCH and PDCCH.

    % Copyright 2022-2023 The MathWorks, Inc.
    

    genPdschRefGrid = strcmp(phyChannelType,'PDSCH');
    genPdcchRefGrid = strcmp(phyChannelType,'PDCCH');
    genPuschRefGrid = strcmp(phyChannelType,'PUSCH');

    if isempty(configArray)
        refGrid = [];
        idealGrid = [];
        return;
    end

    nSubcarriers = carrier.NSizeGrid * 12;
    L = carrier.SymbolsPerSlot*nSlots; % Number of OFDM symbols in the reference grid
    if genPdschRefGrid
        nLayers = configArray(1).PDSCH.NumLayers;
    elseif genPuschRefGrid
        nLayers = configArray.PUSCH.NumLayers;
    else
        nLayers = 1;
    end
    bwpStart = bwpCfg.NStartBWP-carrier.NStartGrid;
    bwpLen = bwpCfg.NSizeBWP;

    if genPdschRefGrid
        slotRange = [];
        pdschConfigLen = numel(configArray);
        for pIdx = 1:pdschConfigLen
            if ~isempty(configArray(pIdx).Resources)
                slotRange = [slotRange configArray(pIdx).Resources.NSlot]; %#ok<*AGROW>
            end
        end
        slotRange = unique(slotRange);
    else
        slotRange = [configArray.Resources.NSlot];
    end

    if genPdschRefGrid || genPdcchRefGrid
        waveformPeriod = inf;
        if nargin == 6
            waveformPeriod = varargin{1};
        end
    end
    refGrid = zeros(nSubcarriers,L,nLayers); % empty grid
    bwpGrid = zeros(bwpLen*12,L,nLayers);
    resPerSlot = bwpLen*12*carrier.SymbolsPerSlot;

    % Create 'idealGrid'
    if genPdschRefGrid || genPuschRefGrid
    	idealGrid = zeros(nSubcarriers,L,nLayers);
        idealBwpGrid = zeros(bwpLen*12,L,nLayers);
    end

    % Populate the DM-RS symbols in the reference grid for all slots
    % Place the bwpGrid in a carrier grid (at an appropriate
    % location) in case the BWP size is not the same as the carrier grid
    for slotIdx=carrier.NSlot+(0:nSlots-1)
        isDataSlot = ~isempty(find(slotIdx == slotRange,1));
        if genPdschRefGrid
            [dataInd,dataSym,dmrsIndices,dmrsSymbols] = hSlotResources(configArray,mod(slotIdx,waveformPeriod));
        elseif genPuschRefGrid
            configArray.PUSCH.TransmissionScheme = 'nonCodeBook';
            [~,dmrsIndices,dmrsSymbols,~,~,~] = ...
                nr5g.internal.wavegen.PXSCHResources(bwpCfg, configArray.PUSCH, carrier.NCellID, slotIdx, [], []);
            if isDataSlot
                dIdx = find(slotIdx == slotRange,1);
                dataInd = configArray.Resources(dIdx).ChannelIndices;
                dataSym = configArray.Resources(dIdx).ChannelSymbols;
            end
        elseif genPdcchRefGrid
            loc = find(mod(slotIdx,waveformPeriod) == slotRange);
            if isempty(loc)
                continue;
            end
            dmrsIndices = configArray.Resources(loc).DMRSIndices;
            dmrsSymbols = configArray.Resources(loc).DMRSSymbols;
        end

        if ~isempty(dmrsIndices)
            for layerIdx = 1:nLayers
                if layerIdx <= size(dmrsIndices,2)
                    dmrsIndices(:,layerIdx) = dmrsIndices(:,layerIdx) - resPerSlot*(layerIdx -1)+ (L*bwpLen*12*(layerIdx-1));

                    % Power adjustment is needed only for PUSCH case 
                    pwr = 1;
                    if genPuschRefGrid
                        pwr = db2mag(configArray.PUSCH.DMRSPower);
                    end
                    bwpGrid(dmrsIndices(:,layerIdx)+(slotIdx-carrier.NSlot)*resPerSlot) = dmrsSymbols(:,layerIdx)*pwr;
                    if isDataSlot && (genPdschRefGrid || genPuschRefGrid)
                        dataInd(:,layerIdx) = dataInd(:,layerIdx) - resPerSlot*(layerIdx -1)+ (L*bwpLen*12*(layerIdx-1));
                        idealBwpGrid(dataInd(:,layerIdx)+(slotIdx-carrier.NSlot)*resPerSlot) = dataSym(:,layerIdx);
                    end
                end
            end
            refGrid(12*bwpStart+1:12*(bwpStart+bwpLen),:,:) = bwpGrid;
            if genPdschRefGrid || genPuschRefGrid
                idealGrid(12*bwpStart+1:12*(bwpStart+bwpLen),:,:) = idealBwpGrid;
            end
        end
    end
end
