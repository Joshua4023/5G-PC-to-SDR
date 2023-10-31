% Generated by MATLAB(R) 23.2 (R2023b) and 5G Toolbox 23.2 (R2023b).
% Generated on: 31-Oct-2023 14:56:17

%% Generating NR Test Models waveform
% NR Test Models configuration
cfgDLTM = nrDLCarrierConfig;
cfgDLTM.Label = 'NR-FR1-TM1.1';
cfgDLTM.FrequencyRange = 'FR1';
cfgDLTM.ChannelBandwidth = 40;
cfgDLTM.NCellID = 1;
cfgDLTM.NumSubframes = 20;
cfgDLTM.InitialNSubframe = 0;
cfgDLTM.WindowingPercent = 0;
cfgDLTM.SampleRate = [];
cfgDLTM.CarrierFrequency = 0;

%% SCS specific carriers
scscarrier = nrSCSCarrierConfig;
scscarrier.SubcarrierSpacing = 15;
scscarrier.NSizeGrid = 216;
scscarrier.NStartGrid = 0;

cfgDLTM.SCSCarriers = {scscarrier};

%% Bandwidth Parts
bwp = nrWavegenBWPConfig;
bwp.BandwidthPartID = 1;
bwp.Label = 'BWP1';
bwp.SubcarrierSpacing = 15;
bwp.CyclicPrefix = 'normal';
bwp.NSizeBWP = 216;
bwp.NStartBWP = 0;

cfgDLTM.BandwidthParts = {bwp};

%% Synchronization Signals Burst
ssburst = nrWavegenSSBurstConfig;
ssburst.Enable = false;
ssburst.Power = 0;
ssburst.BlockPattern = 'Case A';
ssburst.TransmittedBlocks = [1 0 0 0];
ssburst.Period = 10;
ssburst.NCRBSSB = [];
ssburst.KSSB = 0;
ssburst.DataSource = 'MIB';
ssburst.DMRSTypeAPosition = 2;
ssburst.CellBarred = false;
ssburst.IntraFreqReselection = false;
ssburst.PDCCHConfigSIB1 = 0;
ssburst.SubcarrierSpacingCommon = 30;

cfgDLTM.SSBurst = ssburst;

%% CORESET and Search Space Configuration
coreset = nrCORESETConfig;
coreset.CORESETID = 1;
coreset.Label = 'CORESET1';
coreset.FrequencyResources = 1;
coreset.Duration = 2;
coreset.CCEREGMapping = 'noninterleaved';
coreset.REGBundleSize = 2;
coreset.InterleaverSize = 2;
coreset.ShiftIndex = 0;
coreset.PrecoderGranularity = 'sameAsREG-bundle';
coreset.RBOffset = [];

cfgDLTM.CORESET = {coreset};

% Search Spaces
searchspace = nrSearchSpaceConfig;
searchspace.SearchSpaceID = 1;
searchspace.Label = 'SearchSpace1';
searchspace.CORESETID = 1;
searchspace.SearchSpaceType = 'common';
searchspace.StartSymbolWithinSlot = 0;
searchspace.SlotPeriodAndOffset = [1 0];
searchspace.Duration = 1;
searchspace.NumCandidates = [8 8 0 0 0];

cfgDLTM.SearchSpaces = {searchspace};

%% PDCCH Instances Configuration
pdcch = nrWavegenPDCCHConfig;
pdcch.Enable = true;
pdcch.Label = 'PDCCH1';
pdcch.Power = 0;
pdcch.BandwidthPartID = 1;
pdcch.SearchSpaceID = 1;
pdcch.AggregationLevel = 1;
pdcch.AllocatedCandidate = 1;
pdcch.CCEOffset = [];
pdcch.SlotAllocation = 0:3;
pdcch.Period = 5;
pdcch.Coding = false;
pdcch.DataBlockSize = 20;
pdcch.DataSource = 'PN23';
pdcch.RNTI = 0;
pdcch.DMRSScramblingID = [];
pdcch.DMRSPower = 0;

cfgDLTM.PDCCH = {pdcch};

%% PDSCH Instances Configuration
% PDSCH 1
pdsch1 = nrWavegenPDSCHConfig;
pdsch1.Enable = true;
pdsch1.Label = 'Partial band PDSCH sequence with QPSK modulation scheme (target, RNTI = 0) (Full downlink slots)';
pdsch1.Power = 0;
pdsch1.BandwidthPartID = 1;
pdsch1.Modulation = 'QPSK';
pdsch1.NumLayers = 1;
pdsch1.MappingType = 'A';
pdsch1.ReservedCORESET = [];
pdsch1.SymbolAllocation = [0 14];
pdsch1.SlotAllocation = 0:2;
pdsch1.Period = 5;
pdsch1.PRBSet = 3:215;
pdsch1.PRBSetType = 'VRB';
pdsch1.VRBToPRBInterleaving = false;
pdsch1.VRBBundleSize = 2;
pdsch1.NID = [];
pdsch1.RNTI = 0;
pdsch1.Coding = false;
pdsch1.TargetCodeRate = 0.4785;
pdsch1.TBScaling = 1;
pdsch1.XOverhead = 0;
pdsch1.RVSequence = 0;
pdsch1.DataSource = 'PN23';
pdsch1.DMRSPower = 0;
pdsch1.EnablePTRS = false;
pdsch1.PTRSPower = 0;

% PDSCH Reserved PRB
pdsch1ReservedPRB = nrPDSCHReservedConfig;
pdsch1ReservedPRB.PRBSet = 0:2;
pdsch1ReservedPRB.SymbolSet = [0 1];
pdsch1ReservedPRB.Period = 1;

pdsch1.ReservedPRB = {pdsch1ReservedPRB};

% PDSCH DM-RS
pdsch1DMRS = nrPDSCHDMRSConfig;
pdsch1DMRS.DMRSConfigurationType = 1;
pdsch1DMRS.DMRSReferencePoint = 'CRB0';
pdsch1DMRS.DMRSTypeAPosition = 2;
pdsch1DMRS.DMRSAdditionalPosition = 1;
pdsch1DMRS.DMRSLength = 1;
pdsch1DMRS.CustomSymbolSet = [];
pdsch1DMRS.DMRSPortSet = [];
pdsch1DMRS.NIDNSCID = [];
pdsch1DMRS.NSCID = 0;
pdsch1DMRS.NumCDMGroupsWithoutData = 1;
pdsch1DMRS.DMRSDownlinkR16 = false;

pdsch1.DMRS = pdsch1DMRS;

% PDSCH PT-RS
pdsch1PTRS = nrPDSCHPTRSConfig;
pdsch1PTRS.TimeDensity = 4;
pdsch1PTRS.FrequencyDensity = 2;
pdsch1PTRS.REOffset = '00';
pdsch1PTRS.PTRSPortSet = [];

pdsch1.PTRS = pdsch1PTRS;

% PDSCH 2
pdsch2 = nrWavegenPDSCHConfig;
pdsch2.Enable = true;
pdsch2.Label = 'Partial band PDSCH sequence with QPSK modulation scheme (target, RNTI = 2) (Full downlink slots)';
pdsch2.Power = 0;
pdsch2.BandwidthPartID = 1;
pdsch2.Modulation = 'QPSK';
pdsch2.NumLayers = 1;
pdsch2.MappingType = 'A';
pdsch2.ReservedCORESET = [];
pdsch2.SymbolAllocation = [2 12];
pdsch2.SlotAllocation = 0:2;
pdsch2.Period = 5;
pdsch2.PRBSet = 0:2;
pdsch2.PRBSetType = 'VRB';
pdsch2.VRBToPRBInterleaving = false;
pdsch2.VRBBundleSize = 2;
pdsch2.NID = [];
pdsch2.RNTI = 2;
pdsch2.Coding = false;
pdsch2.TargetCodeRate = 0.4785;
pdsch2.TBScaling = 1;
pdsch2.XOverhead = 0;
pdsch2.RVSequence = 0;
pdsch2.DataSource = 'PN23';
pdsch2.DMRSPower = 0;
pdsch2.EnablePTRS = false;
pdsch2.PTRSPower = 0;

pdschreserved = nrPDSCHReservedConfig;
pdschreserved.PRBSet = 0:2;
pdschreserved.SymbolSet = [0 1];
pdschreserved.Period = 1;

pdsch2.ReservedPRB = {pdschreserved};

pdschdmrs = nrPDSCHDMRSConfig;
pdschdmrs.DMRSConfigurationType = 1;
pdschdmrs.DMRSReferencePoint = 'CRB0';
pdschdmrs.DMRSTypeAPosition = 2;
pdschdmrs.DMRSAdditionalPosition = 1;
pdschdmrs.DMRSLength = 1;
pdschdmrs.CustomSymbolSet = [];
pdschdmrs.DMRSPortSet = [];
pdschdmrs.NIDNSCID = [];
pdschdmrs.NSCID = 0;
pdschdmrs.NumCDMGroupsWithoutData = 1;
pdschdmrs.DMRSDownlinkR16 = false;

pdsch2.DMRS = pdschdmrs;

pdschptrs = nrPDSCHPTRSConfig;
pdschptrs.TimeDensity = 4;
pdschptrs.FrequencyDensity = 2;
pdschptrs.REOffset = '00';
pdschptrs.PTRSPortSet = [];

pdsch2.PTRS = pdschptrs;

% PDSCH 3
pdsch3 = nrWavegenPDSCHConfig;
pdsch3.Enable = true;
pdsch3.Label = 'Partial band PDSCH sequence with QPSK modulation scheme (target, RNTI = 0) (Partial downlink slots)';
pdsch3.Power = 0;
pdsch3.BandwidthPartID = 1;
pdsch3.Modulation = 'QPSK';
pdsch3.NumLayers = 1;
pdsch3.MappingType = 'A';
pdsch3.ReservedCORESET = [];
pdsch3.SymbolAllocation = [0 10];
pdsch3.SlotAllocation = 3;
pdsch3.Period = 5;
pdsch3.PRBSet = 3:215;
pdsch3.PRBSetType = 'VRB';
pdsch3.VRBToPRBInterleaving = false;
pdsch3.VRBBundleSize = 2;
pdsch3.NID = [];
pdsch3.RNTI = 0;
pdsch3.Coding = false;
pdsch3.TargetCodeRate = 0.4785;
pdsch3.TBScaling = 1;
pdsch3.XOverhead = 0;
pdsch3.RVSequence = 0;
pdsch3.DataSource = 'PN23';
pdsch3.DMRSPower = 0;
pdsch3.EnablePTRS = false;
pdsch3.PTRSPower = 0;

pdschreserved = nrPDSCHReservedConfig;
pdschreserved.PRBSet = 0:2;
pdschreserved.SymbolSet = [0 1];
pdschreserved.Period = 1;

pdsch3.ReservedPRB = {pdschreserved};

pdschdmrs = nrPDSCHDMRSConfig;
pdschdmrs.DMRSConfigurationType = 1;
pdschdmrs.DMRSReferencePoint = 'CRB0';
pdschdmrs.DMRSTypeAPosition = 2;
pdschdmrs.DMRSAdditionalPosition = 1;
pdschdmrs.DMRSLength = 1;
pdschdmrs.CustomSymbolSet = [];
pdschdmrs.DMRSPortSet = [];
pdschdmrs.NIDNSCID = [];
pdschdmrs.NSCID = 0;
pdschdmrs.NumCDMGroupsWithoutData = 1;
pdschdmrs.DMRSDownlinkR16 = false;

pdsch3.DMRS = pdschdmrs;

pdschptrs = nrPDSCHPTRSConfig;
pdschptrs.TimeDensity = 4;
pdschptrs.FrequencyDensity = 2;
pdschptrs.REOffset = '00';
pdschptrs.PTRSPortSet = [];

pdsch3.PTRS = pdschptrs;

% PDSCH 4
pdsch4 = nrWavegenPDSCHConfig;
pdsch4.Enable = true;
pdsch4.Label = 'Partial band PDSCH sequence with QPSK modulation scheme (target, RNTI = 2) (Partial downlink slots)';
pdsch4.Power = 0;
pdsch4.BandwidthPartID = 1;
pdsch4.Modulation = 'QPSK';
pdsch4.NumLayers = 1;
pdsch4.MappingType = 'A';
pdsch4.ReservedCORESET = [];
pdsch4.SymbolAllocation = [2 8];
pdsch4.SlotAllocation = 3;
pdsch4.Period = 5;
pdsch4.PRBSet = 0:2;
pdsch4.PRBSetType = 'VRB';
pdsch4.VRBToPRBInterleaving = false;
pdsch4.VRBBundleSize = 2;
pdsch4.NID = [];
pdsch4.RNTI = 2;
pdsch4.Coding = false;
pdsch4.TargetCodeRate = 0.4785;
pdsch4.TBScaling = 1;
pdsch4.XOverhead = 0;
pdsch4.RVSequence = 0;
pdsch4.DataSource = 'PN23';
pdsch4.DMRSPower = 0;
pdsch4.EnablePTRS = false;
pdsch4.PTRSPower = 0;

pdschreserved = nrPDSCHReservedConfig;
pdschreserved.PRBSet = 0:2;
pdschreserved.SymbolSet = [0 1];
pdschreserved.Period = 1;

pdsch4.ReservedPRB = {pdschreserved};

pdschdmrs = nrPDSCHDMRSConfig;
pdschdmrs.DMRSConfigurationType = 1;
pdschdmrs.DMRSReferencePoint = 'CRB0';
pdschdmrs.DMRSTypeAPosition = 2;
pdschdmrs.DMRSAdditionalPosition = 1;
pdschdmrs.DMRSLength = 1;
pdschdmrs.CustomSymbolSet = [];
pdschdmrs.DMRSPortSet = [];
pdschdmrs.NIDNSCID = [];
pdschdmrs.NSCID = 0;
pdschdmrs.NumCDMGroupsWithoutData = 1;
pdschdmrs.DMRSDownlinkR16 = false;

pdsch4.DMRS = pdschdmrs;

pdschptrs = nrPDSCHPTRSConfig;
pdschptrs.TimeDensity = 4;
pdschptrs.FrequencyDensity = 2;
pdschptrs.REOffset = '00';
pdschptrs.PTRSPortSet = [];

pdsch4.PTRS = pdschptrs;

cfgDLTM.PDSCH = {pdsch1,pdsch2,pdsch3,pdsch4};

%% CSI-RS Instances Configuration
csirs = nrWavegenCSIRSConfig;
csirs.Enable = false;
csirs.Label = 'CSIRS1';
csirs.Power = 0;
csirs.BandwidthPartID = 1;
csirs.CSIRSType = 'nzp';
csirs.CSIRSPeriod = 'on';
csirs.RowNumber = 1;
csirs.Density = 'three';
csirs.SymbolLocations = 0;
csirs.SubcarrierLocations = 0;
csirs.NumRB = 216;
csirs.RBOffset = 0;
csirs.NID = 1;

cfgDLTM.CSIRS = {csirs};

% Generation
[waveform,info] = nrWaveformGenerator(cfgDLTM);

Fs = info.ResourceGrids(1).Info.SampleRate; 								 % Specify the sample rate of the waveform in Hz

%% Visualize
% Time Scope
timeScope = timescope('SampleRate', Fs, ...
    'TimeSpanOverrunAction', 'scroll', ...
    'TimeSpanSource', 'property', ...
    'TimeSpan', 4.8828e-07);
timeScope(waveform);
release(timeScope);


