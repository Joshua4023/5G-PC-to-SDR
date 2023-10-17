% Generated by MATLAB(R) 23.2 (R2023b) and 5G Toolbox 23.2 (R2023b).
% Generated on: 11-Oct-2023 19:37:03

%% Generating Uplink FRC waveform
% Uplink FRC configuration
cfgULFRC = nrULCarrierConfig;
cfgULFRC.Label = 'G-FR1-A2-1';
cfgULFRC.FrequencyRange = 'FR1';
cfgULFRC.ChannelBandwidth = 40;
cfgULFRC.NCellID = 1;
cfgULFRC.NumSubframes = 10;
cfgULFRC.InitialNSubframe = 0;
cfgULFRC.WindowingPercent = 0;
cfgULFRC.SampleRate = [];
cfgULFRC.CarrierFrequency = 0;

%% SCS specific carriers
scscarrier = nrSCSCarrierConfig;
scscarrier.SubcarrierSpacing = 15;
scscarrier.NSizeGrid = 216;
scscarrier.NStartGrid = 0;

cfgULFRC.SCSCarriers = {scscarrier};

%% Bandwidth Parts
bwp = nrWavegenBWPConfig;
bwp.BandwidthPartID = 1;
bwp.Label = 'BWP1';
bwp.SubcarrierSpacing = 15;
bwp.CyclicPrefix = 'normal';
bwp.NSizeBWP = 216;
bwp.NStartBWP = 0;

cfgULFRC.BandwidthParts = {bwp};

intracellguardbands = nrIntraCellGuardBandsConfig;
intracellguardbands.GuardBandSize = zeros(0,2);
intracellguardbands.SubcarrierSpacing = 15;

cfgULFRC.IntraCellGuardBands = {intracellguardbands};

%% PUSCH Instances Configuration
pusch = nrWavegenPUSCHConfig;
pusch.Enable = true;
pusch.Label = 'PUSCH sequence for G-FR1-A2-1';
pusch.Power = 0;
pusch.BandwidthPartID = 1;
pusch.Modulation = '16QAM';
pusch.NumLayers = 1;
pusch.MappingType = 'A';
pusch.SymbolAllocation = [0 14];
pusch.SlotAllocation = 0:9;
pusch.Period = 10;
pusch.PRBSet = 95:119;
pusch.TransformPrecoding = false;
pusch.TransmissionScheme = 'codebook';
pusch.NumAntennaPorts = 1;
pusch.TPMI = 0;
pusch.FrequencyHopping = 'neither';
pusch.SecondHopStartPRB = 0;
pusch.Interlacing = false;
pusch.RBSetIndex = 0;
pusch.InterlaceIndex = 0;
pusch.NID = [];
pusch.RNTI = 1;
pusch.NRAPID = [];
pusch.Coding = true;
pusch.TargetCodeRate = 0.642578125;
pusch.XOverhead = 0;
pusch.RVSequence = 0;
pusch.DataSource = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
pusch.EnableACK = false;
pusch.NumACKBits = 10;
pusch.BetaOffsetACK = 20;
pusch.DataSourceACK = 'PN9-ITU';
pusch.EnableCSI1 = false;
pusch.NumCSI1Bits = 10;
pusch.BetaOffsetCSI1 = 6.25;
pusch.DataSourceCSI1 = 'PN9-ITU';
pusch.EnableCSI2 = false;
pusch.NumCSI2Bits = 10;
pusch.BetaOffsetCSI2 = 6.25;
pusch.DataSourceCSI2 = 'PN9-ITU';
pusch.EnableCGUCI = false;
pusch.NumCGUCIBits = 7;
pusch.BetaOffsetCGUCI = 20;
pusch.DataSourceCGUCI = 'PN9-ITU';
pusch.EnableULSCH = true;
pusch.UCIScaling = 1;
pusch.DMRSPower = 3;
pusch.EnablePTRS = false;
pusch.PTRSPower = 0;

% PUSCH DM-RS
puschDMRS = nrPUSCHDMRSConfig;
puschDMRS.DMRSConfigurationType = 1;
puschDMRS.DMRSTypeAPosition = 2;
puschDMRS.DMRSAdditionalPosition = 1;
puschDMRS.DMRSLength = 1;
puschDMRS.CustomSymbolSet = [];
puschDMRS.DMRSPortSet = 0;
puschDMRS.NIDNSCID = [];
puschDMRS.NSCID = 0;
puschDMRS.GroupHopping = false;
puschDMRS.SequenceHopping = false;
puschDMRS.NRSID = [];
puschDMRS.NumCDMGroupsWithoutData = 2;
puschDMRS.DMRSUplinkR16 = false;
puschDMRS.DMRSUplinkTransformPrecodingR16 = false;

pusch.DMRS = puschDMRS;

% PUSCH PT-RS
puschPTRS = nrPUSCHPTRSConfig;
puschPTRS.TimeDensity = 1;
puschPTRS.FrequencyDensity = 2;
puschPTRS.NumPTRSSamples = 2;
puschPTRS.NumPTRSGroups = 2;
puschPTRS.REOffset = '00';
puschPTRS.PTRSPortSet = 0;
puschPTRS.NID = [];

pusch.PTRS = puschPTRS;

cfgULFRC.PUSCH = {pusch};

%% PUCCH Instances Configuration
pucch = nrWavegenPUCCH0Config;
pucch.Enable = false;
pucch.Label = 'PUCCH format 0';
pucch.Power = 0;
pucch.BandwidthPartID = 1;
pucch.SymbolAllocation = [13 1];
pucch.SlotAllocation = 0:9;
pucch.Period = 10;
pucch.PRBSet = 0;
pucch.FrequencyHopping = 'neither';
pucch.SecondHopStartPRB = 1;
pucch.Interlacing = false;
pucch.RBSetIndex = 0;
pucch.InterlaceIndex = 0;
pucch.GroupHopping = 'neither';
pucch.HoppingID = [];
pucch.InitialCyclicShift = 0;
pucch.NumUCIBits = 1;
pucch.DataSourceUCI = 'PN9-ITU';
pucch.DataSourceSR = 0;

cfgULFRC.PUCCH = {pucch};

%% SRS Instances Configuration
srs = nrWavegenSRSConfig;
srs.Enable = false;
srs.Label = 'SRS1';
srs.Power = 0;
srs.BandwidthPartID = 1;
srs.NumSRSPorts = 1;
srs.SymbolStart = 13;
srs.NumSRSSymbols = 1;
srs.SlotAllocation = 0:9;
srs.Period = 10;
srs.FrequencyStart = 0;
srs.NRRC = 0;
srs.CSRS = 0;
srs.BSRS = 0;
srs.BHop = 0;
srs.Repetition = 1;
srs.KTC = 2;
srs.KBarTC = 0;
srs.FrequencyScalingFactor = 1;
srs.StartRBIndex = 0;
srs.EnableStartRBHopping = false;
srs.CyclicShift = 0;
srs.GroupSeqHopping = 'neither';
srs.NSRSID = 0;
srs.SRSPositioning = false;

cfgULFRC.SRS = {srs};

% Generation
[waveform,info] = nrWaveformGenerator(cfgULFRC);

Fs = info.ResourceGrids(1).Info.SampleRate; 								 % Specify the sample rate of the waveform in Hz

%% Visualize
% Time Scope
timeScope = timescope('SampleRate', Fs, ...
    'TimeSpanOverrunAction', 'scroll', ...
    'TimeSpanSource', 'property', ...
    'TimeSpan', 4.8828e-07);
timeScope(waveform);
release(timeScope);

% Spectrum Analyzer
spectrum = spectrumAnalyzer('SampleRate', Fs);
spectrum(waveform);
release(spectrum);




