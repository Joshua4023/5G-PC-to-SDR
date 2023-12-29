% hNRReferenceWaveformGenerator 5G NR reference waveform generation
%   hNRReferenceWaveformGenerator is a class which enables the creation of 
%   unencrypted, golden reference 5G NR baseband waveforms for 
%   simulation, verification and test and measurement applications. The
%   generator supports the NR test model (NR-TM) waveforms, defined in 
%   TS 38.141-1 Section 4.9.2 (FR1) and TS 38.141-2 Section 4.9.2 (FR2).
%   The generator also supports the DL fixed reference channels (FRC), 
%   defined in TS 38.101-1 Annex A.3 (FR1) and TS 38.101-2 Annex A.3 (FR2).
%   For the uplink, the generator supports the UL FRC defined in TS 38.104
%   Annex A. By default the generated waveform length is 10ms for FDD and 
%   20ms for TDD.
% 
%   WG = hNRReferenceWaveformGenerator() creates a default object, WG, 
%   configured for the NR-FR1-TM1.1 NR test model, 10MHz bandwidth,
%   15kHz subcarrier spacing, FDD duplexing, and NCellID = 1. The test 
%   model specification used is taken from TS 38.141-1 v17.8.0.
% 
%   WG = hNRReferenceWaveformGenerator(RC,BW,SCS,DM,NCELLID,SV,CS,OCNG) creates
%   a generator object, WG, given the reference waveform identifier RC, the 
%   channel bandwidth BW, and the subcarrier spacing SCS. All parameters 
%   are optional. For the NR-TM and downlink FRC, the default bandwidth and
%   subcarrier spacing values are 10MHz and 15kHz respectively for FR1, and
%   100MHz and 120kHz for FR2. For the uplink FRC, the default bandwidth 
%   and subcarrier spacing values depends on the individual FRC, as listed
%   in TS 38.104 Annex A. These values can be overridden when they are 
%   supplied as inputs. Note that, in all cases, the combination of BW and
%   SCS must be a valid combination from the associated bandwidth 
%   configuration table (FR1 or FR2). The duplex mode defaults to FDD for
%   FR1 and TDD for FR2 (FR2 NR-TMs are only formally defined for TDD). 
%   The cell identity, NCELLID, defaults to 1 for the downlink and 0 for 
%   the uplink. This is also used to set the various NCellID controllable
%   scrambling identities, where appropriate. If any of the BW, SCS, DM and
%   NCELLID inputs are supplied but empty then they take their default
%   values. The additional optional inputs, SV, CS and OCNG, are described below.
%  
%   The RC identifier should be a char vector or string from the set of 
%   FR1 NR-TM  ('NR-FR1-TM1.1','NR-FR1-TM1.2','NR-FR1-TM2','NR-FR1-TM2a',
%               'NR-FR1-TM2b','NR-FR1-TM3.1','NR-FR1-TM3.1a','NR-FR1-TM3.1b',
%               'NR-FR1-TM3.2','NR-FR1-TM3.3')
%   FR2 NR-TM  ('NR-FR2-TM1.1','NR-FR2-TM2','NR-FR2-TM2a','NR-FR2-TM3.1','NR-FR2-TM3.1a')
%   FR1 DL FRC ('DL-FRC-FR1-QPSK','DL-FRC-FR1-64QAM','DL-FRC-FR1-256QAM','DL-FRC-FR1-1024QAM')
%   FR2 DL FRC ('DL-FRC-FR2-QPSK','DL-FRC-FR2-16QAM','DL-FRC-FR2-64QAM')
%   FR1 UL FRC ('G-FR1-Ax-y', where x and y values come from TS 38.104 Annex A)
%   FR2 UL FRC ('G-FR2-Ax-y', where x and y values come from TS 38.104 Annex A)
%   
%   The BW input is the channel bandwidth as a char vector or string or 
%   numerical input, for example '10MHz', "10MHz" or 10. For backward
%   compatibility, it is also possible to prefix the text inputs with 'BW_'.
%   If it is empty then the default value is selected.
% 
%   The SCS input is the subcarrier spacing as a char vector or string or
%   numerical input, for example '15kHz', "15kHz" or 15. If it is empty 
%   then the default value is selected.
% 
%   The DM input is the duplexing mode and should be a char vector or
%   string and either 'FDD' or 'TDD'. If it is empty then the default value
%   is selected.
% 
%   The NCELLID input is the physical-layer cell identity. For NR-TM, the
%   3GPP standard defines NCELLID = 1 for the lowest configured carrier, 
%   NCELLID = 2 for 2nd lowest configured carrier, etc. Early version of
%   TS 38.141-1 and TS 38.141-2 specified NCELLID = 0 however this exhibits
%   high PAPR. The NCELLID value is also used to set the various scrambling 
%   identities that are allowed to be set to the cell identity in the 
%   NR specification.
%    
%   The SV input is applicable to NR-TM only and specifies the version
%   of TS 38.141-1 and TS 38.141-2 used ('15.1.0','15.2.0','15.7.0','16.7.0','17.8.0').
%   Note that V15.0.0 of the NR-TM is the same as the V15.1.0 specification
%   but with NCellID = 0 which can lead to a high PAPR. The logical CS input, 
%   if set to true, enables the use of a legacy MATLAB structure to contain
%   the parameters within the Config property. By default, the parameters
%   are represented by a nrDLCarrierConfig or nrULCarrierConfig object.
%
%   The OCNG input only applies to DL FRCs. The standard TS 38.101 Annex
%   A.5 defines OCNG patterns for downlink reference measurement channels
%   (DL RMCs) as a way of ensuring uniform power distribution in those
%   waveforms. If set to true, all the unused REs are filled in with new
%   PDSCHs. If set to false (default), no new PDSCHs are generated.
%
%   hNRReferenceWaveformGenerator properties:
%
%   Bandwidth configuration tables and reference model lists (Read-only)
%   FR1BandwidthTable - FR1 transmission bandwidth configurations (TS 38.104 Table 5.3.2-1)
%   FR2BandwidthTable - FR2 transmission bandwidth configurations (TS 38.104 Tables 5.3.2-2 and 5.3.2-3)
%   FR1TestModels     - FR1 NR test model (NR-TM) list (TS 38.141-1 Section 4.9.2)
%   FR2TestModels     - FR2 NR test model (NR-TM) list (TS 38.141-2 Section 4.9.2)
%   FR1DownlinkFRC    - FR1 downlink fixed reference channels (TS 38.101-1 Annex A.3)
%   FR2DownlinkFRC    - FR2 downlink fixed reference channels (TS 38.101-2 Annex A.3)
%   FR1UplinkFRC      - FR1 uplink fixed reference channels (TS 38.104 Annex A)
%   FR2UPlinkFRC      - FR2 uplink fixed reference channels (TS 38.104 Annex A)
% 
%   Configuration parameters (Read-only)
%   ConfiguredModel - Configured reference model parameters
%   TargetRNTI      - RNTI values of the target reference measurement channels
%   Config          - Generator parameter object (nrDLCarrierConfig) or structure. Read-only by default but can be made writable using makeConfigWritable function
%   IsReadOnly      - Indication of whether the Config parameter property is read-only or writable
%
%   hNRReferenceWaveformGenerator methods:
%
%   hNRReferenceWaveformGenerator - Class constructor
%   generateWaveform              - Generate baseband waveform, resource grids and info structure
%   createWaveformSource        -   Create a waveform source to stream the waveform per sample
%   displayResourceGrid           - Display image of resource grid across all subframes
%   makeConfigWritable            - Make the Config parameter object or structure writable
%
%   Examples:
%   % Display the standard-defined bandwidth configuration tables for 
%   % FR1 and FR2.
% 
%   hNRReferenceWaveformGenerator.FR1BandwidthTable
%   hNRReferenceWaveformGenerator.FR2BandwidthTable  
% 
%   % Display the Release 16 NR-TMs, defined for FR1.
%   
%   hNRReferenceWaveformGenerator.FR1TestModels
% 
%   % Create a waveform generator for FDD NR-FR1-TM1.2, visualize the
%   % resource element layout and generate the associated baseband 
%   % waveform. This particular NR-TM is used for tests on base station
%   % unwanted emissions (ACLR and operating band unwanted emissions) and
%   % has bands of boosted and deboosted PRB which can be seen in the
%   % subcarrier magnitude plot.
%
%   wavegen = hNRReferenceWaveformGenerator('NR-FR1-TM1.2','40MHz','15kHz','FDD');
%   displayResourceGrid(wavegen);
%   [waveform,waveinfo] = generateWaveform(wavegen);
%
%   % The sample rate of the waveform is part of the optional info 
%   % output structure, which also contains the carrier/BWP resource
%   % element grid matrices
%
%   waveinfo.Info.SampleRate
% 
%   % The full parameter set which configures the reference model is 
%   % contained in the Config property. This is read-only, by default, to
%   % prevent inadvertent modification. Use the makeConfigWritable function 
%   % to allow these parameters to be modified. For example, to turn the 
%   % DL-SCH transport coding on in the test models, set the 'Coding'
%   % parameters to true.
%   
%   mywavegen = makeConfigWritable(wavegen);        % Make the Config property writable
%   pdscharray = [mywavegen.Config.PDSCH{:}];       % Extract all PDSCH configs into an array
%   [pdscharray.Coding] = deal(1);                  % Set all Coding parameters to turn transport coding on
%   mywavegen.Config.PDSCH = num2cell(pdscharray);  % Reassign the updated PDSCH configs
%   codedwaveform = generateWaveform(mywavegen);    % Generate the new waveform
%
%   See also nrDLCarrierConfig, nrULCarrierConfig, nrWaveformGenerator.

%   Copyright 2019-2023 The MathWorks, Inc.

classdef hNRReferenceWaveformGenerator
    
    % Static private properties
    properties (Constant, Access=private)
        DLRMConfigs = getRMConfigs();                               % Internal set of downlink reference model definitions
        DLModels = [hNRReferenceWaveformGenerator.DLRMConfigs.Name];% Internal list of downlink reference model names
        ULModels = getULModels();                                   % Internal list of uplink reference model names
    end
 
    % Static public properties 
    properties (Constant)
        FR1BandwidthTable = getFR1BandwidthTable(); % FR1 transmission bandwidth configurations (TS 38.104 Table 5.3.2-1)
        FR2BandwidthTable = getFR2BandwidthTable(); % FR2 transmission bandwidth configurations (TS 38.104 Tables 5.3.2-2 and 5.3.2-3)
        FR1TestModels = hNRReferenceWaveformGenerator.DLModels(contains(hNRReferenceWaveformGenerator.DLModels,'FR1-TM'))';    % FR1 test model names (TS 38.141-1 Section 4.9.2)
        FR2TestModels = hNRReferenceWaveformGenerator.DLModels(contains(hNRReferenceWaveformGenerator.DLModels,'FR2-TM'))';    % FR2 test model names (TS 38.141-2 Section 4.9.2)
        FR1DownlinkFRC = hNRReferenceWaveformGenerator.DLModels(contains(hNRReferenceWaveformGenerator.DLModels,'FRC-FR1'))';  % FR1 DL FRC names (TS 38.101-1 Annex A.3)
        FR2DownlinkFRC = hNRReferenceWaveformGenerator.DLModels(contains(hNRReferenceWaveformGenerator.DLModels,'FRC-FR2'))';  % FR2 DL FRC names (TS 38.101-2 Annex A.3)
        FR1UplinkFRC = hNRReferenceWaveformGenerator.ULModels(contains(hNRReferenceWaveformGenerator.ULModels,'G-FR1'))';      % FR1 UL FRC names (TS 38.104 Annex A)
        FR2UplinkFRC = hNRReferenceWaveformGenerator.ULModels(contains(hNRReferenceWaveformGenerator.ULModels,'G-FR2'))';      % FR2 UL FRC names (TS 38.104 Annex A)
    end

    properties (Dependent,Access=public)
        Config;             % Full parameter set for the reference waveform generation. Read-only unless made writable using makeConfigWritable       
    end
               
    properties (SetAccess=private)
        IsReadOnly = 1;     % Read only. Indicator of whether Config parameters are read only or modifiable
        ConfiguredModel;    % Read only. Top-level parameters defining the reference waveform
        TargetRNTI;         % Read only. RNTI values of the target reference measurement channels
    end
    
    properties (Access=private)
        LinkDirection;      % Uplink or downlink
    end
    
    % Methods associated with Config property runtime access
    methods
        function obj = makeConfigWritable(obj)
        %  makeConfigWritable Make the Config parameter object or structure writable
        %   WG = makeConfigWritable(WG) allows the Config parameter object or structure
        %   to be modified. By default, the Config property is read-only since it 
        %   specifies the configuration of the standard-defined reference waveform. 
        %   This function makes it writable so that the waveform can be further customized. 
        
            obj.IsReadOnly = 0; 
        end
        
        % Config property setter
        function obj = set.Config(obj,c)
            if obj.IsReadOnly
                error('The waveform configuration is set to read only. Use makeConfigWritable to make the Config property writable.');
            end
            obj.ConfigValue = c;
        end
        
        % Config property getter
        function c = get.Config(obj)
            c = obj.ConfigValue;
        end
           
    end
   
    properties (Access=private)
         ConfigValue;          % Internal 'backing store' of the public Config parameter object/structure property
    end
    
    % Methods associated with main waveform API
    methods        
        
        % Class constructor
        function obj = hNRReferenceWaveformGenerator(rc,bw,scs,duplexmode,ncellid,sv,cs,ocng)
        % hNRReferenceWaveformGenerator 5G NR reference waveform generator
        %   WG = hNRReferenceWaveformGenerator() creates a default object, WG, 
        %   configured for the NR-FR1-TM1.1 NR test model, 10MHz bandwidth and
        %   15kHz subcarrier spacing. The generated waveform length is
        %   10ms for FDD and 20ms for TDD. The version of TS 38.141-1 and
        %   TS 38.141-2 used is v17.8.0.
        %
        %   WG = hNRReferenceWaveformGenerator(RC,BW,SCS,DM,NCELLID,SV,CS,OCNG) creates
        %   a generator object, WG, given the reference waveform identifier RC,
        %   the channel bandwidth BW, and the subcarrier spacing SCS.   
        %   The RC identifier should be a char vector or string from the set of 
        %   FR1 NR-TM  ('NR-FR1-TM1.1','NR-FR1-TM1.2','NR-FR1-TM2','NR-FR1-TM2a',
        %               'NR-FR1-TM2b','NR-FR1-TM3.1','NR-FR1-TM3.1a','NR-FR1-TM3.1b',
        %               'NR-FR1-TM3.2','NR-FR1-TM3.3')
        %   FR2 NR-TM  ('NR-FR2-TM1.1','NR-FR2-TM2','NR-FR2-TM2a','NR-FR2-TM3.1','NR-FR2-TM3.1a')
        %   FR1 DL FRC ('DL-FRC-FR1-QPSK','DL-FRC-FR1-64QAM','DL-FRC-FR1-256QAM','DL-FRC-FR1-1024QAM')
        %   FR2 DL FRC ('DL-FRC-FR2-QPSK','DL-FRC-FR2-16QAM','DL-FRC-FR2-64QAM')
        %   The BW input is the channel bandwidth as a char vector or string or 
        %   numerical input, for example '10MHz', "10MHz" or 10. It is also 
        %   possible to prefix the text inputs with 'BW_' to be compatible with 
        %   the column names in the bandwidth configuration table properties.
        %   The SCS input is the subcarrier spacing as a char vector or string or
        %   numerical input, for example '15kHz', "15kHz" or 15.
        %   The DM input is the duplexing mode and should be a char vector or
        %   string and either 'FDD' or 'TDD'.
        %   The NCELLID input is the cell identity and should be 1,2,...N
        %   starting at 1 for the lowest configured carrier.
        %   The SV input is applicable to NR-TM only and specifies the version
        %   of TS 38.141-1 and TS 38.141-2 required ('15.1.0', '15.2.0',
        %   '15.7.0','16.7.0' or '17.8.0').
        %   The CS input only applies to downlink. If set to true, the 
        %   Config property is a structure. If set to false (default), the
        %   Config property is an nrDLCarrierConfig object. This input 
        %   provides limited backwards compatibility with the structure-based
        %   parameterization of previous 5G Toolbox examples.
        %   The OCNG input only applies to DL FRCs. If set to true, all the
        %   unused REs are filled in with new PDSCHs. If set to false
        %   (default), no new PDSCHs are generated.
        %   
        %   Note that the combination of BW and SCS must be a valid combination 
        %   from the associated bandwidth configuration table (FR1 or FR2).    
            
            if nargin < 8           % OCNG not defined
                ocng = 0;   
                if nargin < 7           % Use of structure config not defined
                    cs = [];     
                    if nargin < 6           % Version not defined
                        sv = "17.8.0";
                        if nargin < 5           % NCellID not defined
                            ncellid = [];
                            if nargin < 4           % Duplexing not defined
                                duplexmode = [];
                                if nargin < 3           % SCS not defined
                                   scs = []; % "15kHz";
                                   if nargin < 2            % BW not defined
                                      bw = []; % "10MHz";
                                      if nargin < 1             % Nothing defined
                                         rc = hNRReferenceWaveformGenerator.DLModels(1);
                                      end
                                   end
                                end
                            end
                        end
                    end
                end
            end
                   
            % Set up some variables on the link direction
            if startsWith(rc,'G-','IgnoreCase',true)
                obj.LinkDirection = "uplink";
                refnames = hNRReferenceWaveformGenerator.ULModels;
            else
                obj.LinkDirection = "downlink";
                refnames = hNRReferenceWaveformGenerator.DLModels;                
                % Downlink specific defaulting choices
                frselect = 1+contains(rc,'FR2','IgnoreCase',true);
                if isempty(scs)
                    scsdefs = ["15kHz","120kHz"];   % SCS defaults for FR1,FR2
                    scs = scsdefs(frselect);
                end
                if isempty(bw)
                    bwdefs = ["10MHz","100MHz"];    % BW defaults for FR1,FR2
                    bw = bwdefs(frselect);
                end
            end
            
            % By default, do not use the structure representation for the parameter configuration
            if isempty(cs)
                cs = 0;
            end
                
            % By default, use NCellID to control the various scrambling identities
            % using NCellID=1 for the downlink and 0 for the uplink
            if isempty(ncellid)
               ncellid = double(obj.LinkDirection == "downlink");
            end
            
            % Validate the reference name against the uplink/downlink name set
            selected = strcmpi(rc,refnames);
            if ~any(selected)
                error('The %s reference model name (%s) must be one of the set (%s).',obj.LinkDirection,rc,join(refnames,', '));
            end
                 
            % Get internal reference channel definition for the selected model
            if obj.LinkDirection == "uplink"
                tdef = getULFRCDefinition(rc,bw,scs);
            else
                tdef = hNRReferenceWaveformGenerator.DLRMConfigs(selected);
            end
                       
            % Adjust for TDD overridding 
            % Check the input duplex mode, if presented, then use it to override
            % the value from the internal definition
            if ~isempty(duplexmode)
                duplexingmodes = ["FDD","TDD"];
                selected = strcmpi(string(duplexmode),duplexingmodes);
                if ~any(selected)
                    error('The duplexing mode (%s) must be one of the set (%s).',string(duplexmode),join(duplexingmodes,', '));
                end
                tdef.DuplexMode = duplexingmodes(selected);
            end
            
            % Turn the reference channel definition into a full parameter set
            if obj.LinkDirection == "uplink"
                % Create the parameter set from the UL FRC reference definition
                obj.ConfigValue = getFRCPUSCHParameters(tdef,ncellid);
                obj.ConfiguredModel = {rc,bw,scs,tdef.DuplexMode,ncellid};  % Store the top-level configuration parameters
            else    
                % Build the full configuration parameter structure
                if contains(rc,'FRC','IgnoreCase',true) 
                    % Create the parameter set from the DL FRC reference definition
                    waveconfig = getCommonParameters(tdef,bw,scs,ncellid);      % Get the overall common downlink parameter set, excluding the PDSCH        
                    waveconfig = getFRCPDSCHParameters(tdef,waveconfig);        % Add in the PDSCH part into the parameter set

                    % Update the CORESET sequence so that it aligns with the active PDSCH FRC instances
                    waveconfig.CORESET.AllocatedSlots = waveconfig.PDSCH.AllocatedSlots;
                    waveconfig.CORESET.AllocatedPeriod = waveconfig.PDSCH.AllocatedPeriod;

                    obj.ConfigValue = waveconfig;
                    obj.ConfiguredModel = {rc,bw,scs,tdef.DuplexMode,ncellid};  % Store the top-level configuration parameters
                else
                    % Create the parameter set from the NR-TM reference definition
                    waveconfig = getCommonParameters(tdef,bw,scs,ncellid);              % Get the overall common parameter set, excluding the PDSCH
                    obj.ConfigValue = getTestModelPDSCHParameters(tdef,waveconfig,sv);  % Add in the PDSCH part and store parameters
                    obj.ConfiguredModel = {rc,bw,scs,tdef.DuplexMode,ncellid,sv};       % Store the top-level configuration parameters
                end              
            end

            % If not using parameter structures then convert preset definition to an object
            if ~cs
                obj.ConfigValue = getGenConfigObj(obj);
                % OCNG capability in DL FRCs is only available when using objects
                if contains(rc,'FRC','IgnoreCase',true) && obj.LinkDirection == "downlink" && ocng==1
                    obj.ConfigValue.PDSCH = [obj.ConfigValue.PDSCH nr5g.internal.wavegen.getOCNGPDSCHs(obj.ConfigValue)];
                end
            end  
                      
            % Configure the RNTI of the channels intended for reference measurement
            if contains(rc,["TM3.1", "TM3.1a", "TM3.1b"],'IgnoreCase',true)        % "NR-FR1-TM3.1", "NR-FR1-TM3.1a", "NR-FR1-TM3.1b", "NR-FR2-TM3.1", "NR-FR2-TM3.1a"
                targets = [0 2];    
            elseif contains(rc,["TM3.2", "TM3.3"],'IgnoreCase',true)               % "NR-FR1-TM3.2", "NR-FR1-TM3.3"
                targets = 1;
            elseif contains(rc,["TM1.2", "TM2", "TM2a", "TM2b"],'IgnoreCase',true) % "NR-FR1-TM1.2", "NR-FR1-TM2", "NR-FR1-TM2a", "NR-FR1-TM2b", "NR-FR2-TM2", "NR-FR2-TM2a"
                targets = 2;
            else
                targets = 0; 
            end 
            obj.TargetRNTI = targets;
        end
        
        function varargout = generateWaveform(obj,numsf)
        % generateWaveform Generate baseband waveform, resource grids and info structure
        %  [WAVEFORM,INFO] = generateWaveform(WG) generates the reference baseband 
        %  waveform WAVEFORM and information structures INFO, for the default 
        %  lengths of 10ms for FDD and 20ms for TDD.
        %  
        %  [WAVEFORM,INFO] = generateWaveform(WG,NUMSF) is the same as above but
        %  the generated waveform is of NUMSF subframes in length.
            
            nsf = obj.ConfigValue.NumSubframes;
            if nargin > 1
                obj.ConfigValue.NumSubframes = numsf;
            end
                
            cv = getGenConfigObj(obj);
            [wave,winfo] = nrWaveformGenerator(cv);

            % Add a name label into the waveform resources structure
            winfo.WaveformResources.Label = obj.ConfiguredModel{1};

            % Replicate the legacy 'SamplingRate', 'SamplesPerSubframe' etc field names for local backwards compatibility
            rg = winfo.ResourceGrids;
            for n=1:length(rg)
                rg(n).Info.NSubcarriers = size(rg(n).ResourceGridInCarrier,1);
                rg(n).Info.SymbolsPerSubframe = rg(n).Info.SymbolsPerSlot * rg(n).Info.SlotsPerSubframe;
                rg(n).Info.SamplesPerSubframe = rg(n).Info.SampleRate/1000;
                rg(n).Info.SamplingRate = rg(n).Info.SampleRate;
            end
               
            % Bundle into the output
            varargout = {wave,rg,winfo};
          
            obj.ConfigValue.NumSubframes = nsf;
        end

        function [source,info] = createWaveformSource(obj,varargin)
        % createWaveformSource Create a baseband waveform source object
        %  [SRC,INFO] = createWaveformSource(WG) creates a waveform source
        %  object SRC and information structures INFO, for the default
        %  lengths of 10ms for FDD and 20ms for TDD. The source is a dsp.SourceSource
        %  object and waveform is output one sample at a time. It is cyclically
        %  extended beyond the default length.
        %  
        %  [SRC,INFO] = createWaveformSource(WG,NUMSF) is the same as above but
        %  the repeated waveform vector is of NUMSF subframes in length.
        
            [waveform,info] = generateWaveform(obj,varargin{:});
            source = dsp.SignalSource(waveform,SignalEndAction='cyclic repetition');
        end

        function  displayResourceGrid(obj,numsf)
        % displayResourceGrid Display plots of the underlying resource grid 
        %   displayResourceGrid(WG) displays the resource grid for the
        %   default lengths of 10ms for FDD and 20ms for TDD.
        %
        %   displayResourceGrid(WG,NUMSF) displays the resource grid for a 
        %   period of NUMSF subframes.
                
            nsf = obj.ConfigValue.NumSubframes;
            if nargin > 1
                obj.ConfigValue.NumSubframes = numsf;
            end

            cv = getGenConfigObj(obj);

            % Build a structure of the parameters needed for the plotting operations      
            sv.Name = obj.ConfiguredModel{1};
            sv.ChannelBandwidth = cv.ChannelBandwidth;
            sv.FrequencyRange = cv.FrequencyRange;

            % Map the SCS carrier definitions into structure form
            sv.Carriers(length(cv.SCSCarriers)) = struct();
            for n=1:length(sv.Carriers)
                sv.Carriers(n).SubcarrierSpacing = cv.SCSCarriers{n}.SubcarrierSpacing;
                sv.Carriers(n).NRB = cv.SCSCarriers{n}.NSizeGrid;
                sv.Carriers(n).RBStart = cv.SCSCarriers{n}.NStartGrid;
            end
            % Get list of all the defined SCS carrier spacings
            carrierscs = [sv.Carriers.SubcarrierSpacing];

            %  Map the BWP definitions into structure form
            sv.BWP(length(cv.BandwidthParts)) = struct();
            for n=1:length(sv.BWP)
                bwpscs = cv.BandwidthParts{n}.SubcarrierSpacing;
                sv.BWP(n).SubcarrierSpacing = bwpscs;                   
                sv.BWP(n).RBOffset = cv.BandwidthParts{n}.NStartBWP - sv.Carriers( find(bwpscs==carrierscs,1) ).RBStart;  % Need to find carrier associated with BWP                        
            end

            % Remove any of the SCS carriers that are not referenced by the BWP
            sv.Carriers(~ismember(carrierscs, [sv.BWP.SubcarrierSpacing])) = [];

            % Create a structure containing PRB grids in ResourceGridPRB field 
            figure;
            gridPRB = wirelessWaveformGenerator.internal.computeResourceGridPRB(cv);

            plotResourceGrid(gca, sv, cv, gridPRB, obj.LinkDirection == "downlink");

            % Create a new figure to display the SCS carrier alignment plot
            figure;
            wirelessWaveformGenerator.internal.plotCarriers(gca, cv);

            % Get the underlying resource grids required for the plotting
            [~,winfo] = nrWaveformGenerator(cv);
            rg = winfo.ResourceGrids;
            
            % Create a new figure to display the subcarrier grid plots
            figure;
            cmap = parula(64);
            for bp = 1:length(rg)
                % Plot the resource element grid (scaled complex magnitude)
                subplot(length(rg),1,bp)
                im = image(40*abs(rg(bp).ResourceGridInCarrier(:,:,1))); axis xy;     
                colormap(im.Parent,cmap);
                title(sprintf('BWP %d in Carrier (SCS=%dkHz)',bp,sv.BWP(bp).SubcarrierSpacing)); xlabel('Symbols'); ylabel('Subcarriers');
            end

            obj.ConfigValue.NumSubframes = nsf;

            % Augment the last carrier plot title
            h=get(gca,'Title');
            titletext = get(h,'String');
            set(h,'String',sprintf('%s: %s',obj.ConfiguredModel{1},titletext)); 
        end
    
    end
end

%% File local functions

% Get transmission bandwidth configuration table for FR1
function table = getFR1BandwidthTable()
    % TS 38.104 Table 5.3.2-1:   Transmission bandwidth configuration NRB for FR1
    % TS 38.101-1 Table 5.3.2-1: Maximum transmission bandwidth configuration NRB (FR1)
    % NRB, for BW and SCS
    % BW MHz    5   10 15 20  25  30  35  40  45  50  60  70  80  90  100
    nrbtable = [25  52 79 106 133 160 188 216 242 270 NaN NaN NaN NaN NaN;     % 15 kHz
                11  24 38 51  65  78  92  106 119 133 162 189 217 245 273;     % 30 kHz
                NaN 11 18 24  31  38  44  51  58  65  79  93  107 121 135];    % 60 kHz

    % Package NRB array into a table
    table = array2table(nrbtable,"RowNames",["15kHz","30kHz","60kHz"],"VariableNames",["5MHz", "10MHz", "15MHz", "20MHz", "25MHz", "30MHz", "35MHz", "40MHz", "45MHz", "50MHz", "60MHz", "70MHz", "80MHz", "90MHz", "100MHz"]);
    table.Properties.Description = 'TS 38.104 Table 5.3.2-1: Transmission bandwidth configuration NRB for FR1';
end

% Get transmission bandwidth configuration table for FR2
function table = getFR2BandwidthTable()
    % TS 38.104 Table 5.3.2-2:   Transmission bandwidth configuration NRB for FR2-1
    % TS 38.104 Table 5.3.2-3:   Transmission bandwidth configuration NRB for FR2-2
    % TS 38.101-2 Table 5.3.2-1: Maximum transmission bandwidth configuration NRB (FR2)
    % NRB, for BW and SCS
    % BW MHz    50  100 200 400 800 1600 2000
    nrbtable = [66  132 264 NaN NaN NaN  NaN;    % 60 kHz
                32  66  132 264 NaN NaN  NaN;    % 120 kHz
                NaN NaN NaN 66  124 248  NaN;    % 480 kHz
                NaN NaN NaN 33  62  124  148];   % 960 kHz

    % Package NRB array into a table
    table = array2table(nrbtable,"RowNames",["60kHz","120kHz","480kHz","960kHz"],"VariableNames",["50MHz", "100MHz", "200MHz", "400MHz", "800MHz", "1600MHz", "2000MHz"]);
    table.Properties.Description = 'TS 38.104 Tables 5.3.2-2 and 5.3.2-3: Transmission bandwidth configuration NRB for FR2';
end

% Get supported uplink FRC waveform definition list
function names = getULModels()

    % Expand s  
    expa = @(s,e,n)reshape(s + e + string(n(:)),1,[]);

    nr = {getNRBUplinkFRC(1),getNRBUplinkFRC(2)};                           % Get the allocation sizes for all the supported FRCs, grouped by FR and Annex section
    lens = cellfun(@(x)cellfun('length',x),nr,'uniformoutput',false);       % From this, establish the number of FRC in each FR and Annex section
    lens = reshape([lens{:}],[],2);                                         % Turn into columns for each FR
    tn = size(lens,1);                                                      % Number of sections
    zn = arrayfun(@(f,x,y)expa(expa(expa("G-FR","",f),"-A",x),"-",1:y),...  % Turn this information into a set of FRC names
            [1 2].*ones(tn,1),((1:tn).*ones(2,1))',lens,'uniformoutput',0);
    names = [zn{:}];

end

% Get supported waveform definition list
function p = getRMConfigs()

    % TS 38.141-1 Section 4.9.2 NR FR1 test models
    % 'NR-FR1-TM1.1','NR-FR1-TM1.2','NR-FR1-TM2','NR-FR1-TM2a','NR-FR1-TM2b','NR-FR1-TM3.1','NR-FR1-TM3.1a','NR-FR1-TM3.1b','NR-FR1-TM3.2','NR-FR1-TM3.3' 
    ptm.Name = "NR-FR1-TM1.1";
    ptm.FR = "FR1";
    ptm.DuplexMode = "FDD";
    ptm.BoostedPercent = 100;
    ptm.BoostedPower = 0;
    ptm.Modulation = "QPSK";
    
    ptm(end+1).Name = "NR-FR1-TM1.2";
    ptm(end).FR = "FR1";
    ptm(end).DuplexMode = "FDD";
    ptm(end).BoostedPercent = 40;
    ptm(end).BoostedPower = 3;
    ptm(end).Modulation = ["QPSK" "QPSK"];
    
    ptm(end+1).Name = "NR-FR1-TM2";
    ptm(end).FR = "FR1";
    ptm(end).DuplexMode = "FDD";
    ptm(end).BoostedPercent = 1;      % Indicates single PRB case
    ptm(end).BoostedPower = 0;
    ptm(end).Modulation = "64QAM";
    
    ptm(end+1).Name = "NR-FR1-TM2a";
    ptm(end).FR = "FR1";
    ptm(end).DuplexMode = "FDD";
    ptm(end).BoostedPercent = 1;      % Indicates single PRB case
    ptm(end).BoostedPower = 0;
    ptm(end).Modulation = "256QAM";
    
    ptm(end+1).Name = "NR-FR1-TM2b";
    ptm(end).FR = "FR1";
    ptm(end).DuplexMode = "FDD";
    ptm(end).BoostedPercent = 1;      % Indicates single PRB case
    ptm(end).BoostedPower = 0;
    ptm(end).Modulation = "1024QAM";

    ptm(end+1).Name = "NR-FR1-TM3.1";
    ptm(end).FR = "FR1";
    ptm(end).DuplexMode = "FDD";
    ptm(end).BoostedPercent = 100;
    ptm(end).BoostedPower = 0;
    ptm(end).Modulation = "64QAM";
    
    ptm(end+1).Name = "NR-FR1-TM3.1a";
    ptm(end).FR = "FR1";
    ptm(end).DuplexMode = "FDD";
    ptm(end).BoostedPercent = 100;
    ptm(end).BoostedPower = 0;
    ptm(end).Modulation = "256QAM";
    
    ptm(end+1).Name = "NR-FR1-TM3.1b";
    ptm(end).FR = "FR1";
    ptm(end).DuplexMode = "FDD";
    ptm(end).BoostedPercent = 100;
    ptm(end).BoostedPower = 0;
    ptm(end).Modulation = "1024QAM";

    ptm(end+1).Name = "NR-FR1-TM3.2";
    ptm(end).FR = "FR1";
    ptm(end).DuplexMode = "FDD";
    ptm(end).BoostedPercent = 60;
    ptm(end).BoostedPower = -3;
    ptm(end).Modulation = ["16QAM" "QPSK"];   % Target and power balancing
    
    ptm(end+1).Name = "NR-FR1-TM3.3";
    ptm(end).FR = "FR1";
    ptm(end).DuplexMode = "FDD";
    ptm(end).BoostedPercent = 50;
    ptm(end).BoostedPower = -6;
    ptm(end).Modulation = ["QPSK" "16QAM"];   % Target and power balancing (changed to QPSK in later std versions/downstream code)
    
    % TS 38.141-2 Section 4.9.2 NR FR2 test models
    % 'NR-FR2-TM1.1','NR-FR2-TM2','NR-FR2-TM2a','NR-FR2-TM3.1','NR-FR2-TM3.1a'
    ptm(end+1).Name = "NR-FR2-TM1.1";
    ptm(end).FR = "FR2";
    ptm(end).DuplexMode = "TDD";
    ptm(end).BoostedPercent = 100;
    ptm(end).BoostedPower = 0;
    ptm(end).Modulation = "QPSK";
    
    ptm(end+1).Name = "NR-FR2-TM2";
    ptm(end).FR = "FR2";
    ptm(end).DuplexMode = "TDD";
    ptm(end).BoostedPercent = 1;       % Indicates single PRB case
    ptm(end).BoostedPower = 0;
    ptm(end).Modulation = "64QAM";
    
    ptm(end+1).Name = "NR-FR2-TM2a";
    ptm(end).FR = "FR2";
    ptm(end).DuplexMode = "TDD";
    ptm(end).BoostedPercent = 1;      % Indicates single PRB case
    ptm(end).BoostedPower = 0;
    ptm(end).Modulation = "256QAM";
   
    ptm(end+1).Name = "NR-FR2-TM3.1";
    ptm(end).FR = "FR2";
    ptm(end).DuplexMode = "TDD";
    ptm(end).BoostedPercent = 100;
    ptm(end).BoostedPower = 0;
    ptm(end).Modulation = "64QAM";

    ptm(end+1).Name = "NR-FR2-TM3.1a";
    ptm(end).FR = "FR2";
    ptm(end).DuplexMode = "TDD";
    ptm(end).BoostedPercent = 100;
    ptm(end).BoostedPower = 0;
    ptm(end).Modulation = "256QAM";
        
    % Common to all NR-TM
    % CORESET/PDCCH specification
    [ptm.ControlNRB] = deal(6);       % 6 PRB assigned to CORESET
    [ptm.NumCCE] = deal(1);           % Single CCE for the PDCCH
    [ptm.MCSIndex] = deal(-1);        % Transport/coding is off
    
    % TS 38.101-1 A.3 DL RMC (FRC) (TDD definitions)  (Both FDD & TDD definitions)
    pfrc.Name = "DL-FRC-FR1-QPSK";
    pfrc.FR = "FR1";
    pfrc.DuplexMode = "FDD";
    pfrc.BoostedPercent = 100;
    pfrc.BoostedPower = 0;
    pfrc.Modulation = "QPSK"; 
    pfrc.MCSIndex = 4;  % FRC TCR = 1/3,  MCS 4 in table 1 (64QAM)
    
    pfrc(end+1).Name = "DL-FRC-FR1-64QAM";
    pfrc(end).FR = "FR1";
    pfrc(end).DuplexMode = "FDD";
    pfrc(end).BoostedPercent = 100;
    pfrc(end).BoostedPower = 0;
    pfrc(end).Modulation = "64QAM";
    pfrc(end).MCSIndex = 24;  % FRC TCR = 3/4,  MCS 24 in table 1 (64QAM)
       
    pfrc(end+1).Name = "DL-FRC-FR1-256QAM";
    pfrc(end).FR = "FR1";
    pfrc(end).DuplexMode = "FDD";
    pfrc(end).BoostedPercent = 100;
    pfrc(end).BoostedPower = 0;
    pfrc(end).Modulation = "256QAM";
    pfrc(end).MCSIndex = 23;        % FRC TCR = 4/5,  MCS 23 in table 2 (256QAM)
     
    pfrc(end+1).Name = "DL-FRC-FR1-1024QAM";
    pfrc(end).FR = "FR1";
    pfrc(end).DuplexMode = "FDD";
    pfrc(end).BoostedPercent = 100;
    pfrc(end).BoostedPower = 0;
    pfrc(end).Modulation = "1024QAM";
    pfrc(end).MCSIndex = 23;        % FRC TCR = 0.78,  MCS 23 in table 4 (1024QAM)

    % TS 38.101-2 A.3 DL RMC (FRC) (TDD definitions) 
    pfrc(end+1).Name = "DL-FRC-FR2-QPSK";
    pfrc(end).FR = "FR2";
    pfrc(end).DuplexMode = "TDD";
    pfrc(end).BoostedPercent = 100;
    pfrc(end).BoostedPower = 0;
    pfrc(end).Modulation = "QPSK";
    pfrc(end).MCSIndex = 4;         % FRC TCR = 1/3, MCS 4 in table 1 (64QAM)
        
    pfrc(end+1).Name = "DL-FRC-FR2-16QAM";  % Section empty in TS 38.101-2 v15.6.0 (RAN#84)
    pfrc(end).FR = "FR2";
    pfrc(end).DuplexMode = "TDD";
    pfrc(end).BoostedPercent = 100;
    pfrc(end).BoostedPower = 0;
    pfrc(end).Modulation = "16QAM";
    pfrc(end).MCSIndex = 13;        % FRC TCR = 1/2, MCS 13 (16QAM - 10-16 MCS) in table 1 (64QAM)
        
    pfrc(end+1).Name = "DL-FRC-FR2-64QAM";
    pfrc(end).FR = "FR2";
    pfrc(end).DuplexMode = "TDD";
    pfrc(end).BoostedPercent = 100;
    pfrc(end).BoostedPower = 0;
    pfrc(end).Modulation = "64QAM";
    pfrc(end).MCSIndex = 19;        % FRC TCR = 1/2, MCS 19 in table 1 (64QAM)
    
    % FRC common parameters
    [pfrc.ControlNRB] = deal(-1);       % Indicate that CORESET is not formally defined in the PDSCH FRC
    [pfrc.NumCCE] = deal(-1);           % Indicate that PDCCH is not formally defined in the PDSCH FRC
    
    % Combine sets of definitions into a single list
    p = [ptm pfrc];
    
end

% Get NRB and numerical BW/SCS from the FR and BW/SCS
function [nrb,bw,scs] = getValidNRB(fr,bw,scs)

    % Form property name for FR bandwidth table
    bwtable = [char(fr) 'BandwidthTable'];

    % Get standardized set of SCS values (text and numerical) for the FR
    scslist = reshape(hNRReferenceWaveformGenerator.(bwtable).Properties.RowNames,1,[]);
    scsnumbers = str2double(extractBefore(scslist,'k'));

    % Get standardized set of BW values (text and numerical) for the FR
    bwlist = hNRReferenceWaveformGenerator.(bwtable).Properties.VariableNames;
    bwnumbers = str2double(extractBefore(bwlist,'M'));

    % Validate and process the BW input
    bwsel = ones(1,length(bwnumbers),'logical');
    if ~isempty(bw)
        if isstring(bw) || ischar(bw)
            bwsel = strcmpi(erase(lower(bw),'bw_'),bwlist);  % Remove any 'BW_' prefix to support previous table column headings
        else
            bwlist = bwnumbers;
            bwsel = (bw == bwlist);
        end
        if ~any(bwsel)
            error('For %s, the channel bandwidth (%s) must be one of the set (%s).',fr,string(bw),join(string(bwlist),','));
        end
    end

    % Validate and process the SCS input
    scssel = ones(1,length(scsnumbers),'logical');
    if ~isempty(scs)
        if isstring(scs) || ischar(scs)
            scssel = strcmpi(scs,scslist);
        else
            scslist = scsnumbers;
            scssel = (scs == scslist);
        end
        if ~any(scssel)
            error('For %s, the subcarrier spacing (%s) must be one of the set (%s).',fr,string(scs),join(string(scslist),','));
        end
    end

    % Look-up the selected NRB from the BW table
    nrb = hNRReferenceWaveformGenerator.(bwtable){scssel,bwsel};

    if nrb == 0 | isnan(nrb)  %#ok<OR2> % Operand may be an array
        error('For %s, the combination of BW (%s) and SCS (%s) is not defined (see <a href="matlab:hNRReferenceWaveformGenerator.%sBandwidthTable">%sBandwidthTable</a>)',fr,string(bw),string(scs),fr,fr);
    end

    % Get numerical BW/SCS values
    bw = bwnumbers(bwsel);
    scs = scsnumbers(scssel);

end

% FR1 and FR2 NR-TM and FRC TDD-UL-DL-ConfigCommon configurations
function [config,txpattern] = getTDDConfiguration(fr,scs,link,dlfrc)

    % Default to downlink selections and NR-TM
    if nargin < 4
        dlfrc = 0;
        % Default to NR-TM configs, if (DL) FRC flag is not supplied
        if nargin < 3
            link = "downlink"; 
        end   
    end
    
    % TDD-UL-DL-ConfigCommon structures
    persistent patterns;
    if isempty(patterns)
        
        % TDD reference SCS must not be larger than any SCS of configured BWP, 
        % therefore TDD reference SCS must be less than or equal to the preset waveform SCS

        % NR-TM FR1 / DL FRC FR1
        % FR1 (15kHz,30kHz,60kHz)
        %
        % 3D1S1U, S=10D:2G:2U (5 slots @ 15kHz = 5ms)
        patterns(1).referenceSubcarrierSpacing = 15;
        patterns(1).dl_UL_TransmissionPeriodicity = 5;  % 5ms
        patterns(1).nrofDownlinkSlots = 3;
        patterns(1).nrofDownlinkSymbols = 10;
        patterns(1).nrofUplinkSlots = 1;
        patterns(1).nrofUplinkSymbols = 2;
        %   
        % 7D1S2U, S=6D:4G:4U (10 slots @ 30kHz = 5ms)
        patterns(2).referenceSubcarrierSpacing = 30;
        patterns(2).dl_UL_TransmissionPeriodicity = 5;  % 5ms
        patterns(2).nrofDownlinkSlots = 7;
        patterns(2).nrofDownlinkSymbols = 6;
        patterns(2).nrofUplinkSlots = 2;
        patterns(2).nrofUplinkSymbols = 4;  
        %
        % 14D2S4U, S=12D:4G:0U & S=0D:6G:8U (20 slots @ 60kHz = 5ms)
        patterns(3).referenceSubcarrierSpacing = 60;
        patterns(3).dl_UL_TransmissionPeriodicity = 5;  % 5ms
        patterns(3).nrofDownlinkSlots = 14;
        patterns(3).nrofDownlinkSymbols = 12;
        patterns(3).nrofUplinkSlots = 4;
        patterns(3).nrofUplinkSymbols = 8;  
        
        % NR-TM FR2 (FR2-1 and FR2-2)
        % FR2 (60kHz,120kHz,480kHz,960kHz)
        %
        % 60 kHz ref SCS
        % 3D1S1U, S=10D:2G:2U (5 slots @ 60kHz = 1.25ms)
        patterns(4).referenceSubcarrierSpacing = 60;
        patterns(4).dl_UL_TransmissionPeriodicity = 1.25;  % 1.25ms
        patterns(4).nrofDownlinkSlots = 3;
        patterns(4).nrofDownlinkSymbols = 10;
        patterns(4).nrofUplinkSlots = 1;
        patterns(4).nrofUplinkSymbols = 2;
        % 120 kHz ref SCS
        % 7D1S2U, S=6D:4G:4U (10 slots @ 120 kHz = 1.25ms)
        patterns(5).referenceSubcarrierSpacing = 120;
        patterns(5).dl_UL_TransmissionPeriodicity = 1.25;  % 1.25ms
        patterns(5).nrofDownlinkSlots = 7;
        patterns(5).nrofDownlinkSymbols = 6;
        patterns(5).nrofUplinkSlots = 2;
        patterns(5).nrofUplinkSymbols = 4;
        
        % DL FRC FR2 (60kHz,120kHz)
        % TS 38.101-2 Table A.3.3.1-1 (DL FRC general)
        % 
        % 60 kHz ref SCS
        % 3D1S1U, S=4D:6G:4U (5 slots @ 60kHz = 1.25ms)
        patterns(6).referenceSubcarrierSpacing = 60;
        patterns(6).dl_UL_TransmissionPeriodicity = 1.25;  % 1.25ms
        patterns(6).nrofDownlinkSlots = 3;
        patterns(6).nrofDownlinkSymbols = 4;
        patterns(6).nrofUplinkSlots = 1;
        patterns(6).nrofUplinkSymbols = 4;
        % 120 kHz ref SCS
        % 3D1S1U, S=10D:2G:2U (5 slots @ 120kHz = 0.625ms)
        patterns(7).referenceSubcarrierSpacing = 120;
        patterns(7).dl_UL_TransmissionPeriodicity = 0.625;  % 0.625ms
        patterns(7).nrofDownlinkSlots = 3;
        patterns(7).nrofDownlinkSymbols = 10;
        patterns(7).nrofUplinkSlots = 1;
        patterns(7).nrofUplinkSymbols = 2;
        
        % For UL FRC, TS 38.104 
        % 
        % Section 8 Conducted performance requirements      
        % Table 8.2.1.1-1 Test parameters for testing PUSCH (transform precoding disabled)
        % 15 kHz SCS: 3D1S1U, S=10D:2G:2U (5 slots = 5ms)
        % 30 kHz SCS: 7D1S2U, S=6D:4G:4U  (10 slots = 5ms)
        % Table 8.2.2.1-1: Test parameters for testing PUSCH (transform precoding enabled)
        % 15 kHz SCS: 3D1S1U, S=10D:2G:2U (5 slots = 5ms)
        % 30 kHz SCS: 7D1S2U, S=6D:4G:4U  (10 slots = 5ms)
        % 
        % 15 kHz -> pattern 1 above
        % 30 kHz -> pattern 2 above
        %
        % Section 11 Radiated performance requirements
        % FR1 cases as above (BS type 1-O)
        % FR2 cases (BS type 2-O)
        % Table 11.2.2.1.1-1 Test parameters for testing PUSCH (transform precoding disabled)
        % 60 kHz/120 kHz SCS: 3D1S1U, S=10D:2G:2U (5 slots = 1.25ms or 0.625ms, respectively for SCS)
        % Table 11.2.2.2.1-1: Test parameters for testing PUSCH (transform precoding enabled)
        % 60 kHz/120 kHz SCS: 3D1S1U, S=10D:2G:2U (5 slots = 1.25ms or 0.625ms, respectively for SCS)
        % 
        % 60 kHz -> pattern 7 above
        % 120 kHz -> pattern 7 above
        
    end

    % Select the config for the SCS, FR and NR-TM/DL FRC/Ul FRC type
    refscs = min(scs,120);   % The reference SCS here do not exceed 120 kHz so map any larger target SCS down to 120
    scsidx = 1+log2(refscs/15);
    if link == "downlink"
        index = scsidx + strcmpi(fr,'FR2')*(1 + 2*logical(dlfrc)); % Stagger the FR2 index start in the above table, for the 60kHz duplication in both FR1 and FR2
        config = patterns(index);
    else
        ulselection = [1 2 7 7];  % Patterns to be used for uplink FRC @ 15,30,60,120 kHz SCS respectively
        index = ulselection(scsidx);
        config = patterns(index);
        config.referenceSubcarrierSpacing = refscs;
        config.dl_UL_TransmissionPeriodicity = (1+double(scsidx==3))*config.dl_UL_TransmissionPeriodicity;  % Re-adjust the period (x 2) for 60kHz case
    end
    
    % Convert reference TDD config into one for target SCS

    % Scale the number of slots by the difference between the reference numerology and the target numerology    
    targetconfig = config;
    scsratio = scs/config.referenceSubcarrierSpacing; % Ratio of number of slots/symbols in target SCS numerology per slot/symbol in reference SCS numerology
    symperslot = 14;   % Assume we are mapping to normal CP numerology for calculation the mapped slots/symbols

    targetconfig.referenceSubcarrierSpacing = scs;

    % Map the DL part
    sdls = scsratio*config.nrofDownlinkSymbols;
    targetconfig.nrofDownlinkSlots = scsratio*config.nrofDownlinkSlots + fix(sdls/symperslot);
    targetconfig.nrofDownlinkSymbols = mod(sdls,symperslot);
    % Map the UL part
    suls = scsratio*config.nrofUplinkSymbols; 
    targetconfig.nrofUplinkSlots = scsratio*config.nrofUplinkSlots + fix(suls/symperslot);
    targetconfig.nrofUplinkSymbols = mod(suls,symperslot);

    config = targetconfig;

    % DL-UL period pattern, marked out in a vector (0=DL, 1=S containing DL, 2=UL)
    slotspertddperiod = config.dl_UL_TransmissionPeriodicity*fix(config.referenceSubcarrierSpacing/15);  % Number of slots in the TDD DL-UL period (1 slot in 1ms subframe, when @ 15kHz)
    txpattern = 2*ones(1,slotspertddperiod);                                                             % Mark all slots with 2 (to indicate 'not' the link direction of interest)
    
    if link == "downlink"
        txpattern(config.nrofDownlinkSlots+(1:config.nrofDownlinkSymbols~=0)) = 1;   % Mark downlink partial slots (with 1)
        txpattern(1:config.nrofDownlinkSlots) = 0;                                   % Mark downlink slots (with 0)
    else
        txpattern(end-config.nrofUplinkSlots+(0:config.nrofUplinkSymbols==0)) = 1;   % Mark uplink partial slots (with 1)
        txpattern(end-config.nrofUplinkSlots+1:end) = 0;                             % Mark uplink slots (with 0)
    end

end

% Get basic common downlink waveform parameter set (excludes any PDSCH definition)
function waveconfig = getCommonParameters(tdef,bw,scs,ncellid)

    % Use the name to identify the ref configuration as an FRC or NR-TM 
    frc = contains(tdef.Name,"FRC");
    % FRC
    %   SSB in slot 0 in a frame 
    %   CORESET PRB equals all NRB 
    %   PDCCH 'OCNG'
    % 
    % NR-TM 
    %   No SSB
    %   CORESET PRB equals 2
    %   PDCCH single instance

    % Get the max NRB from the BW and SCS
    [nrb,bw,scs] = getValidNRB(tdef.FR,bw,scs);
    
    tdd = strcmpi(tdef.DuplexMode,'TDD');   % Local bool for TDD or FDD
   
    % Carrier-level parameters
    waveconfig.Name = tdef.Name;         % Label the configuration
    waveconfig.NCellID = ncellid;        % Cell identity
    waveconfig.ChannelBandwidth = bw;    % Channel bandwidth (MHz)
    waveconfig.FrequencyRange = tdef.FR; % 'FR1' or 'FR2'
    waveconfig.NumSubframes = 10+10*tdd; % Number of 1ms subframes in generated waveform (1,2,4,8 slots per 1ms subframe, depending on SCS, and 10 subframes in a 10ms frame)  
    waveconfig.Windowing = 0;            % Windowing amount (set it to 0 in these presets i.e. no windowing)
    waveconfig.DisplayGrids = 0;         % Display the resource grids after signal generation

    % SCS carrier configuration
    % Define a single SCS-specific carrier for the waveform
    carriers(1).SubcarrierSpacing = scs;
    carriers(1).NRB = nrb;
    carriers(1).RBStart = 0;

    % SS burst configuration
    % This is not enabled for the NR test models since the DM-RS are expected
    % to be used for synchronization
    % 
    % Case A (FR1)
    %  f <= 3 GHz and 3 < f <= 6 GHz
    %  SSB SCS = 15 kHz
    % Case B (FR1)
    %  f <= 3 GHz and 3 < f <= 6 GHz
    %  SSB SCS = 30 kHz
    % Case C (FR1)
    %  f <= 3 GHz and 3 < f <= 6 GHz
    %  SSB SCS = 30 kHz
    % Case D (FR2)
    %  f > 6 GHz
    %  SSB SCS = 120 kHz
    % Case E (FR2)
    %  f > 6 GHz
    %  SSB SCS = 240 kHz
    % Case F (FR2)
    %  f > 6 GHz
    %  SSB SCS = 480 kHz
    % Case G (FR2)
    %  f > 6 GHz
    %  SSB SCS = 960 kHz
    %
    % FR1: 'Case B' = 30 kHz SSB SCS, FR1 carriers are 15 / 30 / 60 kHz SCS
    % FR2: 'Case D' = 120 kHz SSB SCS, FR2 carriers are 60 / 120 / 240 kHz SCS
    
    % Selected mapping table
    % scs      =  15, 30, 60, 120, 240, 480, 960           % If 60kHz and FR2 then "D" (adjusted below)
    % ssbscs   =  15, 30, 30, 120, 240, 480, 960
    ssbcases   = ["A","B","B","D","E","F","G"];
    bmplengths = [ 4,  4,  4,  64, 64, 64, 64];
           
    scsind = log2(scs/15)+1;
    % If 5MHz carrier @ 30kHz SCS on the data then we need to go back to SSB case A @ 15 kHz
    % since case B (30kHz) would not fit into PRB in the 5MHz carrier
    if nrb==11 && scs==30
        scsind = scsind-1;
    end
    % If FR2 and 60kHz data SCS then we need to go up to the 120kHz SSB case
    if tdef.FR == "FR2" && scs==60
        scsind = scsind+1;
    end
    
    ssbbp = strcat("Case ",ssbcases(scsind));
    ssbbmap = zeros(1,bmplengths(scsind));   % Create bitmap of the appropriate length
    ssbbmap(1) = 1;                          % Enable the first SSB occurrence 
    
    ssburst = [];
    ssburst.Enable = frc;                    % Enable SS Burst if an FRC
    ssburst.BlockPattern = ssbbp;            % Case B (30kHz) subcarrier spacing
    ssburst.SSBTransmitted = ssbbmap;        % Bitmap indicating blocks transmitted in a 5ms half-frame burst
    ssburst.SSBPeriodicity = 10;             % SS burst set periodicity in ms (5, 10, 20, 40, 80, 160)
    ssburst.FrequencySSB = 0*5000;           % Frequency offset of SS burst (Hz), relative to waveform center (multiples of 5kHz)
    ssburst.Power = 0;                       % Power scaling in dB

    % Bandwidth part configuration
    % Define a single BWP with spans the whole SCS carrier
    bwp(1).SubcarrierSpacing = scs;          % BWP Subcarrier Spacing
    bwp(1).CyclicPrefix = 'Normal';          % BWP Cyclic prefix
    bwp(1).NRB = nrb;                        % Size of BWP
    bwp(1).RBOffset = 0;                     % Position of BWP in SCS carrier

    % Downlink carrying slots for the PDCCH
    if strcmpi(tdef.DuplexMode,'FDD')
        % For FDD, every slot is DL
        dlallocatedslots = 0;
        dlallocatedperiod = 1; 
    else
        % For TDD, use downlink and downlink partial slots across the TDD DL-UL repetition period
        [~,dltxpattern] = getTDDConfiguration(tdef.FR,scs,"downlink",frc);
        dlallocatedslots = find(dltxpattern<(2-frc))-1;      % Schedule in DL slots only for FRC, and both DL and special slots for the NR-TM
        dlallocatedperiod = length(dltxpattern);     
    end
    
    % Expand the CORESET/PDCCH related parameterization
    ncontrolsymbs = 2;                    % Number of OFDM symbols used for control channel/CORESET
    ncontrolrb = tdef.ControlNRB;         % Number of RB used for CORESET
    if ncontrolrb < 0                     % If specified as negative then use the entire bandwidth
        ncontrolrb = nrb;
    end
    ncontrolcce = tdef.NumCCE;            % Number of CCE used for PDCCH instance in the CORESET
    if ncontrolcce < 0                    % If specified as negative then it's not a test model definition
       ncontrolcce = 1;                   % One CCE equals 6 REG (which is one RB during one OFDM symbol)
       % If trying to span the CORESET then ncontrolcce = fix(ncontrolrb/6)*ncontrolsymbs
    end
    
    % CORESET/search space configuration
    coreset(1).Duration = ncontrolsymbs;            % CORESET symbol duration (1,2,3)
    coreset(1).AllocatedSymbols = 0;                % First symbol of each CORESET monitoring opportunity in a slot
    coreset(1).AllocatedSlots =  dlallocatedslots;  % Allocated slots within a period
    coreset(1).AllocatedPeriod = dlallocatedperiod; % Allocated slot period (empty implies no repetition)
    coreset(1).AllocatedPRB = 0:ncontrolrb-1;       % 6 REG sized indices, relative to BWP
    coreset(1).CCEREGMapping = 'noninterleaved';    % CCE-to-REG mapping, 'interleaved' or 'noninterleaved'
    coreset(1).REGBundleSize = 2;                   % REG bundle size, L (2,6) or (3,6)
    coreset(1).InterleaverSize = 2;                 % Interleaver size, R (2,3,6)
    coreset(1).ShiftIndex = 0;                      % Shift index (0...274)
      
    % PDCCH Configuration
    pdcch(1).Enable = tdef.NumCCE>0;        % Enable PDCCH config when defined in the test waveform
    pdcch(1).BWP = 1;                       % Bandwidth part ID for the PDCCH sequence
    pdcch(1).CORESET = 1;                   % CORESET ID which carries this PDCCH sequence
    pdcch(1).Power = 0;                     % Power scaling in dB
    pdcch(1).AllocatedSearchSpaces = 0;     % Index of monitoring occasion in the overall CORESET/search space pattern
    pdcch(1).AllocatedPeriod = 1;           % Allocation slot period (empty implies no repetition of pattern)
    pdcch(1).NumCCE = ncontrolcce;          % Number of CCE used by PDCCH
    pdcch(1).StartCCE = 0;                  % Starting CCE of PDCCH
    pdcch(1).RNTI = 0;                      % RNTI
    pdcch(1).NID = ncellid;                 % PDCCH and DM-RS scrambling NID
    pdcch(1).PowerDMRS = 0;                 % Additional power boosting in dB
    pdcch(1).EnableCoding = 0;              % Enable DCI coding for the PDCCH
    pdcch(1).DataBlkSize = 20;              % DCI payload size
    pdcch(1).DataSource = 0;                % DCI data source

    % Break the total number of control elements to use into separate PDCCH
    targetal = 16;
    ntargetdci = fix(ncontrolcce/targetal);  % Number of PDCCH available at the target AL
    naremdci = mod(ncontrolcce,targetal);    % Number of additional PDCCH at AL=1 to fill the CORESET
    
    % Calculate the starting NCCE for instances
    alset = [ones(1,ntargetdci)*targetal ones(1,naremdci)];  
    startset = cumsum([0 alset(1:end-1)]); 
    alset = num2cell(alset); 
    startset = num2cell(startset);
    
    % Create separate PDCCH at the starting positions
    pdcch = repmat(pdcch,1,length(alset));
    [pdcch.NumCCE] = alset{:};
    [pdcch.StartCCE] = startset{:};
    rntiset = num2cell(0:length(alset)-1);
    [pdcch.RNTI] = rntiset{:};
    
    % CSI-RS Configuration
    csirs.Enable = 0;
    csirs.BWP = 1;
    csirs.Power = 0;
    csirs.CSIRSType = 'nzp';
    csirs.RowNumber = 1;
    csirs.Density = 'three';
    csirs.SubcarrierLocations = 0;
    csirs.NumRB = nrb;
    csirs.RBOffset = 0;
    csirs.SymbolLocations = 0;
    csirs.AllocatedSlots = 0;
    csirs.AllocatedPeriod = 1;
    csirs.NID = ncellid;
    
    % Assign channel specific sub-structures
    waveconfig.SSBurst = ssburst;
    waveconfig.Carriers = carriers;
    waveconfig.BWP = bwp;
    waveconfig.CORESET = coreset;
    waveconfig.PDCCH = pdcch;
    waveconfig.CSIRS = csirs;

end

% Add FRC PDSCH parameter definition to the parameter set
function waveconfig = getFRCPDSCHParameters(tdef,waveconfig)
    
    % The parameter set below is geared towards TS 38.101 DL FRC RMCs
    % 
    % TS 38.101-1 User Equipment (UE) radio transmission and reception; Part 1: Range 1 Standalone
    % 
    % A.3 DL reference measurement channels
    % A.3.1 General
    % A.3.2 DL reference measurement channels for FDD
    % A.3.2.1 General
    % A.3.2.2 FRC for receiver requirements for QPSK
    % A.3.2.3 FRC for maximum input level for 64QAM
    % A.3.2.4 FRC for maximum input level for 256 QAM
    % A.3.3 DL reference measurement channels for TDD
    % A.3.3.1 General
    % A.3.3.2 FRC for receiver requirements for QPSK
    % A.3.3.3 FRC for maximum input level for 64QAM
    % A.3.3.4 FRC for maximum input level for 256 QAM
    % 
    % TS 38.101-2 User Equipment (UE) radio transmission and reception; Part 2: Range 2 Standalone
    % 
    % DL reference measurement channels
    % A.3.1 General
    % A.3.2 Void
    % A.3.3 DL reference measurement channels for TDD
    % A.3.3.1 General
    % A.3.3.2 FRC for receiver requirements for QPSK
    % A.3.3.3 FRC for receiver requirements for 16QAM
    % A.3.3.4 FRC for receiver requirements for 64QAM
    
    % Get the NRB/SCS values from the existing parameter set 
    nrb = waveconfig.Carriers.NRB;
    scs = waveconfig.BWP.SubcarrierSpacing;
    
    % Get basic full band PDSCH parameter definition 
    % Pick NCellID up from the overall wavegen parameters
    pdsch = getPDSCHCommonParameters(tdef.FR,nrb,waveconfig.NCellID,tdef.Modulation,tdef.MCSIndex);
        
    % Complete entries from TS 38.101-1 Table A.3.1-1 and TS 38.101-2 Table A.3.1-1 (Common reference channel parameters)
    pdsch.AllocatedSymbols = 2:max(pdsch.AllocatedSymbols);    % Leave the first two symbols for the CORESET
    pdsch.DMRSAdditionalPosition = 2;    % 2 additional DM-RS symbols (single symbol length)
    pdsch.NumCDMGroupsWithoutData = 2;   % Disable any FDM between DM-RS and PDSCH
    pdsch.PowerDMRS = 3;                 % Additional power boosting in dB
    
    % Set the PT-RS and resulting TBS overhead as required
    pdsch.EnablePTRS = (tdef.FR=="FR2") && (tdef.Modulation=="64QAM");
    pdsch.PTRSTimeDensity = 1;
    pdsch.PTRSFrequencyDensity = 2;
    pdsch.Xoh_PDSCH = 6*pdsch.EnablePTRS;   % Overhead to be used if PT-RS enabled
    
    % Identify the starting slot of the reference PDSCH sequence in a frame
    if  waveconfig.SSBurst.BlockPattern <= "Case C" && ...    % SSB SCS <= 30kHz AND
            waveconfig.Carriers.SubcarrierSpacing == 60       % Data SCS == 60kHz
        startslot = 2;
    else
        startslot = 1;
    end
    
    % Finalize the PDSCH slot allocation and period part
    pdsch = setFRCSlotAllocation(pdsch,tdef,scs,startslot);

    % Add PDSCH config to the parameter set
    waveconfig.PDSCH = pdsch;  
    
end

% Set the time allocation part (allocated slots and period, accounting for TDD)
function pxsch = setFRCSlotAllocation(pxsch,tdef,scs,startslot)

    if nargin < 4
        startslot = 0;
    end
    
    % Cue the link direction from the FRC definition name
    links = ["uplink","downlink"];
    linkdirection = links(1+double(startsWith(tdef.Name,'DL')));
     
    spf = 10*scs/15;
    pxsch.AllocatedSlots = startslot:spf-1; % All slots in a frame except 0 (and 1 if FR1 60kHz (since 30kHz SSB))
    pxsch.AllocatedPeriod = spf;            % Slots per frame
           
    % If TDD then update the slot allocations and expand the sequence definitions
    % to deal with the partial DL slots
    if strcmpi(tdef.DuplexMode,'TDD')
     % Get the FRC defined DL-UL TDD configuration
        [~,tp] = getTDDConfiguration(tdef.FR,scs,linkdirection,1);   % Get TDD DL-UL slot period mask
        
        tddmask = repmat(tp,1,ceil(pxsch.AllocatedPeriod/length(tp)));          % Repeat TDD mask period so that it covers the allocation period 
        tddmask = tddmask(1:pxsch.AllocatedPeriod);                             % Extract mask across allocated period only, where the DL/UL only slots = 0 (i.e. not special, or the 'other' link direction)
        tddmask(pxsch.AllocatedSlots+1) = tddmask(pxsch.AllocatedSlots+1)-1;    % Set allocated DL/UL slots only to -1
        
        pxsch.AllocatedSlots = find(tddmask<0)-1;
        
    end
end

% Add PDSCH parameter definition to the parameter set
function waveconfig = getTestModelPDSCHParameters(tdef,waveconfig,sv)

    % Using Rel 15/16/17 badged versions of 38.141: 15.0.0, 15.1.0, 15.2.0... 15.7.0... 16.7.0... 17.8.0
    % Identify the std version of NR-TM required
    supportedversions = ["15.1.0", "15.2.0", "15.7.0", "16.7.0", "17.8.0"];
    vcomp = strcmpi(sv,supportedversions);
    if ~any(vcomp)
        error("The NR-TM standard version (%s) of TS 38.141-1 and TS 38.141-2 must be one of (%s)",string(sv),join(supportedversions,', '));
    end
    stdversion = find(vcomp);

    % The parameter set below is geared to the NR-TMs
    % TS 38.141-1 Section 4.9.2 (FR1) and TS 38.141-2 Section 4.9.2 (FR2)

    % Get the NRB/SCS values from the existing parameter set 
    nrb = waveconfig.Carriers.NRB;
    scs = waveconfig.BWP.SubcarrierSpacing;
    
    % Get basic full band PDSCH parameter definition 
    % Pick NCellID up from the overall wavegen parameters
    pdsch = getPDSCHCommonParameters(tdef.FR,nrb,waveconfig.NCellID,tdef.Modulation(1));
    
    % Later versions use PN23 data for PDCCH/PDSCH rather than 0's
    if stdversion>2
        waveconfig.PDCCH.DataSource = 'PN23';
        pdsch.DataSource = 'PN23';
    end
    
    if stdversion>1
        % From TS 38.141 V15.2.0, an additional PDSCH (RNTI=2) was introduced to define the first three PRB which
        % time interleaves with the PDCCH
        
        % Starting with the full band preset, adjust the frequency resources (NRB-3) to create the RNTI = 0 version
        pdsch.AllocatedPRB = 3:nrb-1;
        pdsch.Name = sprintf("Partial band PDSCH sequence with %s modulation scheme (target, RNTI = %d)",pdsch.Modulation, pdsch.RNTI);
        
        % Create the additional PDSCH for PRB 0:2 (RNTI=2). This is the PDSCH which is time multiplexed with the PDCCH
        pdsch(2) = pdsch;
        pdsch(2).AllocatedPRB = 0:2;
        pdsch(2).AllocatedSymbols = 2:max(pdsch(2).AllocatedSymbols);
        pdsch(2).RNTI = 2;
        pdsch(2).Name = sprintf("Partial band PDSCH sequence with %s modulation scheme (target, RNTI = %d)",pdsch(2).Modulation, pdsch(2).RNTI);   
    end 
    
    % Modify and expand the PDSCH parameter set to enable dynamic aspects 
    % of the NR-TMs (boosting/deboosting and single PRB ramp)
    if tdef.BoostedPercent < 100
        if tdef.BoostedPercent == 1
            % Single PRB ramp case
            slotsperframe = 10*fix(scs/15);       
            % Adjust the allocation 
            pdsch(1).RNTI = 2*(stdversion>1);
            pdsch(1).AllocatedSymbols = 2*(stdversion>1):max(pdsch(1).AllocatedSymbols);
            pdsch(1).AllocatedPRB = 0;
            pdsch(1).AllocatedSlots = 3*(0:ceil(slotsperframe/3)-1);
            pdsch(1).AllocatedPeriod = slotsperframe;
            pdsch(1).Name = "PDSCH sequence for lower PRB (3n slots)"; 
            % Add a second PDSCH using a different PRB  
            pdsch(2) = pdsch(1);
            pdsch(2).AllocatedPRB = fix(nrb/2);
            pdsch(2).AllocatedSlots = 1+3*(0:ceil((slotsperframe-1)/3)-1);
            pdsch(2).Name = "PDSCH sequence for middle PRB (3n+1 slots)"; 
            % Add a third PDSCH using a different PRB  
            pdsch(3) = pdsch(1);
            pdsch(3).AllocatedPRB = nrb-1;
            pdsch(3).AllocatedSlots = 2+3*(0:ceil((slotsperframe-2)/3)-1);
            pdsch(3).Name = "PDSCH sequence for upper PRB (3n+2 slots)"; 
        else
            % Boosted/deboosted 3 PDSCH models
            
            % Later versions always use QPSK for the power balancing PDSCH
            if stdversion>2
                tdef.Modulation(2) = "QPSK";
            end  
            
            % From V15.2.0, exclude the first PRG from the target always
            exprg1 = stdversion>1;
    
            % Boosted/deboosted split
            % Take the default definition, modify it and created an additional one 
            % which interleaves with it in frequency
            bwpalloc = pdsch(1).AllocatedPRB;
            bwpnrb =  length(bwpalloc);    
            P = getPRGSize(bwpnrb,1);      % PRG size, P 
                   
            % Identify the target PRB (PRG) (RNTI=0 if boosted or 1 if deboosted) 
            bpercent = tdef.BoostedPercent/100;                    
            nprgmax =  fix((bwpnrb + mod(bwpalloc(1),P) - P*exprg1)/P);
            nprg = min(fix(bpercent*bwpnrb/P),(nprgmax-mod(nprgmax,2))/2+1);   
                    
            % PRG indices of the target
            lastprg = fix((bwpnrb + mod(bwpalloc(1),P))/P)-1;
            prg = [exprg1+(0:2:2*(nprg-2)),lastprg];        % PRB of first PRG in 'BWP'
            prgprb = P*fix(bwpalloc(1)/P) + P*prg;
            boostedprb = reshape(prgprb + (0:P-1)',1,[]);   % Expand PRG into associated PRB
                
            pdsch(1).AllocatedPRB = boostedprb;
            pdsch(1).Power = tdef.BoostedPower;
            pdsch(1).RNTI = double(tdef.BoostedPower<0);
            pdsch(1).Name = sprintf("Partial band PDSCH sequence with %s modulation scheme (target, RNTI = %d)",pdsch(1).Modulation, pdsch(1).RNTI);
            
            % Set the modulation to QPSK if RNTI=2 PDSCH defined
            if length(pdsch)==2 && pdsch(2).RNTI==2
                pdsch(2).Modulation='QPSK';
            end
                    
            % Add the power-balanced PRB (RNTI=1 if boosted or 0 if deboosted)
            deboostedprb = setdiff(bwpalloc,boostedprb);
            pdsch(end+1) = pdsch(1);
            pdsch(end).AllocatedPRB = deboostedprb;
            pdsch(end).Power = 10*log10((bwpnrb-10^(tdef.BoostedPower/10)*P*nprg)/(bwpnrb-P*nprg));
            pdsch(end).Modulation = tdef.Modulation(2);
            pdsch(end).RNTI = double(tdef.BoostedPower>=0);                  
            pdsch(end).Name = sprintf("Partial band PDSCH sequence with %s modulation scheme (non-target, RNTI = %d)",pdsch(end).Modulation, pdsch(end).RNTI);
            
        end
    end
     
    % If TDD then update the slot allocations and expand the sequence definitions
    % to deal with the partial DL slots
    if strcmpi(tdef.DuplexMode,'TDD')
        
        % Get the NR-TM defined DL-UL TDD configuration
        [pattern,tp] = getTDDConfiguration(tdef.FR,scs);
        % The single varying PRB case needs special handling to get the sequences right
        if tdef.BoostedPercent == 1
            
            % Create the PRB allocation (frequency) schedule across a *frame*
            fp = repmat(0:length(pdsch)-1,1, ceil(10*fix(scs/15)/3));    % Sequence index/PRB type used in consecutive slot numbers
            fp = fp(1:10*fix(scs/15));
            
            % Create an expanded sequence of indices which identify the sequence definition
            % in use in each slot
            p1 = length(tp);
            p2 = length(fp);
            tpl = lcm(p1,p2);  % Length of combined TDD and varying PRB patterns

            % Start with three different PRB allocations in sequence
            % Across the entire pattern, each different allocation will be expanded
            % 
            % The fp pattern indicates original sequence selection associated
            % each consecutive block of three slots
            % Combine the sequences, expanding the slot type number ()
            ipattern = repmat(length(pdsch)*tp,1,tpl/p1)+repmat(fp,1,tpl/p2);
            sequenceslotindices = arrayfun(@(x)find(x==ipattern)-1,0:5,'UniformOutput',0);
     
        else            
            % Downlink slots already created but adjust the allocated slots and period for 
            % the full DL-UL period
            tpl = length(tp);
            sequenceslotindices = [ repmat( {0:pattern.nrofDownlinkSlots-1 },1,length(pdsch)), ...
                   repmat( {pattern.nrofDownlinkSlots+(0:ceil(pattern.nrofDownlinkSymbols/14)-1)},1,length(pdsch))];
        end
 
        % Create a set of sequence specifications for the partial DL slots
        pdschpart = pdsch;
        names = [pdsch.Name] + " (Full downlink slots)";        % Extend names to indicate full slot sequences
        [pdsch.Name] = deal(names{:});
        names = [pdschpart.Name] + " (Partial downlink slots)"; % Extend names to indicate partial slot sequences
        [pdschpart.Name] = deal(names{:});
        % Adjust the end of each allocation for the special slot
        allocsymbs = cellfun( @(x)min(x):pattern.nrofDownlinkSymbols-1, {pdschpart.AllocatedSymbols} , 'UniformOutput', false);
        [pdschpart.AllocatedSymbols] = deal(allocsymbs{:});
        pdsch = [pdsch pdschpart];        
        % Update the sequence periods and slot allocations for the sequences
        [pdsch.AllocatedPeriod] = deal(tpl);
        [pdsch.AllocatedSlots] = deal(sequenceslotindices{:});  % Assignment updates both full and partial DL slot sequences   
    end

    % Add PDSCH config to the parameters
    waveconfig.PDSCH = pdsch;
    
end

% Get basic full band, full slot PDSCH parameter definition
function pdsch = getPDSCHCommonParameters(fr,nrb,ncellid,mod,mcs) %#ok<INUSD>

    % Default to coding disabled (NR-TM)
    if nargin < 5
       mcs = -1;        % -1 signals no coding
    end

    % Basic PDSCH parameters
    pdsch.Name = "Full-band PDSCH sequence";% Description of this PDSCH sequence
    pdsch.Enable = 1;                       % Enable PDSCH config
    pdsch.BWP = 1;                          % Bandwidth part
    pdsch.Power = 0;                        % Power scaling in dB
    pdsch.EnableCoding = mcs>=0;            % Enable DL-SCH coding for the PDSCH
    if pdsch.EnableCoding
        tcr = lookupPDSCHTCR(mcs,mod);
        datasource = "PN9-ITU";
    else
        tcr = 0.4785;
        datasource = 0;
    end
    pdsch.DataSource = datasource;       % Transport block data source
    pdsch.TargetCodeRate = tcr;          % Code rate used to calculate transport block sizes
    pdsch.Xoh_PDSCH = 0;                 % Rate matching overhead
    pdsch.Modulation = mod;              % 'QPSK', '16QAM', '64QAM', '256QAM'
    pdsch.NLayers = 1;                   % Number of PDSCH layers
    pdsch.RVSequence = 0;                % RV sequence to be applied cyclically across the PDSCH allocation sequence (General sequence is [0,2,3,1])
    pdsch.VRBToPRBInterleaving = 0;      % Interleaving enable/disable for PDSCH
    pdsch.VRBBundleSize = 2;             % Bundle size for the interleaving (2,4)

    % Allocation part
    % Full band, full slot
    pdsch.AllocatedSymbols = 0:13;      % Range of symbols in a slot
    pdsch.AllocatedSlots = 0;           % Allocated slot indices in a period
    pdsch.AllocatedPeriod = 1;          % Allocation period in slots (empty implies no repetition of allocated slots pattern)
    pdsch.AllocatedPRB = 0:nrb-1;       % PRB allocation
    pdsch.RNTI = 0;                     % RNTI
    pdsch.NID = [];                     % Scrambling for PDSCH data part

    % Reserved resources 
    % Reserve 3 PRB by 2 symbols in each slot for PDCCH/CORESET across any of the PDSCH
    pdsch.RateMatch(1).CORESET = [];                     % Rate matching pattern, defined by one CORESET
    pdsch.RateMatch(1).Pattern.AllocatedPRB = 0:2;       % Rate matching pattern, defined by set of 'bitmaps'
    pdsch.RateMatch(1).Pattern.AllocatedSymbols = 0:1;
    pdsch.RateMatch(1).Pattern.AllocatedSlots = 0;
    pdsch.RateMatch(1).Pattern.AllocatedPeriod = 1;

    % DM-RS and antenna configuration (TS 38.211 section 7.4.1.1)
    pdsch.PortSet = [];                   % DM-RS ports to use for the layers
    pdsch.PDSCHMappingType = 'A';         % PDSCH mapping type ('A'(slot-wise),'B'(non slot-wise))
    pdsch.DMRSTypeAPosition = 2;          % Mapping type A only. First DM-RS symbol position (l0=2,3)
    pdsch.DMRSLength = 1;                 % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
    pdsch.DMRSAdditionalPosition = double(fr=="FR1");% Additional DM-RS symbol positions (max range 0...3)
    pdsch.DMRSConfigurationType = 1;      % DM-RS configuration type (1,2)
    pdsch.NumCDMGroupsWithoutData = 1;    % CDM groups without data
    pdsch.NIDNSCID = [];                  % PDSCH DM-RS scrambling identity (0...65535)
    pdsch.NSCID = 0;                      % PDSCH DM-RS scrambling initialisation (0,1)
    pdsch.PowerDMRS = 0;                  % Additional power boosting in dB
 
    % For reference, earlier standard names (R2018b) as below
    % pdsch.DL_DMRS_typeA_pos = 2;          % Mapping type A only. First DM-RS symbol position (l0=2,3)
    % pdsch.DL_DMRS_max_len = 1;            % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
    % pdsch.DL_DMRS_add_pos = (fr=="FR1");  % Additional DM-RS symbol positions (max range 0...3)
    % pdsch.DL_DMRS_config_type = 1;        % DM-RS configuration type (1,2)
    
    % PT-RS configuration
    pdsch.EnablePTRS = (fr=="FR2");         % Enable or disable the PT-RS (1 or 0)
    pdsch.PTRSTimeDensity = 4;              % Time density (L_PT-RS) of PT-RS (1,2,4)
    pdsch.PTRSFrequencyDensity = 2;         % Frequency density (K_PT-RS) of PT-RS (2,4)
    pdsch.PTRSREOffset = "00";              % PT-RS resource element offset ('00','01','10','11')
    pdsch.PTRSPortSet = min(pdsch.PortSet); % PT-RS antenna ports must be a subset of DM-RS ports
    pdsch.PowerPTRS = 0;                    % Additional power boosting in dB
        
end

% TS 38.214 Section 5.1.2.2.1 Downlink resource allocation type 0
function p = getPRGSize(nrb,config)
    
    % Bandwidth Part Size Configuration 1  Configuration 2
    % 1-36                2                4
    % 37-72               4                8
    % 73-144              8                16
    % 145-275             16               16
    
    p = min(16,(1+(config==2))*2^find(nrb <= [36 72 144 275],1));
    
end

% Get the nearest signalable target coding rate to the input rate, using
% TS 38.214 Table 5.1.3.1-1 (qam64), Table 5.1.3.1-2 (qam256) and Table 5.1.3.1-4 (qam1024)
function tcr = lookupPDSCHTCR(mcs,mod)
   
    % TS 38.214 Table 5.1.3.1-1: MCS index table 1 for PDSCH (qam64)
    table1 = [...
        0	2	120	0.2344
        1	2	157	0.3066
        2	2	193	0.3770
        3	2	251	0.4902
        4	2	308	0.6016
        5	2	379	0.7402
        6	2	449	0.8770
        7	2	526	1.0273
        8	2	602	1.1758
        9	2	679	1.3262
        10	4	340	1.3281
        11	4	378	1.4766
        12	4	434	1.6953
        13	4	490	1.9141
        14	4	553	2.1602
        15	4	616	2.4063
        16	4	658	2.5703
        17	6	438	2.5664
        18	6	466	2.7305
        19	6	517	3.0293
        20	6	567	3.3223
        21	6	616	3.6094
        22	6	666	3.9023
        23	6	719	4.2129
        24	6	772	4.5234
        25	6	822	4.8164
        26	6	873	5.1152
        27	6	910	5.3320
        28	6	948	5.5547];

    % TS 38.214 Table 5.1.3.1-2: MCS index table 2 for PDSCH (qam256)
    table2 = [...
        0	2	120	0.2344;
        1	2	193	0.3770
        2	2	308	0.6016
        3	2	449	0.8770
        4	2	602	1.1758
        5	4	378	1.4766
        6	4	434	1.6953
        7	4	490	1.9141
        8	4	553	2.1602
        9	4	616	2.4063
        10	4	658	2.5703
        11	6	466	2.7305
        12	6	517	3.0293
        13	6	567	3.3223
        14	6	616	3.6094
        15	6	666	3.9023
        16	6	719	4.2129
        17	6	772	4.5234
        18	6	822	4.8164
        19	6	873	5.1152
        20	8	682.5	5.3320
        21	8	711	5.5547
        22	8	754	5.8906
        23	8	797	6.2266
        24	8	841	6.5703
        25	8	885	6.9141
        26	8	916.5	7.1602
        27	8	948	7.4063];

    % TS 38.214 Table 5.1.3.1-3: MCS index table 3 for PDSCH (qam64LowSE)
    table3 = [...
        0	2	30	0.0586
        1	2	40	0.0781
        2	2	50	0.0977
        3	2	64	0.1250
        4	2	78	0.1523
        5	2	99	0.1934
        6	2	120	 0.2344
        7	2	157	 0.3066
        8	2	193	 0.3770
        9	2	251	 0.4902
        10	2	308	 0.6016
        11	2	379	 0.7402
        12	2	449	 0.8770
        13	2	526	 1.0273
        14	2	602	 1.1758
        15	4	340	 1.3281
        16	4	378	 1.4766
        17	4	434	 1.6953
        18	4	490	 1.9141
        19	4	553	 2.1602
        20	4	616	 2.4063
        21	6	438	 2.5664
        22	6	466	 2.7305
        23	6	517	 3.0293
        24	6	567	 3.3223
        25	6	616	 3.6094
        26	6	666	 3.9023
        27	6	719	 4.2129
        28	6	772	 4.5234]; %#ok<NASGU>
     
    % TS 38.214 Table 5.1.3.1-4: MCS index table 4 for PDSCH (qam1024)
    table4 = [...
        0	2	120	0.2344
        1	2	193	0.3770
        2	2	449	0.8770
        3	4	378	1.4766
        4	4	490	1.9141
        5	4	616	2.4063
        6	6	466	2.7305
        7	6	517	3.0293
        8	6	567	3.3223
        9	6	616	3.6094
        10	6	666	3.9023
        11	6	719	4.2129
        12	6	772	4.5234
        13	6	822	4.8164
        14	6	873	5.1152
        15	8	682.5	5.3320
        16	8	711	5.5547
        17	8	754	5.8906
        18	8	797	6.2266
        19	8	841	6.5703
        20	8	885	6.9141
        21	8	916.5	7.1602
        22	8	948	7.4063
        23	10	805.5	7.8662
        24	10	853	8.3301
        25	10	900.5	8.7939
        26	10	948	9.2578];
        % 27	2	reserved
        % 28	4	reserved
        % 29	6	reserved
        % 30	8	reserved
        % 31	10	reserved

    % Turn the modulation type into BPS for the original table lookup
    bps = sum(strcmpi(mod,["QPSK","16QAM","64QAM","256QAM","1024QAM"]).*[2 4 6 8 10]);  
    if bps == 10  
        actable = table4;   % qam1024 table
    elseif bps == 8
        actable = table2;   % qam256 table
    else
        actable = table1;   % qam64 table
    end
    tcr = actable(mcs+1,3)/1024;
end

function tdef = getULFRCDefinition(rc,bw,scs)
% Uplink 5G NR Fixed Reference Channels (FRC) as per TS 38.104 Annex A
%
%    CFG = getULFRCDefinition(rc,bw,scs) generates a configuration for the
%    Fixed Reference Channel (FRC) RC. RC must follow a G-FRX-AY-Z format,
%    where X is the frequency-range number (1 or 2), Y is the order of the
%    MCS used in Annex A of TS 38.104 (i.e., Y is 1-5) and Z is the index 
%    of the FRC for the given MCS. The range of valid Z values depends on
%    the Frequency Range and the MCS (X and Y).

    % Test requirements...
    % 
    % Section 7 - Conducted receiver characteristics
    %
    % Section 8 - Conducted performance requirements
    % - BS type 1-C or BS type 1-H
    % Section 8.2 - Performance requirements for PUSCH
    % 
    % Section 11 - Radiated performance requirements
    % - BS type 1-O or BS type 2-O
    % Section 11.2 - Performance requirements for PUSCH

    % Parse input, e.g., 'G-FR1-A3-18'
    rc = char(rc);
    FR = rc(3:5);
    frNum = str2double(FR(end));
    secNum = str2double(rc(8));
    frcNum = str2double(rc(10:end));

    % Build the FRC definition
    tdef.Name = rc;
    tdef.FR = string(FR);

    % Subcarrier spacing (SCS) - Row #1
    tdef.SubcarrierSpacing = getSCSUplinkFRC(frNum, secNum, frcNum);
    % Number of allocated PUSCH resource blocks - Row #2
    tdef.AllocatedNRB = getNRBUplinkFRC(frNum, secNum, frcNum);

    % Modulation - Row #4:
    switch secNum
      case {1, 3}
        tdef.Modulation = "QPSK";
      case {2, 4}
        tdef.Modulation = "16QAM";
      otherwise % A5
        tdef.Modulation = "64QAM";
    end

    % Target Code rate - Row#5 + Note#2:
    targetRatesPerSec = [308 658 193 658 567]/1024; % Annex A, TS 38.104
    tdef.TargetCodeRate = targetRatesPerSec(secNum);

    % All tables specify 1 Layer, except for a few that specify 2:
    if (any(secNum == [3 4]) && frNum == 1 && (frcNum >= 15 && frcNum <= 28)) || ...
       (any(secNum == [3 4]) && frNum == 2 && (frcNum >= 6  && frcNum <= 10))
      tdef.NLayers = 2;
    else
      tdef.NLayers = 1;
    end

    % All tables disable TransformPrecoding, except for:
    tdef.TransformPrecoding = (secNum == 3 && frNum == 1 && any(frcNum == [29 30 31 32])) || ...
                                  (secNum == 3 && frNum == 2 && any(frcNum == [11 12])); 

    % full allocation in time domain, except for some specific tables with 9 data symbols per slot (not counting DM-RS containing symbol)
    if frNum == 2 && secNum >= 3
      tdef.NAllocatedSymbols = 10;
    else
      tdef.NAllocatedSymbols = 14;
    end                          

    % Default to non-slotwise type A
    tdef.MappingType = 'A';
    tdef.DMRSTypeAPosition = 2;
     
    %  First DM-RS symbol position (l_0)
    if frNum == 2 && secNum >= 3
      tdef.MappingType = 'B';
    end

    % Additional positions
    if (frNum == 1 && any(secNum == [3 4 5]) && any(frcNum == [1:7 15:21 29 30]))
      tdef.DMRSAdditionalPosition = 0;
    elseif (frNum == 2 && secNum > 1)
      tdef.DMRSAdditionalPosition = 0;
    else
      tdef.DMRSAdditionalPosition = 1; % all other cases
    end
    if frNum==2 && secNum==5 && any(frcNum == 6:10)
        tdef.DMRSAdditionalPosition = 1;
    end   

    % Establish the active SCS
    activescs = tdef.SubcarrierSpacing;
    if nargin > 2 && ~isempty(scs)  % If empty then use the FRC SCS
        activescs = scs;
    end
    [looknrb,lookbw,lookscs] = getValidNRB(tdef.FR,bw,activescs);
    minnrbidx = find(looknrb >= tdef.AllocatedNRB,1);   % Find the smallest BW 
    if isempty(minnrbidx)
        minnrbidx = find(~isnan(looknrb),1,'last');  
        newnrb = looknrb(minnrbidx);
        % This could occur if the BW has been specified but is too small for the allocated NRB, or the SCS has been specified 
        % but allocated NRB exceeds the possible BW for that SCS
        warning('The allocated PUSCH NRB was changed from %d to %d so that it would fit into the specified channel bandwidth (%dMHz) and subcarrier spacing (%dkHz)', tdef.AllocatedNRB,newnrb,lookbw(minnrbidx),lookscs);
        tdef.AllocatedNRB = newnrb;
    end
    tdef.ChannelBandwidth = lookbw(minnrbidx);      % Channel BW (MHz)
    tdef.ChannelNRB = looknrb(minnrbidx);           % Number of RB in carrier in channel (max RB for that channel BW and SCS)
    tdef.SubcarrierSpacing = lookscs;
    
    % Default to TDD for FR2
    duplexingmodes = ["FDD","TDD"];
    tdef.DuplexMode = duplexingmodes(1+double(tdef.FR=="FR2"));
end


function scs = getSCSUplinkFRC(frNum, secNum, frcNum)
    % Annex A, TS 38.104.
    % 2 cell arrays, one for each frequency range (FR).
    % each cell is indexed by section number, e.g., 3 for 'G-FR1-A3-18' to get
    % a vector of SCS that is indexed by the FRC for this section, e.g., 18 for 'G-FR1-A3-18'
    commonFR1 = [15 15 15 30 30 30 30];
    fr1SCS = {[15 30 60 15 30 60 15 30 60], ...         % A1
              [15 30 60 15 30 60], ...                  % A2
              [repmat(commonFR1, 1, 4) 15 30 15 30], ...% A3
               repmat(commonFR1, 1, 4), ...             % A4
               repmat(commonFR1, 1, 2) };               % A5

    commonFR2 = [60 60 120 120 120];
    fr2SCS = {[60 120 120 60 120], ...              % A1
              [], ...                               % A2
              [repmat(commonFR2, 1, 2) 60 120], ... % A3
               repmat(commonFR2, 1, 2), ...         % A4
               repmat(commonFR2, 1, 2)};            % A5

    if nargin == 2
        scs = fr1SCS{secNum};
        return;
    end

    % indexing:
    if frNum == 1
      scs = fr1SCS{secNum}(frcNum);
    else % FR2
      scs = fr2SCS{secNum}(frcNum);
    end

end

function nrb = getNRBUplinkFRC(frNum, secNum, frcNum)
    % Annex A, TS 38.104.
    % 2 cell arrays, one for each frequency range (FR).
    % Each cell is indexed by section number, e.g., 3 for 'G-FR1-A3-18' to get
    % a vector of NRB that is indexed by the FRC for this section, e.g., 18 for 'G-FR1-A3-18'
    commonFR1 = [25 52 106 24 51 106 273];
    fr1NRB = {[25 11 11 106 51 24 15 6 6], ...              % A1
              [25 11 11 106 51 24], ...                     % A2
              [repmat(commonFR1, 1, 4) 25 24 25 24] , ...   % A3
               repmat(commonFR1, 1, 4), ...                 % A4
               repmat(commonFR1, 1, 2)};                    % A5

    commonFR2 = [66 132 32 66 132];  % Represents the standard [60kHz/66RB, 60kHz/132RB, 120kHz/32RB, 120kHz/66RB, 120kHz/132RB] columns
    fr2NRB = {[66 32 66 33 16], ...                         % A1
              [], ...                                       % A2
              [repmat(commonFR2, 1, 2) 30 30], ...          % A3
               repmat(commonFR2, 1, 2), ...                 % A4
               repmat(commonFR2, 1 ,2)};                    % A5

    % If there were no indexing inputs to the function
    % then return the entire table for the FR
    if nargin == 1
        if frNum == 1
            nrb = fr1NRB;
        else % FR2
            nrb = fr2NRB;
        end
        return;
    end   

    % Other use indexing inputs to lookup tables
    if frNum == 1
      nrb = fr1NRB{secNum}(frcNum);
    else % FR2
      nrb = fr2NRB{secNum}(frcNum);
    end

end

function waveconfig = getFRCPUSCHParameters(tdef,ncellid)

    tdd = strcmpi(tdef.DuplexMode,'TDD');   % Local bool for TDD or FDD
   
    % Carrier-level parameters
    waveconfig.Name = tdef.Name;         % Label the configuration
    waveconfig.NCellID = ncellid;        % Cell identity
    waveconfig.ChannelBandwidth = tdef.ChannelBandwidth;    % Channel bandwidth (MHz)
    waveconfig.FrequencyRange = tdef.FR; % 'FR1' or 'FR2'
    waveconfig.NumSubframes = 10+10*tdd; % Number of 1ms subframes in generated waveform (1,2,4,8 slots per 1ms subframe, depending on SCS, and 10 subframes in a 10ms frame)   
    waveconfig.Windowing = 0;            % Windowing amount (set it to 0 in these presets i.e. no windowing)
    waveconfig.DisplayGrids = 0;         % Display the resource grids after signal generation

    % SCS carrier configuration
    % Define a single SCS-specific carrier for the waveform
    carriers(1).SubcarrierSpacing = tdef.SubcarrierSpacing;
    carriers(1).NRB = tdef.ChannelNRB;
    carriers(1).RBStart = 0;

    % Bandwidth part configuration
    % Define a single BWP with spans the whole SCS carrier
    bwp(1).SubcarrierSpacing = carriers(1).SubcarrierSpacing; % BWP Subcarrier Spacing
    bwp(1).CyclicPrefix = 'Normal';                           % BWP Cyclic prefix
    bwp(1).NRB = carriers(1).NRB;                             % Size of BWP
    bwp(1).RBOffset = 0;                                      % Position of BWP in SCS carrier
     
    %  PUCCH - Physical uplink control channel configuration. A structure array with the following fields:                   
    pucch(1).Enable = 0;     % Enable/Disable the PUCCH instance. It is either 0 or 1                                                       

    % PUSCH - Physical uplink shared channel configuration. A structure array with the following fields:
    pusch(1).Name = sprintf("PUSCH sequence for %s",tdef.Name);% Description of this PUSCH sequence           
    pusch(1).Enable  = 1;                  % Enable/Disable the PUSCH configuration. It is either 0 or 1                                
    pusch(1).BWP = 1;                      % Bandwidth part identity. It is indicated to which bandwidth part PUSCH is configured                         
    pusch(1).Power = 0;                    % Power scaling (in dB)
    pusch(1).EnableCoding = 1;             % Enable UL-SCH transport channel coding
    pusch(1).TargetCodeRate = tdef.TargetCodeRate; % Code rate used to calculate transport block sizes                             
    pusch(1).Xoh_PUSCH = 0;                % Overhead parameter. It is one of the set {0,6,12,18}                              
    pusch(1).TxScheme = 'codebook';        % Transmission scheme. It is one of the set {'codebook','nonCodebook'}. When set to codebook, the precoding matrix is applied based on the number of layers and the antenna ports                        
    pusch(1).Modulation = tdef.Modulation; % Modulation scheme. It is one of the set {'pi/2-BPSK','QPSK','16QAM','64QAM','256QAM'}                                
    pusch(1).NLayers = tdef.NLayers;       % Number of layers. It is in range 1 to 4
    pusch(1).NAntennaPorts = pusch(1).NLayers; % Number of antenna ports. It is one of the set {1,2,4}                                
    pusch(1).RVSequence = 0;               % Redundancy version sequence. It is applied cyclically across the PUSCH allocation sequence                             
    pusch(1).IntraSlotFreqHopping = 'disabled'; % Intra-slot frequency hopping configuration. It is one of the set {'enabled',disabled'}                             
    pusch(1).TransformPrecoding = tdef.TransformPrecoding;  % Transform precoding flag. It is either 0 or 1. When set to 1, DFT operation is applied before precoding and it is only used for single layer                              
    pusch(1).TPMI = 0;                     % Transmitted precoding matrix indicator. The range depends on the number of layers and antenna ports                             
    pusch(1).GroupHopping = 'neither';     % Group hopping configuration. It is one of the set {'enable','disable','neither'}                                     
    pusch(1).RBOffset = 0;                 % Resource block offset for second hop
    pusch(1).InterSlotFreqHopping = 'disabled';  % Inter-slot frequency hopping configuration. It is one of the set{'enabled',disabled'}. When the value is set to 'enabled', intra-slot frequency hopping is considered as disabled                               
    pusch(1).NID = [];                     % Scrambling identity. It is in range 0 to 1023. Use empty ([]) to use physical layer cell identity                                                              
    pusch(1).RNTI = 1;                     % Radio network temporary identifier
    pusch(1).DataSource = "PN9-ITU";       % Transport block data source. Use one of the following standard PN sequences: 'PN9-ITU', 'PN9', 'PN11', 'PN15', 'PN23'.                              
                                           % The seed for the generator can be specified using a cell array in the form {'PN9',seed}. If no seed is specified, the generator is initialized with all ones                          
    pusch(1).PUSCHMappingType = tdef.MappingType; % PUSCH mapping type. It is either 'A' or 'B'.     
    pusch(1).AllocatedSymbols = 0:tdef.NAllocatedSymbols-1; % Symbols in a slot. It needs to be a contiguous allocation. For PUSCH mapping type 'A', the start symbol must be zero and the length can be from 4 to 14 (normal CP) and up to 12 (extended CP)                             
    pusch(1).AllocatedSlots = 0;           % Slots in a frame used for PUSCH
    pusch(1).AllocatedPeriod = 1;          % Allocation period in slots. Use empty for no repetition  
    allocprboffset = fix((tdef.ChannelNRB - tdef.AllocatedNRB)/2);     % RB offset to place allocation in the middle of the band
    pusch(1).AllocatedPRB = allocprboffset+(0:tdef.AllocatedNRB-1);    % PRB allocation
    pusch(1).PortSet = 0:tdef.NLayers-1;   % DM-RS ports to use for the layers
    pusch(1).DMRSTypeAPosition = tdef.DMRSTypeAPosition;% DM-RS symbol position for mapping type 'A'. It is either 2 or 3            
    pusch(1).DMRSLength = 1;               % Number of front-loaded DM-RS symbols. It is either 1 (single symbol) or 2 (double symbol)                    
    pusch(1).DMRSAdditionalPosition = tdef.DMRSAdditionalPosition; % Additional DM-RS symbols positions. It is in range 0 to 3. Value of zero indicates no additional DM-RS symbols                       
    pusch(1).DMRSConfigurationType = 1;    % DM-RS configuration type. It is either 1 or 2. The number of subcarriers
    pusch(1).NumCDMGroupsWithoutData = 2;  % Number of DM-RS CDM groups without data. It is value in range 1 to 3                          
    pusch(1).NIDNSCID = 0;                 % Scrambling identity for CP-OFDM DMRS. It is in range 0 to 65535. Use empty ([]) to use physical layer cell identity NCellID                     
    pusch(1).NSCID = 0;                    % CP-OFDM DM-RS scrambling initialisation. It is either 0 or 1
    pusch(1).NRSID = [];                   % Scrambling identity for DFT-s-OFDM DM-RS. It is in range 0 to 1007. Use empty ([]) to use the physical layer cell identity NCellID                                                
    pusch(1).PowerDMRS = 3;                % Power boosting for DM-RS (in dB)
    pusch(1).DisableULSCH = 0;             % Disable UL-SCH on overlapping slots of PUSCH and PUCCH. It is either 0 or 1
    pusch(1).BetaOffsetACK = 0;            % Rate matching offset for HARQ-ACK
    pusch(1).BetaOffsetCSI1 = 0;           % Rate matching offset for CSI part 1
    pusch(1).BetaOffsetCSI2 = 0;           % Rate matching offset for CSI part 2
    pusch(1).ScalingFactor = 1;            % Scaling factor
    pusch(1).EnablePTRS = (tdef.FR=="FR2") && (~tdef.TransformPrecoding);         % Enable or disable the PT-RS (1 or 0)
    pusch(1).PTRSFrequencyDensity = 2;     % Frequency density of PT-RS (2,4)
    pusch(1).PTRSTimeDensity = 1;          % Time density of PT-RS (1,2,4)
    pusch(1).PTRSNumSamples = 2;           % Number of PT-RS samples (2,4). It is used only for DFT-s-OFDM
    pusch(1).PTRSNumGroups = 2;            % Number of PT-RS groups (2,4,8). It is used only for DFT-s-OFDM                  
    pusch(1).PTRSREOffset = '00';          % Resource element offset ('00','01','10','11')
    pusch(1).PTRSPortSet  = 0:tdef.NLayers-1; % Antenna ports of PT-RS
    pusch(1).PTRSNID = [];                 % Scrambling identity of PT-RS (0...1007). It is used only for DFT-s-OFDM                 
    pusch(1).PowerPTRS = 0;                % Power boosting for PT-RS (dB)
    
    % Finalize the PUSCH slot allocation and period part 
    pusch = setFRCSlotAllocation(pusch,tdef,tdef.SubcarrierSpacing);
    
    % Putting everything together
    waveconfig.Carriers = carriers;
    waveconfig.BWP = bwp;
    waveconfig.PUCCH = pucch;
    waveconfig.PUSCH = pusch;
      
end

% Plot the PRB level resource grids
function plotResourceGrid(ax, ~, cfgObj, gridset, isDownlink)

    if nargin < 5
        isDownlink = 1;
    end

    wirelessWaveformGenerator.internal.plotResourceGrid(ax, cfgObj, gridset, isDownlink, 1);

    if isDownlink
        channelList = 'PDSCH and PDCCH';
    else
        channelList = 'PUSCH';
    end
    title(ax, sprintf('BWP 1 in Carrier (SCS=%dkHz). %s location',cfgObj.BandwidthParts{1}.SubcarrierSpacing, channelList));

end

% Get the generator configuration in object form
function cfgObj = getGenConfigObj(gen)

    cfgObj = gen.ConfigValue;
    if isstruct(cfgObj)
        if gen.LinkDirection == "uplink"
            cfgObj = mapStruct2ObjUplink(cfgObj);
        else
            cfgObj = mapStruct2ObjDownlink(cfgObj);
        end
    end
end            

function cfgObj = mapStruct2ObjUplink(cfgS)

    % Map the FRC configurations into object-based ones needed by nrWaveformGenerator
    cfgObj = nrULCarrierConfig;

    % Top-level parameters
    [cfgObj,cfgS] = mapStruct2ObjCCCommon(cfgObj,cfgS);
    
    % Continue with uplink channels and resources 
    cfgObj.PUSCH                 = mapPUSCH(cfgS);

end


function cfgObj = mapStruct2ObjDownlink(cfgS)

    % Map the FRC and TM preset configurations into object-based ones needed by nrWaveformGenerator
    cfgObj = nrDLCarrierConfig;

    % Top-level parameters
    [cfgObj,cfgS] = mapStruct2ObjCCCommon(cfgObj,cfgS);
    
    % Continue with downlink channels and resources 
    cfgObj.SSBurst               = mapSSB(cfgS);
    [cset, ss, pdcch]            = mapControl(cfgS);
    cfgObj.CORESET               = cset;
    cfgObj.SearchSpaces          = ss;
    cfgObj.PDCCH                 = pdcch;
    cfgObj.PDSCH                 = mapPDSCH(cfgS);
    cfgObj.CSIRS                 = mapCSIRS(cfgS);
end

% Component carrier common parameters
function [cfgObj,cfgS] = mapStruct2ObjCCCommon(cfgObj,cfgS)

    % Top-level parameters
    cfgObj.FrequencyRange        = cfgS.FrequencyRange;
    cfgObj.ChannelBandwidth      = cfgS.ChannelBandwidth;
    cfgObj.NCellID               = cfgS.NCellID;
    cfgObj.NumSubframes          = cfgS.NumSubframes;
    if isprop(cfgObj,'WindowingPercent') && isfield(cfgS,'Windowing')
        % Only map to either off or default windowing levels
        if cfgS.Windowing == 0
            w = 0;
        else
            w = [];
        end
        cfgObj.WindowingPercent = w;
    end   
    % Add in a name field should one not exist and copy across string otherwise
    if ~isfield(cfgS,'Name')
        cfgS.Name = '';
    else
        cfgObj.Label = cfgS.Name;
    end
    
    cfgObj.SCSCarriers           = mapCarriers(cfgS);
    [cfgObj.BandwidthParts,cfgS] = mapBWP(cfgS);

end


function scsCfg = mapCarriers(cfgS)
    
    % SCS carrier configuration
    carriers = cfgS.Carriers;
    for idx = 1:length(carriers)
        scsCfg{idx} = nrSCSCarrierConfig;
        % Map across parameters
        scsCfg{idx}.SubcarrierSpacing = carriers(idx).SubcarrierSpacing; %#ok<*AGROW>
        scsCfg{idx}.NSizeGrid         = carriers(idx).NRB;
        scsCfg{idx}.NStartGrid        = carriers(idx).RBStart;
    end
    
    % React to any SS burst parameters present
    if isfield(cfgS,'SSBurst')
        % Create a carrier for the SS Burst, if a relevant one does not yet exist
        if cfgS.SSBurst.Enable 
            burstSCS = nr5g.internal.wavegen.blockPattern2SCS(cfgS.SSBurst.BlockPattern);
            scs = [carriers.SubcarrierSpacing];
            if ~any(burstSCS == scs)

                [~,imax] = max(scs);
                % Center frequencies of lowest and highest subcarriers for highest SCS carrier (relative to point A as 0 Hz)
                fMin = scsCfg{imax}.NStartGrid*12*scsCfg{imax}.SubcarrierSpacing*1e3;
                fMax = fMin + scsCfg{imax}.NSizeGrid*12*scsCfg{imax}.SubcarrierSpacing*1e3;
                waveformCenter = fMin+(fMax-fMin)/2;

                idx = idx + 1;
                scsCfg{idx} = nrSCSCarrierConfig;         % Add a new SCS carrier to the end of the set
                scsCfg{idx}.SubcarrierSpacing = burstSCS; % #ok<*AGROW>
                scsCfg{idx}.NSizeGrid = 20;               % Size it for the SS burst alone
                % Place this carrier at the center of the frequencies spanned by the other set of carriers
                ssbBottom = waveformCenter - (scsCfg{idx}.NSizeGrid/2)*12*burstSCS*1e3;
                start = ssbBottom/(burstSCS*12*1e3);
                if floor(start)~=start
                    scsCfg{idx}.NSizeGrid = scsCfg{idx}.NSizeGrid+1;
                    ssbBottom = waveformCenter - (scsCfg{idx}.NSizeGrid/2)*12*burstSCS*1e3;
                    start = ssbBottom/(burstSCS*12*1e3);
                end
                scsCfg{idx}.NStartGrid = start;
            end
        end
    end
   
end

function [bwp2,cfgS] = mapBWP(cfgS)
  
    bwp2 = {};

    bwp = cfgS.BWP;
    for idx = 1:length(bwp)
        bwp2{idx} = nrWavegenBWPConfig;
        bwp2{idx}.BandwidthPartID = idx;
        bwp2{idx}.SubcarrierSpacing = bwp(idx).SubcarrierSpacing;
        bwp2{idx}.CyclicPrefix      = bwp(idx).CyclicPrefix;
        bwp2{idx}.NSizeBWP          = bwp(idx).NRB;
        % Find carrier for corresponding bwp 
        carrierid = [cfgS.Carriers(:).SubcarrierSpacing] == bwp2{idx}.SubcarrierSpacing;
        bwp2{idx}.NStartBWP         = bwp(idx).RBOffset + cfgS.Carriers(carrierid).RBStart;
        cfgS.BWP(idx).Carrier = find(carrierid,1);     % Link back the SCS carrier index associated with the BWP (linked through SCS value)
    end
end

function ssb2 = mapSSB(cfgS)

    ssb = cfgS.SSBurst;

    ssb2 = nrWavegenSSBurstConfig;
    ssb2.Enable             = ssb.Enable;
    ssb2.Power              = ssb.Power;
    ssb2.BlockPattern       = char(ssb.BlockPattern);
    ssb2.TransmittedBlocks  = ssb.SSBTransmitted;
    ssb2.Period             = ssb.SSBPeriodicity;
    ssb2.NCRBSSB            = []; % Carrier center
    ssb2.DataSource         = 'MIB';
    % SubcarrierSpacingCommon must be 15 or 30 for FR1 (Case A, B, C) and 60 or 120 kHz for FR2 (Case D, E)
    ssb2.SubcarrierSpacingCommon = 30*(1+strcmpi(cfgS.FrequencyRange,'FR2'));
    
end

function [coreset2, ss2, pdcch2] = mapControl(cfgS)
    
    coreset2 = {};
    ss2 = {};
    pdcch2 = {};

    % The structure based parameterization links CORESET and BWP within the PDCCH definitions
    coresetused = [cfgS.PDCCH(:).CORESET];
    bwpused = [cfgS.PDCCH(:).BWP];
    
    % CORESET & search spaces
    coreset = cfgS.CORESET;
    for idx = 1:length(coreset)

        % Each CORESET structure instance will be mapped into a CORESET and search space object pair.
        % The original CORESET/search space parameters were combined into a single 
        % container and need to be separated out now as follows:
        %
        %             Duration: 2                         -> CORESET
        %         AllocatedPRB: [0 1 2 3 4 5 6 7... 24]   -> CORESET
        %        CCEREGMapping: 'noninterleaved'          -> CORESET
        %        REGBundleSize: 2                         -> CORESET
        %      InterleaverSize: 2                         -> CORESET
        %           ShiftIndex: 0                         -> CORESET
        %
        %     AllocatedSymbols: 0                         -> Search space
        %       AllocatedSlots: [1 2 3 4 5 6 7 8 9]       -> Search space
        %      AllocatedPeriod: 10                        -> Search space

        % Map the CORESET part
        coreset2{idx} = nrCORESETConfig;
        coreset2{idx}.CORESETID = idx;
        coreset2{idx}.Label = ['CORESET' num2str(idx)];

        % For the CORESET frequency resource bitmap, each bit is a block of 6 *CRB* in the associated BWP, 
        % with the first bit being associated with first complete block of 6 CRB in the BWP, and any
        % incomplete block of 6 CRB is associated to have the bits set to 0 (RRC).
        % So, when turning a set of PRB in a bitmap, we actually need to know their CRB, and thus the
        % associated BWP containing the CORESET.
        %
        % The original CORESET definition used a set of PRB to define the frequency resources. This 
        % set was allowed to be outside of the blocks of 6 CRB, for ease of specification, and, in use,
        % would only select the legal blocks as part of the actual CORESET. We now need to map for the behavior.
        % 
        % Identify all BWP associated with this CORESET ('idx'), originally defined in terms of PRB,
        % and find the sizing info associated with projecting the CORESET PRB into each BWP
              
        % Find all BWP that are associated with this CORESET   
        bwpusingcoreset = bwpused(coresetused==idx);
        allblocks = []; maxsizes = []; maxnblocks = [];
        for b = bwpusingcoreset
            [coresetprb,nblocks,maxnblocks] = getCORESETPRB(cfgS.Carriers(cfgS.BWP(b).Carrier),cfgS.BWP(b),coreset(idx)); %#ok<ASGLU>
            allblocks = [allblocks nblocks];
            maxsizes = [maxsizes maxnblocks];
        end
        % Create a frequency bitmap which will allow the CORESET to fit into all references BWP
        bmap(1+allblocks( allblocks < min(maxnblocks))) = 1;
        coreset2{idx}.FrequencyResources = bmap;

        % Copy other CORESET parameters across
        coreset2{idx}.Duration = coreset(idx).Duration;  
        if isfield(coreset(idx), 'CCEREGMapping')
          coreset2{idx}.CCEREGMapping   = coreset(idx).CCEREGMapping;
          coreset2{idx}.REGBundleSize   = coreset(idx).REGBundleSize;
          coreset2{idx}.InterleaverSize = coreset(idx).InterleaverSize;
          coreset2{idx}.ShiftIndex      = coreset(idx).ShiftIndex;
        end
        
        % Map the search space part
        ss2{idx} = nrSearchSpaceConfig;
        ss2{idx}.SearchSpaceID = idx;
        ss2{idx}.Label = ['SearchSpace' num2str(ss2{idx}.SearchSpaceID)];
        ss2{idx}.CORESETID = idx;

        ss2{idx}.SearchSpaceType = 'common';    % 'common' or 'ue' (implications for scrambling and candidate placement control) 
        
        % The search space can only identify itself with one symbol starting position in each slot occasion,
        % so, for an input CORESET defining locations at multiple symbols within a each slot occasion,
        % separate search space components would need to be defined. Here only the first occasion is mapped.
        % 
        % Define the search space to occur in every slot and let the PDCCH sequences
        % specify the slots used for their transmission instances
        symbols = coreset(idx).AllocatedSymbols;        % Get first symbols of instances in an active slot
        ss2{idx}.StartSymbolWithinSlot = symbols(1);    % First symbol used in each slot
        ss2{idx}.Duration = 1;                          % Number of slots in each 'block' slots (starting at 'offset' slot)
        ss2{idx}.SlotPeriodAndOffset = [1 0];           % Slot pattern associated with SS pattern (slot period and slot offset associated with a block of 'Duration' slot)
        
        % Avoid validation error occuring if an AL extends beyond the CORESET
        numREGs = 6*numel(bmap)*coreset(idx).Duration;
        crstCCEs = fix(numREGs/6);              % Number of CCE in associated CORESET (6 REG per CCE) 
        ALs = [1 2 4 8 16];                     % Possible aggregation levels 
        ind2del = ALs>crstCCEs;                 % AL which don't fit in the complete CORESET
        ss2{idx}.NumCandidates(ind2del) = 0;

    end
    % End of CORESET/Search space for loop
  
    pdcch = cfgS.PDCCH;
    for idx = 1:length(pdcch)
        pdcch2{idx} = nrWavegenPDCCHConfig;
        pdcch2{idx}.Enable                = pdcch(idx).Enable;
        pdcch2{idx}.Label                 = ['PDCCH' num2str(idx)];
        pdcch2{idx}.Power                 = pdcch(idx).Power;
        pdcch2{idx}.BandwidthPartID       = pdcch(idx).BWP;
        pdcch2{idx}.Coding                = pdcch(idx).EnableCoding;
    
        pdcch2{idx}.SearchSpaceID =  pdcch(idx).CORESET;  %  ID of search space containing this PDCCH
        pdcch2{idx}.AggregationLevel = pdcch(idx).NumCCE; %  PDCCH aggregation level {1,2,4,8,16}
        pdcch2{idx}.AllocatedCandidate = 1 + fix(pdcch(idx).StartCCE/pdcch(idx).NumCCE); % Candidate used for the PDCCH instance
 
        % Original CORESET/search space allocation used by the PDCCH sequence
        coresetidx = pdcch(idx).CORESET;
        sallocation = coreset(coresetidx).AllocatedSlots;
        speriod = coreset(coresetidx).AllocatedPeriod;

        % Original PDCCH allocation
        callocation = pdcch(idx).AllocatedSearchSpaces;
        cperiod = pdcch(idx).AllocatedPeriod;

        % Expand to get the actual slot sequence used
        [eallocation,eperiod] = controlexpansion(sallocation,speriod,callocation,cperiod);

        pdcch2{idx}.SlotAllocation = eallocation;    % Time-domain location of PDCCH (in slots)
        pdcch2{idx}.Period = eperiod;                % Slot allocation period of PDCCH 
    
        % Map additional parameters
        pdcch2{idx}.RNTI                  = pdcch(idx).RNTI;           % Used in the PDCCH data scrambling (NRNTI), and UE RNTI (nrPDCCHResources)
        pdcch2{idx}.DMRSScramblingID      = pdcch(idx).NID;            % Used in the PDCCH data scrambling (NID), and DM-RS symbols (nrPDCCHResources)
        pdcch2{idx}.DMRSPower             = pdcch(idx).PowerDMRS;
        pdcch2{idx}.DataBlockSize         = pdcch(idx).DataBlkSize;
        pdcch2{idx}.DataSource            = pdcch(idx).DataSource;
    end
  
end

function pdsch2 = mapPDSCH(cfgS)

    pdsch2 = {};

    % PDSCH
    pdsch = cfgS.PDSCH;
    for idx = 1:length(pdsch)
        pdsch2{idx} = nrWavegenPDSCHConfig;
        pdsch2{idx}.Enable               = pdsch(idx).Enable;
        pdsch2{idx}.Label                = pdsch(idx).Name;
        pdsch2{idx}.Power                = pdsch(idx).Power;
        pdsch2{idx}.BandwidthPartID      = pdsch(idx).BWP;
        pdsch2{idx}.Coding               = pdsch(idx).EnableCoding;
        pdsch2{idx}.DataSource           = pdsch(idx).DataSource;
        pdsch2{idx}.TargetCodeRate       = pdsch(idx).TargetCodeRate;
        pdsch2{idx}.TBScaling            = 1;
        pdsch2{idx}.XOverhead            = pdsch(idx).Xoh_PDSCH;
        pdsch2{idx}.Modulation           = pdsch(idx).Modulation;
        pdsch2{idx}.NumLayers            = pdsch(idx).NLayers;
        pdsch2{idx}.RVSequence           = pdsch(idx).RVSequence;
        pdsch2{idx}.VRBToPRBInterleaving = pdsch(idx).VRBToPRBInterleaving;
        pdsch2{idx}.VRBBundleSize        = pdsch(idx).VRBBundleSize;
        pdsch2{idx}.SymbolAllocation     = [pdsch(idx).AllocatedSymbols(1) (pdsch(idx).AllocatedSymbols(end)+1-pdsch(idx).AllocatedSymbols(1))];
        pdsch2{idx}.SlotAllocation       = pdsch(idx).AllocatedSlots;
        pdsch2{idx}.Period               = pdsch(idx).AllocatedPeriod;
        pdsch2{idx}.PRBSet               = pdsch(idx).AllocatedPRB;
        pdsch2{idx}.RNTI                 = pdsch(idx).RNTI;
        pdsch2{idx}.NID                  = pdsch(idx).NID;
        pdsch2{idx}.MappingType          = pdsch(idx).PDSCHMappingType;

        for idx2 = 1:length(pdsch(idx).RateMatch)
            pdsch2{idx}.ReservedCORESET     = [pdsch2{idx}.ReservedCORESET pdsch(idx).RateMatch(idx2).CORESET];
            pdsch2{idx}.ReservedPRB{idx2} = nrPDSCHReservedConfig;     
            pdsch2{idx}.ReservedPRB{idx2}.PRBSet    = pdsch(idx).RateMatch(idx2).Pattern.AllocatedPRB;
            pdsch2{idx}.ReservedPRB{idx2}.SymbolSet = pdsch(idx).RateMatch(idx2).Pattern.AllocatedSymbols;
            pdsch2{idx}.ReservedPRB{idx2}.Period    = pdsch(idx).RateMatch(idx2).Pattern.AllocatedPeriod;
        end

        % Antenna port and DM-RS configuration (TS 38.211 section 7.4.1.1)
        % Using DM-RS reference point CRB0
        dmrs{idx} = nrPDSCHDMRSConfig;
        dmrs{idx}.DMRSPortSet             = pdsch(idx).PortSet;
        dmrs{idx}.DMRSTypeAPosition       = pdsch(idx).DMRSTypeAPosition;
        dmrs{idx}.DMRSLength              = pdsch(idx).DMRSLength;
        dmrs{idx}.DMRSAdditionalPosition  = pdsch(idx).DMRSAdditionalPosition;
        dmrs{idx}.DMRSConfigurationType   = pdsch(idx).DMRSConfigurationType;
        dmrs{idx}.NumCDMGroupsWithoutData = pdsch(idx).NumCDMGroupsWithoutData;
        dmrs{idx}.NIDNSCID                = pdsch(idx).NIDNSCID;
        dmrs{idx}.NSCID                   = pdsch(idx).NSCID;
        pdsch2{idx}.DMRSPower             = pdsch(idx).PowerDMRS;
        pdsch2{idx}.DMRS = dmrs{idx};

        ptrs{idx} = nrPDSCHPTRSConfig;
        pdsch2{idx}.EnablePTRS              = pdsch(idx).EnablePTRS;
        ptrs{idx}.TimeDensity               = pdsch(idx).PTRSTimeDensity;
        ptrs{idx}.FrequencyDensity          = pdsch(idx).PTRSFrequencyDensity;
        ptrs{idx}.REOffset                  = pdsch(idx).PTRSREOffset;
        ptrs{idx}.PTRSPortSet               = pdsch(idx).PTRSPortSet;
        pdsch2{idx}.PTRSPower               = pdsch(idx).PowerPTRS;
        pdsch2{idx}.PTRS = ptrs{idx};
    end
end

function [o,cond] = passign(s,o,f,p,ac)
% Assign from s.f->o.f or s.f->o.p conditionally
% Structure, object, fieldname, propertyname, additionalcondition

    cond = isfield(s,f) && ~isempty(s.(f));  % Base assignment condition is that 'f' field be present and its value is non-empty
    if nargin == 5
        cond = cond && ac;                   % Combine with additional input condition if present
    end

    if cond
        if nargin == 3
            o.(f) = s.(f);     % Direct mapping with same source/target name
        else
            o.(p) = s.(f);     % Mapping with change of source/target name
        end
    end
end

function pusch2 = mapPUSCH(cfgS)

    pusch2 = {};

    % PDSCH
    pusch = cfgS.PUSCH;
    for idx = 1:length(pusch)
        pusch2{idx} = nrWavegenPUSCHConfig;
        % General shared parameters
        pusch2{idx}.Enable               = pusch(idx).Enable;
        pusch2{idx}.Label                = pusch(idx).Name;
        pusch2{idx}.Power                = pusch(idx).Power;
        pusch2{idx}.BandwidthPartID      = pusch(idx).BWP;
        pusch2{idx}.Coding               = pusch(idx).EnableCoding;
        pusch2{idx}.DataSource           = pusch(idx).DataSource;
        pusch2{idx}.TargetCodeRate       = pusch(idx).TargetCodeRate;
        pusch2{idx}.XOverhead            = pusch(idx).Xoh_PUSCH;
        pusch2{idx}.Modulation           = pusch(idx).Modulation;
        pusch2{idx}.NumLayers            = pusch(idx).NLayers;
        pusch2{idx}.RVSequence           = pusch(idx).RVSequence;
        pusch2{idx}.SymbolAllocation     = [pusch(idx).AllocatedSymbols(1) (pusch(idx).AllocatedSymbols(end)+1-pusch(idx).AllocatedSymbols(1))];
        pusch2{idx}.SlotAllocation       = pusch(idx).AllocatedSlots;
        pusch2{idx}.Period               = pusch(idx).AllocatedPeriod;
        pusch2{idx}.PRBSet               = pusch(idx).AllocatedPRB;
        pusch2{idx}.RNTI                 = pusch(idx).RNTI;
        pusch2{idx}.NID                  = pusch(idx).NID;
        pusch2{idx}.MappingType          = pusch(idx).PUSCHMappingType;
  
        % PUSCH specific parameters
        pusch2{idx}.TransformPrecoding   = pusch(idx).TransformPrecoding;
        pusch2{idx}.TransmissionScheme   = pusch(idx).TxScheme;
        pusch2{idx}.NumAntennaPorts      = pusch(idx).NAntennaPorts;
        pusch2{idx}.TPMI                 = pusch(idx).TPMI;
        if strcmpi(pusch(idx).InterSlotFreqHopping,'enabled')
            pusch2{idx}.FrequencyHopping = 'interSlot';
        elseif strcmpi(pusch(idx).IntraSlotFreqHopping,'enabled')
            pusch2{idx}.FrequencyHopping = 'intraSlot';
        else
            pusch2{idx}.FrequencyHopping = 'neither';
        end
        if ~strcmpi(pusch2{idx}.FrequencyHopping,'neither')
            if isempty(pusch(idx).RBOffset) % Replace the empty value with 0
                pusch(idx).RBOffset = 0;
            end
            secondHopPRB = mod(pusch(idx).AllocatedPRB(1)+pusch(idx).RBOffset,cfgS.BWP(idx).NRB);
            pusch(idx).RBOffset = secondHopPRB; % Replace RB offset with the calculated secondHopPRB
        end
        pusch2{idx}.SecondHopStartPRB    = pusch(idx).RBOffset;

        % UCI part
        pusch2{idx}.EnableULSCH          = ~pusch(idx).DisableULSCH;
        pusch2{idx} = passign(pusch(idx),pusch2{idx},'BetaOffsetACK','BetaOffsetACK',pusch(idx).BetaOffsetACK);
        pusch2{idx} = passign(pusch(idx),pusch2{idx},'BetaOffsetCSI1','BetaOffsetCSI1',pusch(idx).BetaOffsetCSI1);
        pusch2{idx} = passign(pusch(idx),pusch2{idx},'BetaOffsetCSI2','BetaOffsetCSI2',pusch(idx).BetaOffsetCSI2);
        pusch2{idx}.UCIScaling           = pusch(idx).ScalingFactor;
        % Don't bother mapping any UCI enabling here
        pusch2{idx}.EnableACK            = 0; 
        pusch2{idx}.EnableCSI1           = 0; 
        pusch2{idx}.EnableCSI2           = 0;
        pusch2{idx}.EnableCGUCI          = 0;
        
        % Antenna port and DM-RS configuration (TS 38.211 section 7.4.1.1)
        % Using DM-RS reference point CRB0
        dmrs = nrPUSCHDMRSConfig;
        dmrs = passign(pusch(idx),dmrs,'PortSet','DMRSPortSet');
        dmrs = passign(pusch(idx),dmrs,'DMRSConfigurationType','DMRSConfigurationType',~pusch(idx).TransformPrecoding);
        dmrs = passign(pusch(idx),dmrs,'DMRSTypeAPosition');
        dmrs = passign(pusch(idx),dmrs,'DMRSAdditionalPosition');
        dmrs = passign(pusch(idx),dmrs,'DMRSLength');
        dmrs = passign(pusch(idx),dmrs,'DMRSSymbolSet','CustomSymbolSet');
        dmrs = passign(pusch(idx),dmrs,'NIDNSCID');
        dmrs = passign(pusch(idx),dmrs,'NSCID','NSCID',~pusch(idx).TransformPrecoding);
        dmrs = passign(pusch(idx),dmrs,'NRSID','NRSID',pusch(idx).TransformPrecoding);
        if isfield(pusch(idx),'GroupHopping') && pusch(idx).TransformPrecoding
            if strcmpi(pusch(idx).GroupHopping,'enable')
                dmrs.GroupHopping = 1;
            elseif strcmpi(pusch(idx).GroupHopping,'disable')
                dmrs.SequenceHopping = 1;
            end
        end
        if isfield(pusch(idx),'NumCDMGroupsWithoutData') && ~isempty(pusch(idx).NumCDMGroupsWithoutData) ...
                && ~pusch(idx).TransformPrecoding && pusch(idx).NumCDMGroupsWithoutData
            dmrs.NumCDMGroupsWithoutData = pusch(idx).NumCDMGroupsWithoutData;
        else
            dmrs.NumCDMGroupsWithoutData = 1 + pusch(idx).TransformPrecoding;
        end
        pusch2{idx}.DMRSPower             = pusch(idx).PowerDMRS;
        pusch2{idx}.DMRS = dmrs;

        numDMRSPorts = numel(dmrs.DMRSPortSet);
        if numDMRSPorts
            % Assign the number of layers to the number of DM-RS antenna ports
            pusch2{idx}.NumLayers = numDMRSPorts;
        end  
            
        % PT-RS configuration
        pusch2{idx}.EnablePTRS              = pusch(idx).EnablePTRS;
        ptrs{idx} = nrPUSCHPTRSConfig;
        ptrs{idx}.TimeDensity               = pusch(idx).PTRSTimeDensity;
        ptrs{idx}.FrequencyDensity          = pusch(idx).PTRSFrequencyDensity;
        ptrs{idx}.REOffset                  = pusch(idx).PTRSREOffset;
        ptrs{idx}.PTRSPortSet               = pusch(idx).PTRSPortSet;
        ptrs{idx}.NID                       = pusch(idx).PTRSNID;
        pusch2{idx}.PTRSPower               = pusch(idx).PowerPTRS;
        pusch2{idx}.PTRS = ptrs{idx};
    end
end

function csirs2 = mapCSIRS(cfgS)
    
    csirs2 = {};

    % CSIRS
    csirs = cfgS.CSIRS;
    for idx = 1:length(csirs)
        csirs2{idx} = nrWavegenCSIRSConfig;
        csirs2{idx}.Enable = csirs(idx).Enable;
        csirs2{idx}.Label = ['CSIRS' num2str(idx)];
        csirs2{idx}.Power = csirs(idx).Power;
        % Get the SCS of the BWP

        csirs2{idx}.BandwidthPartID = csirs(idx).BWP; 

        % This version of the mapper assumes a single CSI-RS per individual struct/object
        csirs2{idx}.CSIRSType = csirs(idx).CSIRSType;
        csirs2{idx}.RowNumber = csirs(idx).RowNumber;
        csirs2{idx}.Density = csirs(idx).Density;
        csirs2{idx}.SubcarrierLocations = csirs(idx).SubcarrierLocations;
        csirs2{idx}.NumRB = csirs(idx).NumRB;
        csirs2{idx}.RBOffset = csirs(idx).RBOffset;
        csirs2{idx}.SymbolLocations = csirs(idx).SymbolLocations;
        csirs2{idx}.NID = csirs(idx).NID;
        % Separate AllocatedSlots & AllocatedPeriod parameters need to be combined into a single CSIPeriod parameter
        % Mapper will create a single pattern
        allocperiod = csirs(idx).AllocatedPeriod;
        allocslots = csirs(idx).AllocatedSlots;
        if iscell(allocperiod)
            allocperiod = allocperiod{1};
        end
        if iscell(allocslots)
            allocslots = allocslots{1};
        end    
        if isempty(allocperiod)
            csirs2{idx}.CSIRSPeriod = 'off';
        elseif allocperiod == 1 
            csirs2{idx}.CSIRSPeriod = 'on';
        else
            csirs2{idx}.CSIRSPeriod = [allocperiod(1) allocslots(1)];
        end
    end
end

% Get the valid, usable PRB associated with a CORESET definition This
% function returns the PRB associated with the CORESET PRB which overlap
% the AllocatedPRB parameter. Accounts for a non-multiple-of-6 NStartBWP,
% by offsetting the CORESET to the next multiple-of-6 RB.
function [coresetprb,nblocks,nblockmax] = getCORESETPRB(carrier,bwp,coreset)

    nstartbwp = carrier.RBStart + bwp.RBOffset;     % Starting CRB of the BWP containing the CORESET
    offsetprb = mod(6-nstartbwp,6);                 % Offset from the BWP start to first complete block of 6 PRB

    % Turn each PRB into the associated block number that contains it

    % Block 0 should associate with BWP PRB of first complete block of 6, and
    % then expand the resulting set into the associated PRB
    nblocks = unique(floor(coreset.AllocatedPRB/6));                 % Indices of blocks associated by the CORESET AllocatedPRB parameter
    nblockmax = fix((bwp.NRB - offsetprb)/6);                        % Maximum number of complete block that can be contained in the BWP
    nblocks = nblocks( (nblocks >= 0) & (nblocks < nblockmax) );     % Complete block set contained within the BWP
    coresetprb = expander(6*nblocks,6) + offsetprb;                  % PRB associated with the resulting allocated CORESET
end


function [eallocation,eperiod] = controlexpansion(sallocation,speriod,callocation,cperiod)
    
    % Search space/CORESET sequence
    eperiod = [];
    if ~isempty(speriod)
        sallocation = sallocation(sallocation < speriod);
        if isempty(cperiod)
            srep = ceil(max(callocation)+1/length(sallocation));    % Repeat 's' so that it can be indexed by 'c'
        else
            srep = cperiod;                 % If cperiod is empty here then we may need to extend 's' sequence to cover  'c'
            eperiod = srep * speriod;
        end
    else
        srep = 1;
        speriod = 1;
    end
    
    % PDCCH sequence, wrt above search space sequence
    if ~isempty(cperiod)
        callocation = callocation(callocation < cperiod);
        crep = length(sallocation);
    else
        crep = 1;
        cperiod = 1;
    end
   
    % This will be indexed, so has to cover the indexing PDCCH sequence
    sexp = reshape(sallocation(:) + speriod*(0:srep-1),1,[]);

    % This will index the first search slot sequence 
    cexp = reshape(callocation(:) + cperiod*(0:crep-1),1,[]);
    
    % Select allocation from the search space slot sequence  
    cexp = cexp(cexp<length(sexp));
    eallocation = sexp(1+cexp);
    
end

% Expand 'd' by amount 'e', with optional non-unity strides and exclusion 
function expanded = expander(d,e,s,o,excl)
    if nargin < 5
        excl = 0;
    end
    if nargin < 4
        o = 0;
    end
    if nargin < 3
        s = 1;
    end
    eseq = (o:s:e-1)';
    if excl
        eseq = setdiff((0:e-1)',eseq);
    end
    expanded = reshape(reshape(d,1,[]) + eseq,1,[]);  % Use column expansion
end

