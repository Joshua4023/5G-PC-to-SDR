classdef(Abstract) hSDRBase < handle
%hSDRBase Base SDR helper class for RX/TX SDR object creation

%   Copyright 2022-2023 The MathWorks, Inc.

    properties(Dependent)
        % General Params
        CenterFrequency;
        DeviceAddress;
        ChannelMapping;
        SampleRate;
    end

    properties (Hidden, Dependent)
        MasterClockRate;
        SampleRateFactor;
    end

    properties (SetAccess=protected)
        DeviceName;
        SDRObj;
    end

    properties (Access = protected)
        ProvidedUSRPSampleRate;
    end

    properties(Constant, Hidden)
        ListOfUSRPs = {'N200/N210/USRP2', 'N300', 'N310', 'N320/N321', ...
            'B200', 'B210', 'X300', 'X310'};
        ListOfOtherTransceivers = {'AD936x', 'FMCOMMS5', 'E3xx', 'Pluto'};
    end

    methods

        %MasterClockRate
        function set.MasterClockRate(obj,value)
            if matches(obj.DeviceName, obj.ListOfUSRPs)
                obj.SDRObj.MasterClockRate = value;
            end
        end

        function value = get.MasterClockRate(obj)
            if matches(obj.DeviceName, obj.ListOfUSRPs)
                value = obj.SDRObj.MasterClockRate;
            end
        end



        %CenterFrequency
        function set.CenterFrequency(obj,value)
            obj.SDRObj.CenterFrequency = value;
        end

        function value = get.CenterFrequency(obj)
            value = obj.SDRObj.CenterFrequency;
        end

        %DeviceAddress
        function set.DeviceAddress(obj,value)
            switch obj.DeviceName
                case {'AD936x', 'FMCOMMS5', 'E3xx', 'N200/N210/USRP2', ...
                        'N300', 'N310', 'N320/N321', 'X300', 'X310'}
                    obj.SDRObj.IPAddress = value;
                case 'Pluto'
                    obj.SDRObj.RadioID = value;
                case 'RTL-SDR'
                    obj.SDRObj.RadioAddress = value;
                case {'B200', 'B210'}
                    obj.SDRObj.SerialNum = value;
            end
        end

        function value = get.DeviceAddress(obj)
            switch obj.DeviceName
                case {'AD936x', 'FMCOMMS5', 'E3xx', 'N200/N210/USRP2', ...
                        'N300', 'N310', 'N320/N321', 'X300', 'X310'}
                    value = obj.SDRObj.IPAddress;
                case {'B200', 'B210'}
                    value = obj.SDRObj.SerialNum;
                case {'Pluto'}
                    value = obj.SDRObj.RadioID;
                case 'RTL-SDR'
                    value = obj.SDRObj.RadioAddress;
            end
        end

        %ChannelMapping
        function set.ChannelMapping(obj,value)
            if matches(obj.DeviceName, [obj.ListOfOtherTransceivers obj.ListOfUSRPs])
                obj.SDRObj.ChannelMapping = value;
            else
                error("hSDRBase:setChannelMapping","ChannelMapping Property is read-only for %s", obj.DeviceName);
            end
        end

        function value = get.ChannelMapping(obj)
            if matches(obj.DeviceName, [obj.ListOfOtherTransceivers obj.ListOfUSRPs])
                value = obj.SDRObj.ChannelMapping;
            else
                value = 1;
            end
        end

        %SampleRateFactor
        function set.SampleRateFactor(obj, value)
            if matches(obj.DeviceName, obj.ListOfUSRPs)
                if isa(obj.SDRObj, 'comm.SDRuReceiver')
                    obj.SDRObj.DecimationFactor = value;
                else
                    obj.SDRObj.InterpolationFactor = value;
                end
            end
        end

        function value = get.SampleRateFactor(obj)
            if matches(obj.DeviceName, obj.ListOfUSRPs)
                if isa(obj.SDRObj, 'comm.SDRuReceiver')
                    value = obj.SDRObj.DecimationFactor;
                else
                    value = obj.SDRObj.InterpolationFactor;
                end
            end
        end

        %SampleRate
        function set.SampleRate(obj, value)
            switch obj.DeviceName
                case obj.ListOfOtherTransceivers
                    obj.SDRObj.BasebandSampleRate = value;
                case 'RTL-SDR'
                    obj.SDRObj.SampleRate = value;
                case obj.ListOfUSRPs
                    [obj.MasterClockRate, obj.SampleRateFactor] ...
                        = obj.hValidateUSRPSampleRate(obj.DeviceName, value);
            end
        end

        function value = get.SampleRate(obj)
            switch obj.DeviceName
                case obj.ListOfOtherTransceivers
                    value = obj.SDRObj.BasebandSampleRate;
                case 'RTL-SDR'
                    value = obj.SDRObj.SampleRate;
                case obj.ListOfUSRPs
                    value = obj.MasterClockRate/obj.SampleRateFactor;
            end
        end

        function release(obj)
            release(obj.SDRObj);
        end

        function out = info(obj)
            out = info(obj.SDRObj);
        end

    end


    methods(Access=protected)

        function [mcr, f] = hValidateUSRPSampleRate(~, platform, sampleRate)
            % HGETUSRPRATEINFORMATION function provides the master clock rate and the
            % interpolation/decimation factor given a USRP platform and a desired
            % sampleRate. If the sample rate is not realizable using the provided
            % platform then an error is thrown informing the user of this. See
            % comm.SDRuTranmitter or comm.SDRuReceiver documentation pages for further
            % information on supported master clock rates and interpolation/decimation
            % factors.
            switch platform
                case 'N200/N210/USRP2'
                    masterClockRate = 100e6;
                    factor = [4:128 130:2:256 260:4:512];

                case {'N300', 'N310'}
                    masterClockRate = [122.88e6 125e6 153.6e6];
                    factor = [1:4 6:2:128 130:2:256 260:4:512 520:8:1024];

                case 'N320/N321'
                    masterClockRate = [200e6 245.76e6 250e6];
                    factor = [1:4 6:2:128 130:2:256 260:4:512 520:8:1024];

                case {'B200', 'B210'}
                    minMasterClockRate = 5e6;
                    maxMasterClockRate = 61.44e6;
                    masterClockRate = minMasterClockRate:1e3:maxMasterClockRate;
                    factor = [1:128 130:2:256 260:4:512];

                case {'X300', 'X310'}
                    masterClockRate = [184.32e6 200e6];
                    factor = [1:128 130:2:256 260:4:512];

                otherwise
                    masterClockRate = nan;
                    factor = nan;
            end

            possibleSampleRates = masterClockRate'./factor;
            % do not consider smaller sample rates, to satisfy Nyquist:
            possibleSampleRates(possibleSampleRates<sampleRate) = NaN;

            err = abs(possibleSampleRates - sampleRate);
            minErr = min(err,[],"all");
            if isnan(minErr)
                error("hSDRBase:hValidateUSRPSampleRate","The sample rate %.2g is not realizable using the %s radio.",sampleRate,platform);
            end

            [idx1, idx2] = find(err==minErr);
            mcr = masterClockRate(idx1(1));
            f = factor(idx2(1));
        end

        function findUSRP(~, deviceName,sdrObj)
            foundUSRPs = findsdru;
            deviceStatus = foundUSRPs({foundUSRPs.Platform} == deviceName);
            if ~isempty(deviceStatus)
                if matches(deviceName, {'B200', 'B210'})
                    sdrObj.SerialNum = deviceStatus(1).SerialNum;
                else
                    sdrObj.IPAddress = deviceStatus(1).IPAddress;
                end
            else
                warning("hSDRBase:deviceNotFound","Device %s not found", deviceName);
            end
        end

    end
end

