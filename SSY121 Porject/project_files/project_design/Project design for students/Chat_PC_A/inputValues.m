classdef inputValues
    properties
        fs
        fc
        fsymb
        pulseShape
        constellation
        preamble_carrier
        
    end
    methods
        function inputval = inputValues(fs, fc,fsymb,pulseShape,constellation,preamble_carrier)
            inputval.fs = fs;
            inputval.fc = fc;
            inputval.fsymb = fsymb;
            inputval.pulseShape = pulseShape;
            inputval.constellation = constellation;
            inputval.preamble_carrier = preamble_carrier;
            
        end
    end
end