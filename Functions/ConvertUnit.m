function out_val = ConvertUnit(val, unit_to, varargin)
%CONVERTUNIT converts value val (with units) to units of unit_to
    
    p = inputParser;
    addParameter(p, 'Temperature', 'none')
    parse(p, varargin{:});
    
    if strcmp(p.Results.Temperature, 'none')
        out_val = double(separateUnits(unitConvert(val, unit_to)));
    else
        out_val = double(separateUnits(unitConvert(val, unit_to,...
            'Temperature', p.Results.Temperature)));
    end
end