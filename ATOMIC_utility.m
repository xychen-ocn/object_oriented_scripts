classdef ATOMIC_utility
    properties
        varname
        svdatadir
        svdataFN
    end
    
    methods
        
        function obj = ATOMIC_utility(val)
            if nargin> 0
                PropNames = fieldnames(val);
                for i = 1:length(PropNames)
                    PN = PropNames{i};
                    obj.(PN) = val.(PN);
                end
                
            else
                obj.svdatadir = 'none';
                obj.svdataFN = 'none';
                obj.value = [];
                
            end
        end
        
        
        function savedata(obj)
            absFN = [obj.svdatadir filesep obj.svdataFN];
            varname = obj.varname;
            feval(save, absFN, varname);
        end 
        
    end
end