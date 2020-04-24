classdef MixedState
    properties
        States;          % [psi1 ... psin]
        NumberOfSystems; % [N1   ... Nn]
    end
    methods
        function obj = MixedState(states, numberOfSystems)
            obj.States = states;
            obj.NumberOfSystems = numberOfSystems;
        end
        function Pr = getProbability(MixedState, k)
            Pr = MixedState.NumberOfSystems(k) / sum(MixedState.NumberOfSystems);
        end
    end
end