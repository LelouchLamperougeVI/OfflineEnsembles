classdef ensemble < handle
    % Detection of offline ensembles in two-photon data.
    %
    % Companion code for:
    %
    %   Chang, H., Esteves, I. M., Neumann, A. R., Sun, J., Mohajerani, M.
    %   H., McNaughton, B. L. (2022), Cortical Reactivation of Non-Spatial
    %   and Spatial Memory Representations Coordinate with Hippocampus to
    %   Form a Memory Dialogue. bioRxiv.
    %   https://doi.org/10.1101/2022.12.16.520658
    
    properties (GetAccess = 'public', SetAccess = 'protected')
        twop = struct( 'fs',[], 'ts',[], 'deconv',[]);
        behaviour
        ensembles = struct('R', [], 'clust', [], 'order', [])
        hiepi
        ops
    end
    
    methods
        % instantiate object
        function obj = ensemble(behaviour, deconv, params, varargin)
            obj.behaviour = behaviour;
            obj.twop.deconv = deconv;
            obj.twop.ts = params.ts;
            obj.twop.fs = params.fs;
            obj.set(varargin);
        end
        
        % remove movement epochs
        function rm_mvt(obj)
            idx = running(obj.behaviour.unit_vel);
            idx = ~~conv(idx, ones(round(obj.twop.fs * 2), 1), 'same'); % dilate by 2 sec
            obj.twop.deconv(idx, :) = nan;
        end
        
        set(obj, varargin);
        corr(obj);
        hclust(obj);
        hiepi = hPICA(obj);
    end
end

