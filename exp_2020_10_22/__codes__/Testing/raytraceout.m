function [dBmOut] = raytraceout(txs, rxs, varargin)
%raytrace   Plot propagation paths between sites
%   raytrace(TX,RX) plots propagation paths from transmitter site TX to
%   receiver site RX. The propagation paths are found using ray tracing
%   with the terrain and buildings data defined in the Site Viewer map.
%   Each propagation path is colored according to the received power (dBm)
%   or path loss (dB) along the path, assuming unpolarized rays.
%
%   raytrace(___,Name,Value) plots propagation paths with additional
%   options specified by one or more Name-Value pairs.
%
%   The inputs TX and RX can be scalars or arrays. If specified as arrays,
%   then propagation paths are plotted from each transmitter site in TX to
%   each receiver site in RX.
%
%   raytrace Name-Value pairs:
%
%   Type - Type of quantity to plot, specified as 'power' or 'pathloss'.
%      When Type is 'power', each path is colored according to the received
%      power (dBm) along the path. When Type is 'pathloss', each path is
%      colored according to the path loss (dB) along the path. The default
%      value is 'power'.
%
%   NumReflections - Number of reflections in propagation paths to search
%      for using ray tracing, specified as a numeric row vector whose
%      elements are 0, 1, or 2. The default value is [0 1], which results
%      in a search for a line-of-sight propagation path along with
%      propagation paths that each contain a single reflection.
%
%   ColorLimits - Color limits for colormap, specified as a two-element
%      numeric vector of the form [min max], expressed in dBm when Type is
%      'power' and dB when Type is 'pathloss'. The color limits indicate
%      the values that map to the first and last colors in the colormap.
%      Propagation paths with values below the minimum color limit are not
%      plotted. The default value is [-120 -5] if Type is 'power' and [45
%      160] if Type is 'pathloss'.
%
%   Colormap - Colormap for coloring propagation paths, specified as a 
%      predefined colormap name or an M-by-3 array of RGB (red, green,
%      blue) triplets that define M individual colors. The default value is
%      'jet'.
%
%   ShowLegend - Show color legend on map, specified as true or false. The 
%      default value is true.
%
%   Map - Map for visualization and surface data, specified as a siteviewer
%      object. The default value is the current siteviewer or else a new
%      siteviewer if none is open.
%
%   Notes
%   -----
%   - The ray tracing analysis includes surface reflections but does not
%     include effects from refraction, diffraction, or scattering.
%
%   - Path loss and received power values do not include reflection loss
%     due to material or antenna polarization effects.
%
%   % Example: Show reflected propagation paths
%
%   % Launch Site Viewer with buildings in Chicago
%   viewer = siteviewer("Buildings","chicago.osm");
%
%   % Create transmitter site on a building
%   tx = txsite('Latitude',41.8800, ...
%      'Longitude',-87.6295, ...
%      'TransmitterFrequency',2.5e9);
%   show(tx)
%
%   % Create receiver site near another building
%   rx = rxsite('Latitude',41.881352, ...
%      'Longitude',-87.629771, ...
%      'AntennaHeight',30);
%   show(rx)
%
%   % Show obstruction to line-of-sight
%   los(tx,rx)
%
%   % Show reflected propagation path using ray tracing
%   raytrace(tx,rx)
%
%   % Show propagation paths including first order and second order
%   % reflections
%   raytrace(tx,rx,'NumReflections',[1 2])
%
%   See also los, siteviewer

%   Copyright 2019 The MathWorks, Inc.

% Validate sites
validateattributes(txs,{'txsite'},{'nonempty'},'raytrace','',1);
validateattributes(rxs,{'rxsite'},{'nonempty'},'raytrace','',2);

% Process optional name/value pairs
p = inputParser;
p.addParameter('Animation', '');
p.addParameter('EnableWindowLaunch', true);
p.addParameter('Type', 'power');
p.addParameter('NumReflections', [0 1]);
p.addParameter('ColorLimits', []);
p.addParameter('Colormap', 'jet');
p.addParameter('ShowLegend', true);
p.addParameter('Map', []);
p.parse(varargin{:});

% Get Site Viewer and validate web graphics
viewer = rfprop.internal.Validators.validateMap(p, 'raytrace');
isViewerInitiallyVisible = viewer.Visible;
viewer.validateWebGraphicsSupport;

% Validate and get parameters
[animation, enableWindowLaunch] = rfprop.internal.Validators.validateGraphicsControls(p, 'raytrace');
type = validateType(p);
isPathlossPlot = strcmp(type,'pathloss');
numReflections = validateNumReflections(p);
cmap = rfprop.internal.Validators.validateColorMap(p, 'raytrace');
if isPathlossPlot
    defaultColorLimits = [45 160];    
    legendTitle = message('shared_channel:rfprop:RaytracePathlossLegendTitle').getString;
    
    % Flip colormap to maintain consistency of color meaning (red end of
    % jet means good signal), since path loss is inversely proportional to
    % received power
    cmap = flipud(cmap);
else
    defaultColorLimits = [-120 -5];
    legendTitle = message('shared_channel:rfprop:RaytracePowerLegendTitle').getString;
end
clim = rfprop.internal.Validators.validateColorLimits(p, defaultColorLimits, 'raytrace');
showLegend = rfprop.internal.Validators.validateShowLegend(p, 'raytrace');

% Get site antenna coordinates
txsCoords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(txs, viewer);
rxsCoords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(rxs, viewer);

noPathsFound = true;
noPathsMeetThreshold = true;
dBmOut = [];
% Generate triangulation of 3D environment for tx/rx locations
txsLatLon = txsCoords.LatitudeLongitude;
rxsLatLon = rxsCoords.LatitudeLongitude;
lats = [txsLatLon(:,1); rxsLatLon(:,1)];
lons = [txsLatLon(:,2); rxsLatLon(:,2)];
if all(txsCoords.withinBuildingsLimits(lats, lons))
    % Use previously calculated triangulation for buildings region
    tri = viewer.BuildingsTerrainTriangulation;
else
    % Get triangulation for region, which includes buildings if they exist
    [latmin, latmax] = bounds(lats(:));
    [lonmin, lonmax] = bounds(lons(:));
    tri = viewer.regionTriangulation([latmin latmax], [lonmin lonmax]);
end

% Compute rays from each transmitter to all receivers
startIDs = {};
endIDs = {};
rayGroups = {};
for txInd = 1:numel(txs)
    tx = txs(txInd);
    fq = tx.TransmitterFrequency;
    txAxes = antennaAxes(tx.AntennaAngle);
    txCoords = txsCoords.extract(txInd);
    txLatLon = txCoords.LatitudeLongitude;
    txElevation = txCoords.AntennaElevation;
    txPos = txCoords.enuFromBuildingsCenter;
    
    for rxInd = 1:numel(rxs)
        rx = rxs(rxInd);
        rxAxes = antennaAxes(rx.AntennaAngle);
        rxCoords = rxsCoords.extract(rxInd);
        rxLatLon = rxCoords.LatitudeLongitude;
        rxElevation = rxCoords.AntennaElevation;
        rxPos = rxCoords.enuFromBuildingsCenter;
        
        % Perform ray tracing analysis
        rays = comm.internal.raytrace(tri, txPos', rxPos', numReflections, fq, txAxes, rxAxes);
        
        % Initial data variables to descript rays in this group
        raysLatitudes = {};
        raysLongitudes = {};
        raysElevations = {};
        raysColors = {};
        raysInfo = {};
        
        % Analyze each ray
        numRays = numel(rays);
        for rayInd = 1:numRays
            noPathsFound = false;
            ray = rays(rayInd);
            
            % Ray path starts with tx position
            rayLats = txLatLon(1);
            rayLons = txLatLon(2);
            rayElevs = txElevation;
            
            % Add reflection positions to ray path
            if ~ray.LineOfSight
                rayInterfacePositions = ray.InterfacePositions;
                for reflectionInd = 1:size(rayInterfacePositions,2)
                    rayReflectionPos = rayInterfacePositions(:,reflectionInd);
                    [rayLats(end+1),rayLons(end+1),rayElevs(end+1)] = ...
                        txCoords.geodeticFromBuildingsCenter(rayReflectionPos');
                end
            end
            
            % Ray path ends with rx position
            rayLats(end+1) = rxLatLon(1); %#ok<*AGROW>
            rayLons(end+1) = rxLatLon(2);
            rayElevs(end+1) = rxElevation;
            
            % Compute AoD/AoA angles
            [aodaz, aodel] = localangle(rayLats(1), rayLons(1), rayElevs(1), ...
                rayLats(2), rayLons(2), rayElevs(2));
            [aoaaz, aoael] = localangle(rayLats(end), rayLons(end), rayElevs(end), ...
                rayLats(end-1), rayLons(end-1), rayElevs(end-1));
            
            % Compute data value for ray
            pl = ray.Pathloss;
            if isPathlossPlot
                v = pl;
            else
                % Compute signal strength using Friis equation
                Ptx_db = 10 * log10(1000*tx.TransmitterPower); % Convert W to dBm
                Gtx_db = gain(tx, fq, aodaz, aodel); % Transmit gain
                Grx_db = gain(rx, fq, aoaaz, aoael); % Receive gain
                v = Ptx_db + Gtx_db + Grx_db - pl - tx.SystemLoss - rx.SystemLoss;
                dBmOut = [dBmOut v];
            end
            
            % Discard ray if it does not meet threshold
            if isPathlossPlot
                % Check if path loss is greater than max value
                rayDoesNotMeetThreshold = v > clim(2);
            else
                % Check if power is less than min value
                rayDoesNotMeetThreshold = v < clim(1);
            end
            if rayDoesNotMeetThreshold
                continue
            else
                noPathsMeetThreshold = false;
            end
            
            % Compute color value for ray
            rayColorRGB = rfprop.internal.ColorUtils.colorcode(v, cmap, clim);
            raysColors{end+1} = rfprop.internal.ColorUtils.rgb2css(rayColorRGB);
            raysLatitudes{end+1} = rayLats;
            raysLongitudes{end+1} = rayLons;
            raysElevations{end+1} = rayElevs;

            % Build infobox description
            powDecPlaces = rfprop.Constants.MaxPowerInfoDecimalPlaces;
            distDecPlaces = rfprop.Constants.MaxDistanceInfoDecimalPlaces;
            phaseDecPlaces = rfprop.Constants.MaxPhaseInfoDecimalPlaces;
            angDecPlaces = rfprop.Constants.MaxAngleInfoDecimalPlaces;
            if isPathlossPlot
                valueDescription = message('shared_channel:rfprop:RaytraceDescriptionPathloss',  mat2str(round(v, powDecPlaces))).getString;
            else
                valueDescription = message('shared_channel:rfprop:RaytraceDescriptionPower', mat2str(round(v, powDecPlaces))).getString;
            end
            dist = ray.TimeOfArrival*physconst('lightspeed');
            desc = [...
                message('shared_channel:rfprop:RaytraceDescriptionNumReflections', mat2str(ray.NumBounces)).getString, '<br>', ...
                valueDescription, '<br>', ...
                message('shared_channel:rfprop:RaytraceDescriptionPhaseChange', mat2str(round(ray.PhaseShift, phaseDecPlaces))).getString, '<br>', ...
                message('shared_channel:rfprop:RaytraceDescriptionDistance', mat2str(round(dist, distDecPlaces))).getString, '<br>', ...
                message('shared_channel:rfprop:RaytraceDescriptionAoD', mat2str(round(aodaz, angDecPlaces)), mat2str(round(aodel, angDecPlaces))).getString, '<br>', ...
                message('shared_channel:rfprop:RaytraceDescriptionAoA', mat2str(round(aoaaz, angDecPlaces)), mat2str(round(aoael, angDecPlaces))).getString];
            raysInfo{end+1} = desc;
        end
        
        % Add ray group data
        %startIDs{end+1} = tx.UID;
        %endIDs{end+1} = rx.UID;
        rayGroups{end+1} = struct(...
            'RaysLatitudes', {raysLatitudes}, ...
            'RaysLongitudes', {raysLongitudes}, ...
            'RaysElevations', {raysElevations}, ...
            'RaysColors', {raysColors}, ...
            'RaysInfo', {raysInfo});
    end
end

% Show sites
if isViewerInitiallyVisible && viewer.Visible
    return % Abort if Site Viewer has been closed
end
%show(txs,'AntennaSiteCoordinates',txsCoords,'Map',viewer, ...
%    'ShowAntennaHeight',true,'Animation','none','EnableWindowLaunch',false);
%show(rxs,'AntennaSiteCoordinates',rxsCoords,'Map',viewer, ...
%    'ShowAntennaHeight',true,'Animation','none','EnableWindowLaunch',false);

% Warn if no paths to plot and return early
if noPathsFound
    warning(message('shared_channel:rfprop:NoPathsMeetRequiredNumReflections'));
    return;
elseif noPathsMeetThreshold
    if isPathlossPlot
        warning(message('shared_channel:rfprop:NoPathsMeetRequiredPathLoss', mat2str(clim(2))));
    else
        warning(message('shared_channel:rfprop:NoPathsMeetRequiredPower', mat2str(clim(1))));
    end
    return;
end

% Get legend data
if showLegend
    [legendColors, legendColorValues] = rfprop.internal.ColorUtils.colormaplegend(cmap, clim);
else
    legendColors = "";
    legendColorValues = "";
end
if isPathlossPlot
    legendColors = fliplr(legendColors);
    legendColorValues = fliplr(legendColorValues);
end

% Plot rays
data = struct(...
    'StartIDs', {startIDs}, ...
    'EndIDs', {endIDs}, ...
    'RayGroups', {rayGroups}, ...
    'Animation', animation, ...
    'EnableWindowLaunch', enableWindowLaunch, ...
    'ShowLegend', showLegend, ...
    'LegendTitle', legendTitle, ...
    'LegendColors', legendColors, ...
    'LegendColorValues', legendColorValues);
if isViewerInitiallyVisible && ~viewer.Visible
    return % Abort if Site Viewer has been closed
end
viewer.rays(data);
close(viewer)
end

function type = validateType(p)

try
    type = p.Results.Type;
    type = validatestring(type, {'pathloss','power'}, 'raytrace', 'Type');
catch e
    throwAsCaller(e);
end
end

function numReflections = validateNumReflections(p)

try
    numReflections = p.Results.NumReflections;
    validateattributes(numReflections, {'numeric'}, ...
        {'real','finite','nonsparse','row','integer','nonnegative','nonempty','<=',2}, ...
        'raytrace', 'NumReflections');
    numReflections = unique(numReflections);
catch e
    throwAsCaller(e);
end
end

function [az,el] = localangle(lat0, lon0, ht0, lat1, lon1, ht1)
%localangle   Return angles between locations using local ENU coordinate system

% Get ENU position of second location relative to first
[X,Y,Z] = rfprop.internal.MapUtils.geodetic2enu(...
    lat0, lon0, ht0, lat1, lon1, ht1);

% Convert Cartesian position to az/el angles
[azrad,elrad] = cart2sph(X,Y,Z);
az = rad2deg(azrad);
el = rad2deg(elrad);

% Guarantee azimuth in range [-180,180]
az = wrapTo180(az);
end

function axes = antennaAxes(antAngle)
%antennaAxes   Return antenna axes

phi = antAngle(1);
if numel(antAngle) > 1
    theta = antAngle(2);
else
    theta = 0;
end
axes = ([1 0 0; 0 cosd(theta) sind(theta); 0 -sind(theta) cosd(theta)] * ...
    [cosd(phi) sind(phi) 0; -sind(phi) cosd(phi) 0; 0 0 1])';
end


