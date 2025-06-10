%-- X class for GESLA-X data

% Initialise

% gX_obj = sls.gesla3.X(metadata="meta path",directory="data path");

% gX_obj contains the following functions:

% site2file: Get one or more GESLA-X file name/s from a site name or part of a site name
% load: Load all data, a group, or a dataset from a GESLA-X file
% load_multiple: Load data from multiple GESLA-X files
% load_bbox: GESLA-X data from within lon/lat bounding box, bbox -> [northern extent, southern extent, western extent, eastern extent]
% load_nearest: N-number of GESLA-X files with nearest lon/lat to given coordinates in *Euclidean distance*, ties can overwrite 'nsites' input
% load_country: GESLA-X files from specific country, country -> 3 letter country code used in GESLA, e.g. 'GBR' or 'JPN'
% change_fieldnames: Change the field names in an output data struct to different GESLA-X format

% e.g.,

% gX_obj.site2file("sitename","site name")

% gX_obj.load(file="your file",group="timeseries",dataset="water_levels")

% gX_obj.load_multiple(files=["file name 1","file name 2"],group="tides",dataset="tides_15min_res")

% gX_obj.load_bbox(bbox=[northern extent, southern extent, western extent, eastern extent],group="skew_surge",dataset="time_of_actual_high_water")

% gX_obj.load_nearest(lonlat=[longitude, latitude],nsites=5,group="annual_maxima",dataset="non_tidal_residual")

% gX_obj.load_country(country="GBR",group="return_levels",dataset="gev_CI")

% these functions can accept calls with or without group= and dataset= , etc.
% e.g.,
% gX_obj.load(file="your file") %-- loads all data from file
% gX_obj.load(file="your file",group="timeseries") %-- loads all timeseries data from file
% gX_obj.load(file="your file",group="timeseries",dataset="water_levels") %-- loads the water levels from the timeseries group in file

% Luke Jenkins February 2024

classdef gX < handle & matlab.mixin.SetGet
    properties
        metadata
        data_directory string
    end
    methods
        function obj = gX(props)
            arguments
                props.?gX
            end
            set(obj,props)
            obj.metadata = readtable(obj.metadata);
        end
        function filename = site2file(obj,nv)
            % Get one or more GESLA-X file name/s from a site name or part of a site name
            arguments
                obj
                nv.sitename string
            end
            filename = string(obj.metadata.FILENAME(contains(obj.metadata.SITENAME,nv.sitename,'IgnoreCase',true)));
            if isempty(filename) % no sites found
                error("No file name found for site name: " + nv.sitename)
            elseif ~isStringScalar(filename)
                prompt = {strcat("-Site: ",string(obj.metadata.SITENAME(contains(obj.metadata.FILENAME,filename,'IgnoreCase',true))), " -Country: ",...
                    string(obj.metadata.COUNTRY(contains(obj.metadata.FILENAME,filename,'IgnoreCase',true))), " -File: ",filename)};
                dlgtitle = "Multiple files found: Select which file names you wish to output (Place '1' in dialogue box)";
                answer = inputdlg(prompt{:},dlgtitle,[1 length(char(dlgtitle))+15]);
                filename = filename(~cellfun(@isempty, answer));
            end
        end
        function data = load(obj,nv)
            % Load all data, a group, or a dataset from a GESLA-X file
            arguments
                obj
                nv.file string
                nv.group string = strings(0)
                nv.dataset string = strings(0)
            end
            data = gh5load(obj.data_directory + nv.file,nv.group,nv.dataset);
            if isempty(nv.group) && isempty(nv.dataset); data.sitename = string(obj.metadata.SITENAME(ismember(obj.metadata.FILENAME,nv.file))); end
        end
        function data = load_multiple(obj,nv)
            % Load data from multiple GESLA-X files
            arguments
                obj
                nv.files string
                nv.group string = strings(0)
                nv.dataset string = strings(0)
            end
            data = gh5load_multiple(obj,nv.files,nv.group,nv.dataset);
        end
        function data = load_bbox(obj,nv)
            % GESLA-X data from within lon/lat bounding box, bbox -> [northern extent, southern extent, western extent, eastern extent]
            arguments
                obj
                nv.bbox (1,4) double 
                nv.group string = strings(0)
                nv.dataset string = strings(0)
            end
            filenames = string(obj.metadata.FILENAME((obj.metadata.LATITUDE >= nv.bbox(2)  &  obj.metadata.LATITUDE <= nv.bbox(1)) &...
                (obj.metadata.LONGITUDE >= nv.bbox(3)  &  obj.metadata.LONGITUDE <= nv.bbox(4))));
            data = gh5load_multiple(obj,filenames,nv.group,nv.dataset);
        end
        function data = load_nearest(obj,nv)
            % N-number of GESLA-X files with nearest lon/lat to given coordinates in *Euclidean distance*, ties can overwrite 'nsites' input
            arguments
                obj
                nv.lonlat (1,2) double 
                nv.nsites (1,1) double = 1
                nv.group string = strings(0)
                nv.dataset string = strings(0)
            end
            coords = [obj.metadata.LONGITUDE obj.metadata.LATITUDE];
            i = cell2mat(knnsearch(coords,nv.lonlat,'K',nv.nsites,'IncludeTies',true))';
            filenames = string(obj.metadata.FILENAME(i));
            if length(filenames) == 1
                data = gh5load(obj.data_directory + filenames,nv.group,nv.dataset);
                if isempty(nv.group) && isempty(nv.dataset); data.sitename = string(obj.metadata.SITENAME(ismember(obj.metadata.FILENAME,filenames))); end
            else
                data = gh5load_multiple(obj,filenames,nv.group,nv.dataset);
            end
        end
        function data = load_country(obj,nv)
            % GESLA-X files from specific country, country -> 3 letter country code used in GESLA, e.g. 'GBR' or 'JPN'
            arguments
                obj
                nv.country string
                nv.group string = strings(0)
                nv.dataset string = strings(0)
            end
            filenames = string(obj.metadata.FILENAME(contains(obj.metadata.COUNTRY,string(nv.country),'IgnoreCase',true)));
            data = gh5load_multiple(obj,filenames,nv.group,nv.dataset);
        end
        function data = change_fieldnames(obj,nv)
            % Change the field names in an output data struct to different GESLA-X format
            arguments
                obj
                nv.struct struct
            end
            dfields = string(strrep(fieldnames(nv.struct),'_',''));
            files = string(replace(obj.metadata.FILENAME,{'_','-'},{''}));
            options = {string(obj.metadata.SITENAME(contains(files, dfields)));...
                string(obj.metadata.COUNTRY(contains(files, dfields)));...
                string(obj.metadata.CONTRIBUTOR_ABBREVIATED_(contains(files, dfields)));...
                string(obj.metadata.SITECODE(contains(files, dfields)))};
            prompt = {"site name","country","contributor (abbreviated)","site code (cannot be used on own)"};
            dlgtitle = "Select new field names (Place '1' in dialogue box)";
            answer = inputdlg(prompt,dlgtitle,[1 length(char(dlgtitle))+35]);
            options = options(~cellfun(@isempty, answer));
            new_names = strings();
            for i = 1:length(options)
                switch i
                    case 1
                        new_names = append(new_names,string(options{i}));
                    otherwise
                        new_names = append(new_names,"_",string(options{i}));
                end
            end
            data = cell2struct(struct2cell(nv.struct), new_names);
        end
    end
end

function data = gh5load_multiple(obj,filenames,group,dataset)
    for fname = filenames'
        data.(erase(strrep(fname,'-','_'),'.h5')) = gh5load(obj.data_directory + fname,group,dataset);
        if isempty(group) && isempty(dataset); data.(strrep(fname,'-','_')).sitename = string(obj.metadata.SITENAME(ismember(obj.metadata.FILENAME,fname))); end
    end
end

function data = gh5load(file,group,dataset)
    arguments
        file string
        group string = strings(0)
        dataset string = strings(0)
    end
    
    if ~endsWith(file,".h5"); file = file + ".h5"; end

    if isempty(group) || strcmpi(group, "gauge_info")
        x = "/tide_gauge_info/";
        data.datum = h5read(file, x + "datum");
        data.gauge_type = h5read(file, x + "gauge_type");
        data.contributor_flag_information = h5read(file, x + "contributor_flag_information");
        cf = h5read(file, x + "contributor_flags");
        gf = h5read(file, x + "gesla_flags");
        data.flags = table(cf, gf, 'VariableNames', ["contributor_flags", "gesla_flags"]);
        data.latitude = h5read(file, x + "latitude");
        data.longitude = h5read(file, x + "longitude");
    end
    
    if isempty(group) || strcmpi(group, "timeseries")
        % time
        x = "/times/";
        time = datetime(h5read(file, x + "time"), 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
        % water levels
        x = "/water_levels/";
        wl = h5read(file, x + "water_level");
        wld = h5read(file, x + "water_level_detrended");
        % tides
        x = "/tides/";
        tide = h5read(file, x + "tide");
        % surge
        x = "/surge/";
        ntr = h5read(file, x + "non_tidal_residual");
        if ~isempty(group)
            % timeseries
            data = timetable(wl, wld, tide, ntr, 'RowTimes', time, ...
                'VariableNames', ["water_level", "water_level_detrended", ...
                "tide", "non_tidal_residual"]);
        else
            data.timeseries = timetable(wl, wld, tide, ntr, 'RowTimes', time, ...
                'VariableNames', ["water_level", "water_level_detrended", ...
                "tide", "non_tidal_residual"]);
        end
    end
    
    if isempty(group) || strcmpi(group, "mean_sea_level")
        x = "/water_levels/";
        data.mean_sea_level_trend = h5read(file, x + "mean_sea_level_trend");
        data.mean_sea_level = h5read(file, x + "mean_sea_level");
    end
    
    if isempty(group) || strcmpi(group, "tides")
        x = "/tides/";
        data.harmonic_analysis = harmonic2table(h5read(file, x + "harmonic_analysis"));
        t15 = h5read(file, x + "tide_15min_res");
        t15ts = datetime(h5read(file, x + "tide_15min_res_times"), 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
        data.tide_15min_res = timetable(t15, 'RowTimes', t15ts, 'VariableNames', "tide");
    end
    
    if isempty(group) || strcmpi(group, "skew_surge")
        x = "/skew_surge/";
        sk = h5read(file, x + "skew_surge");
        ttp = h5read(file, x + "time_to_predicted_from_actual_high_water_(minutes)");
        pt = h5read(file, x + "predicted_tidal_high_water");
        top = datetime(h5read(file, x + "time_of_predicted_tidal_high_water"), 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
        aw = h5read(file, x + "actual_high_water");
        taw = datetime(h5read(file, x + "time_of_actual_high_water"), 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
        if ~isempty(group)
            data = table(sk, ttp, pt, top, aw, taw, 'VariableNames', ...
                ["skew_surge", "time_to_predicted_from_actual_high_water_(minutes)", ...
                "predicted_tidal_high_water", "time_of_predicted_tidal_high_water", ...
                "actual_high_water", "time_of_actual_high_water"]);
        else
            data.skew_surge = table(sk, ttp, pt, top, aw, taw, 'VariableNames', ...
                ["skew_surge", "time_to_predicted_from_actual_high_water_(minutes)", ...
                "predicted_tidal_high_water", "time_of_predicted_tidal_high_water", ...
                "actual_high_water", "time_of_actual_high_water"]);
        end
    end
    
    if isempty(group) || strcmpi(group, "annual_maxima")
        x = "/annual_maxima/";
        yr = h5read(file, x + "year");
        wl = h5read(file, x + "water_level");
        twl = datetime(h5read(file, x + "time_of_water_level"), 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
        wld = h5read(file, x + "water_level_detrended");
        twld = datetime(h5read(file, x + "time_of_water_level_detrended"), 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
        ntr = h5read(file, x + "non_tidal_residual");
        tntr = datetime(h5read(file, x + "time_of_non_tidal_residual"), 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
        sk = h5read(file, x + "skew_surge");
        tsk = datetime(h5read(file, x + "time_of_skew_surge"), 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
        if ~isempty(group)
            data = table(yr, wl, twl, wld, twld, ntr, tntr, sk, tsk, 'VariableNames', ...
                ["year", "water_level", "time_of_water_level", "water_level_detrended", "time_of_water_level_detrended", ...
                "non_tidal_residual", "time_of_non_tidal_residual", "skew_surge", "time_of_skew_surge"]);
        else
            data.annual_maxima = table(yr, wl, twl, wld, twld, ntr, tntr, sk, tsk, 'VariableNames', ...
                ["year", "water_level", "time_of_water_level", "water_level_detrended", "time_of_water_level_detrended", ...
                "non_tidal_residual", "time_of_non_tidal_residual", "skew_surge", "time_of_skew_surge"]);
        end
    end
    
    if isempty(group) || strcmpi(group, "return_levels")
        x = "/return_levels/";
        return_periods = h5read(file,x + "return_periods");
        gev = h5read(file,x + "gev");
        gum = h5read(file,x + "gum");
        gp = h5read(file,x + "gp");
        gev_CI = h5read(file,x + "gev_CI");
        gum_CI = h5read(file,x + "gum_CI");
        gp_CI = h5read(file,x + "gp_CI");
        data.return_levels = table(return_periods,gev,gum,gp,gev_CI,gum_CI,gp_CI);
        % exceedances per year
        year = datetime(h5read(file,x + "year"), 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
        gevX = h5read(file,x + "gev_exceedances_per_year");
        gumX = h5read(file,x + "gum_exceedances_per_year");
        gpX = h5read(file,x + "gp_exceedances_per_year");
        data.exceedances_per_year = array2timetable([gevX gumX gpX],'RowTimes',year, ...
            'VariableNames',["gev-" + string(return_periods) "gum-" + string(return_periods) ...
            "gp-" + string(return_periods)]);
    end

    if ~isempty(dataset)
        if istable(data) || istimetable(data)
            if istimetable(data) && strcmp(dataset,'time')
                data = data.Properties.RowTimes;
            else
                data = data(:,dataset);
            end
        elseif isstruct(data)
            fields = string(fieldnames(data));
            fields(ismember(fields,dataset)) = [];
            data = rmfield(data,fields);
            if length(fieldnames(data)) == 1
                data = data.(string(fieldnames(data)));
            end
        end
    end
end

function ht = harmonic2table(hin)
    c = hin(1,:);
    r1 = hin(2:end,1);
    hin = str2double(hin(2:end,2:end));
    ht = array2table(hin);
    ht = addvars(ht, r1, 'Before', 1);
    ht.Properties.VariableNames = c;
end

% fini