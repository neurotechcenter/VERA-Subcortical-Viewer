function fig = vera_subcortical_viewer(segFile, electrodeFile, space, lutFile, smoothSurfaces)
% Interactively display a FreeSurfer segmentation and implanted contacts.
% Run without arguments for file pickers, or supply paths and scanner/tkr.
% Requires MATLAB R2021a or later and FreeSurfer's MATLAB directory on path.
if nargin < 5, smoothSurfaces = true; end
fig = [];
if nargin < 1
    [f,p] = uigetfile({'*.mgz;*.mgh','FreeSurfer segmentation'});
    if isequal(f,0), return; end
    segFile = fullfile(p,f);
end
if nargin < 2
    [f,p] = uigetfile('*.tsv','Electrode coordinates (name, x, y, z)');
    if isequal(f,0), return; end
    electrodeFile = fullfile(p,f);
end
if nargin < 3
    space = questdlg('Which coordinate system are the electrodes in?', ...
        'Electrode coordinates','Scanner RAS','tkregister RAS','Cancel','Cancel');
    if isempty(space) || strcmp(space,'Cancel'), return; end
    if strcmp(space,'Scanner RAS'), space = 'scanner'; else, space = 'tkr'; end
end
space = validatestring(space,{'scanner','tkr'});
if nargin < 4
    lutFile = fullfile(getenv('FREESURFER_HOME'),'FreeSurferColorLUT.txt');
    if ~isfile(lutFile), lutFile = ''; end
end
if exist('MRIread','file') ~= 2
    error('Add the matlab folder of your FreeSurfer installation to the MATLAB path.');
end
[~,~,ext] = fileparts(segFile);
if ~ismember(lower(string(ext)),[".mgz",".mgh"])
    error('Select a FreeSurfer MGZ or MGH segmentation.');
end
mri = MRIread(char(segFile),0,1); % Rows, columns, slices; affine expects column, row, slice.
if isempty(mri), error('FreeSurfer could not read the segmentation.'); end
vol = mri.vol;
ids = double(unique(vol(:)));
if ndims(vol) ~= 3 || any(~isfinite(ids) | ids < 0 | ids ~= fix(ids))
    error('Expected a 3-D nonnegative integer segmentation.');
end
ids(ids == 0) = [];
[names,colors] = labelInfo(ids,lutFile);
opts = detectImportOptions(electrodeFile,'FileType','text','Delimiter','\t');
if ~all(ismember({'name','x','y','z'},opts.VariableNames))
    error('The electrode TSV must include name, x, y and z columns.');
end
opts = setvartype(opts,'name','string');
opts = setvartype(opts,{'x','y','z'},'double');
contacts = readtable(electrodeFile,opts);
xyz = contacts{:,{'x','y','z'}};
if any(~isfinite(xyz(:))) || any(ismissing(contacts.name) | strlength(strtrim(contacts.name)) == 0)
    error('Every contact needs a name and three finite coordinates in millimeters.');
end
if strcmp(space,'tkr')
    xyz = transformPoints(mri.vox2ras0 / mri.tkrvox2ras, xyz);
end

fig = uifigure('Name','VERA Subcortical Viewer','Position',[100 100 980 650]);
layout = uigridlayout(fig,[1 2]);
layout.ColumnWidth = {260,'1x'};
side = uigridlayout(layout,[4 1]);
side.RowHeight = {24,'2x',24,'1x'};
uilabel(side,'Text','Surfaces','FontWeight','bold');
surfaces = uitree(side,'checkbox');
uilabel(side,'Text','Electrode contacts','FontWeight','bold');
electrodes = uitree(side,'checkbox');
ax = uiaxes(layout);
hold(ax,'on'); axis(ax,'equal'); grid(ax,'on'); view(ax,3);
xlabel(ax,'R (mm)'); ylabel(ax,'A (mm)'); zlabel(ax,'S (mm)');
ax.Interactions = [rotateInteraction zoomInteraction panInteraction];
axtoolbar(ax,{'rotate','pan','zoomin','zoomout','restoreview'});
meshes = gobjects(numel(ids),1);
for k = 1:numel(ids)
    uitreenode(surfaces,'Text',char(names(k)),'NodeData',k);
end
contactNodes = gobjects(height(contacts),1);
for k = 1:height(contacts)
    contactNodes(k) = uitreenode(electrodes,'Text',char(contacts.name(k)),'NodeData',k);
end
points = scatter3(ax,xyz(:,1),xyz(:,2),xyz(:,3),38,[0.65 0.17 0.26],'filled');
points.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Contact',contacts.name);
electrodes.CheckedNodes = contactNodes;
surfaces.CheckedNodesChangedFcn = @updateSurfaces;
electrodes.CheckedNodesChangedFcn = @updateContacts;

    function updateSurfaces(~,event)
        selected = false(numel(ids),1);
        for node = reshape(event.CheckedNodes,1,[])
            selected(node.NodeData) = true;
        end
        for j = 1:numel(ids)
            if selected(j) && ~isgraphics(meshes(j))
                idx = find(vol == ids(j));
                [r,c,s] = ind2sub(size(vol),idx);
                lo = [min(c) min(r) min(s)];
                crop = vol(min(r):max(r),min(c):max(c),min(s):max(s)) == ids(j);
                mask = zeros(size(crop,1)+2,size(crop,2)+2,size(crop,3)+2,'single');
                mask(2:end-1,2:end-1,2:end-1) = crop;
                fv = isosurface(mask,0.5);
                % Local padded vertex 2 corresponds to global zero-based lo-1.
                fv.vertices = transformPoints(mri.vox2ras0,double(fv.vertices)+lo-3);
                if smoothSurfaces
                    fv.vertices = smoothDisplayMesh(fv.vertices,fv.faces);
                end
                meshes(j) = patch(ax,'Faces',fv.faces,'Vertices',fv.vertices, ...
                    'FaceColor',colors(j,:),'FaceAlpha',0.4,'EdgeColor','none');
            elseif isgraphics(meshes(j))
                meshes(j).Visible = matlab.lang.OnOffSwitchState(selected(j));
            end
        end
    end

    function updateContacts(~,event)
        selected = false(height(contacts),1);
        for node = reshape(event.CheckedNodes,1,[])
            selected(node.NodeData) = true;
        end
        % Retain contact indices for the data tips while hiding unchecked points.
        shown = xyz;
        shown(~selected,:) = NaN;
        set(points,'XData',shown(:,1),'YData',shown(:,2),'ZData',shown(:,3));
    end
end

function xyz = transformPoints(matrix,xyz)
xyz = xyz * matrix(1:3,1:3).' + matrix(1:3,4).';
end

function [names,colors] = labelInfo(ids,lutFile)
names = compose("Label %d",ids);
colors = lines(numel(ids));
if strlength(string(lutFile)) == 0, return; end
rows = readlines(lutFile);
for k = 1:numel(rows)
    fields = split(strtrim(extractBefore(rows(k) + "#","#")));
    if numel(fields) < 6, continue; end
    [matched,j] = ismember(str2double(fields(1)),ids);
    if matched
        names(j) = compose("%s [%d]",fields(2),ids(j));
        colors(j,:) = reshape(str2double(fields(3:5)),1,3)/255;
    end
end
end

function vertices = smoothDisplayMesh(vertices,faces)
% Taubin smoothing changes display vertices only.
edges = [faces(:,[1 2]); faces(:,[2 3]); faces(:,[3 1])];
n = size(vertices,1);
adj = sparse(edges(:,1),edges(:,2),1,n,n);
adj = spones(adj + adj.');
weights = spdiags(1 ./ max(sum(adj,2),1),0,n,n) * adj;
for pass = 1:5
    for step = [0.5 -0.53]
        vertices = vertices + step * (weights*vertices - vertices);
    end
end
end
