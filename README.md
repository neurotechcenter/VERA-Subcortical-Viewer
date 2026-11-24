# VERA Subcortical Viewer

Explore electrode contacts and labeled brain structures together in an interactive 3-D view. The viewer reads a FreeSurfer segmentation directly. It generates a standalone HTML file that opens locally in a browser.

## MATLAB desktop viewer

Open this repository in MATLAB R2021a or newer. Add the FreeSurfer MATLAB utilities to the path, then launch the viewer:

```matlab
addpath(fullfile(getenv('FREESURFER_HOME'), 'matlab'))
vera_subcortical_viewer
```

Choose an MGZ/MGH segmentation, a contact TSV, and the contact coordinate convention. Check a surface to draw it, or uncheck individual contacts to hide them. Drag to rotate, and use data tips to identify contacts. Surfaces are generated on first selection. To open a known dataset directly:

```matlab
vera_subcortical_viewer('aparc+aseg.mgz', 'electrodes.tsv', 'scanner')
```

Use an optional fourth argument for a custom FreeSurfer color lookup table. The MATLAB and Python viewers share the input conventions below.

## Python standalone HTML viewer

Install Python dependencies with `python -m pip install nibabel numpy scipy scikit-image plotly` (Plotly 5.5 or newer). From the repository folder, run:

```sh
python vera_subcortical_viewer.py /path/to/subject/mri/aparc+aseg.mgz electrodes.tsv --space scanner --out viewer.html
```

Open `viewer.html`. Click the checkboxes to show structures and electrode contacts together. Type part of a name to filter the lists. Drag the scene to rotate and scroll to zoom. Contacts start visible; surfaces start hidden. The output embeds the selected subject's geometry and all supplied coordinates.

## Inputs and coordinate system

- **Segmentation:** A 3-D integer-label volume such as `aseg.mgz` or `aparc+aseg.mgz`. Each nonzero integer becomes a selectable surface. For individual thalamic nuclei, supply a segmentation that distinguishes those nuclei.
- **Label names and colors:** Set `FREESURFER_HOME` or supply `--lut /path/to/FreeSurferColorLUT.txt`. Use the table corresponding to your segmentation. Unlisted labels appear by numeric ID.
- **Electrode positions:** A tab-separated text file with a header and one contact per row. The coordinates below illustrate the format.

```text
name	x	y	z
A1	12.0	-15.0	8.0
A2	13.5	-14.0	9.0
```

Use `--space scanner` for scanner RAS millimeters in the segmentation's physical space. Use `--space tkr` for tkregister (surface) RAS millimeters referenced to the **same voxel geometry** as the input MGZ/MGH. A cropped or resampled segmentation can have a different tkregister frame. Register CT, template-space, or other coordinates to the subject first. The viewer does not estimate registration.

## Original VERA / MATLAB workflow

For a VERA export, the [LabelVolume2Surface component](LabelVolume2Surface%20component) and [MatOutput component](matoutput%20component) illustrate how to save a labeled surface. The [original MATLAB example](example_plot_3D) selects three thalamic labels. The Python viewer can be run independently of that workflow.
