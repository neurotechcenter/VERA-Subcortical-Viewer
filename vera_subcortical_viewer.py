"""Build a local, interactive viewer from a FreeSurfer segmentation and electrode coordinates."""
import argparse
import csv
import html
import os
from pathlib import Path

import nibabel as nib
import numpy as np
import plotly.graph_objects as go
from scipy.ndimage import find_objects
from skimage.measure import marching_cubes


def read_lut(path):
    entries = {}
    if path:
        for line in path.read_text().splitlines():
            fields = line.partition("#")[0].split()
            if len(fields) >= 6:
                label, name, r, g, b = fields[:5]
                entries[int(label)] = (name, f"rgb({int(r)},{int(g)},{int(b)})")
    return entries


def contacts(path, transform):
    with path.open(encoding="utf-8-sig", newline="") as source:
        rows = csv.DictReader(source, delimiter="\t")
        if not {"name", "x", "y", "z"}.issubset(rows.fieldnames or []):
            raise ValueError("The electrode TSV needs name, x, y and z columns.")
        for line, row in enumerate(rows, 2):
            xyz = np.array([float(row[k]) for k in ("x", "y", "z")])
            if not row["name"].strip() or not np.isfinite(xyz).all():
                raise ValueError(f"Invalid electrode name or coordinate on TSV line {line}.")
            yield row["name"].strip(), nib.affines.apply_affine(transform, xyz)


def build(segmentation, electrode_file, space, lut):
    image = nib.load(str(segmentation))
    voxels = np.asanyarray(image.dataobj)
    if voxels.ndim != 3 or not np.issubdtype(voxels.dtype, np.integer) or voxels.min() < 0:
        raise ValueError("Expected a 3-D segmentation with nonnegative integer labels.")
    affine = image.affine
    if space == "tkr":
        if not hasattr(image.header, "get_vox2ras_tkr"):
            raise ValueError("Use the matching FreeSurfer MGZ/MGH for tkregister RAS coordinates.")
        transform = affine @ np.linalg.inv(image.header.get_vox2ras_tkr())
    else:
        transform = np.eye(4)
    # Compact arbitrary atlas IDs before computing bounding boxes.
    ids, dense = np.unique(voxels, return_inverse=True)
    boxes = find_objects(dense.reshape(voxels.shape).astype(np.int32) + 1)
    fig = go.Figure()
    groups = []
    for label, box in zip(ids, boxes):
        if label == 0:
            continue
        mask = np.pad(voxels[box] == label, 1)
        vertices, faces, _, _ = marching_cubes(mask.astype(np.float32), level=0.5)
        origin = np.array([axis.start for axis in box]) - 1
        xyz = nib.affines.apply_affine(affine, vertices + origin)
        name, color = lut.get(int(label), (f"Label {label}", None))
        fig.add_mesh3d(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2],
                       i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
                       color=color, opacity=0.4, name=f"{name} [{label}]",
                       hoverinfo="name", visible=False)
        groups.append("Surfaces")
    for name, xyz in contacts(electrode_file, transform):
        fig.add_scatter3d(x=[xyz[0]], y=[xyz[1]], z=[xyz[2]], mode="markers",
                          marker=dict(size=4, color="#a52c43"), name=name,
                          hoverinfo="name+x+y+z", visible=True)
        groups.append("Electrodes")
    fig.update_layout(showlegend=False, margin=dict(t=0, b=0, l=0, r=0),
                      scene=dict(aspectmode="data", xaxis_title="R (mm)",
                                 yaxis_title="A (mm)", zaxis_title="S (mm)"),
                      uirevision="keep-camera")
    return fig, groups


def write_viewer(fig, groups, path):
    panels = []
    for group in ("Surfaces", "Electrodes"):
        controls = []
        for index, (trace, category) in enumerate(zip(fig.data, groups)):
            if category == group:
                checked = " checked" if trace.visible else ""
                controls.append(f'<label><input type="checkbox" data-index="{index}"{checked}>'
                                f'{html.escape(trace.name)}</label>')
        panels.append(f"<section><h2>{group}</h2>{''.join(controls)}</section>")
    plot = fig.to_html(full_html=False, include_plotlyjs=True, div_id="scene",
                       config={"displaylogo": False, "responsive": True})
    path.write_text("""<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>VERA Subcortical Viewer</title><style>
*{box-sizing:border-box}body{margin:0;display:flex;height:100vh;font:15px Arial,sans-serif;color:#183441}
aside{width:clamp(200px,26vw,340px);overflow:auto;padding:20px;background:#f3f6f7}
h1{font-size:23px}h2{font-size:16px}label{display:block;padding:3px}main{flex:1;min-width:0}
input[type=search]{width:100%;padding:7px} input{accent-color:#126d75}
</style><aside><h1>VERA Subcortical Viewer</h1>
<input type="search" aria-label="Filter structures and contacts" placeholder="Filter structures and contacts">
""" + "".join(panels) + "</aside><main>" + plot + """</main><script>
document.querySelector('aside').addEventListener('change', e => {
  if (e.target.type === 'checkbox') Plotly.restyle('scene', {visible:e.target.checked}, [+e.target.dataset.index]);
});
document.querySelector('input[type=search]').addEventListener('input', e => {
  for (const item of document.querySelectorAll('label')) {
    item.style.display = item.textContent.toLowerCase().includes(e.target.value.toLowerCase()) ? '' : 'none';
  }
});
</script></html>""", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("segmentation", type=Path)
    parser.add_argument("electrodes", type=Path, help="TSV with name, x, y, z (mm)")
    parser.add_argument("--space", required=True, choices=["scanner", "tkr"])
    default_lut = Path(os.environ.get("FREESURFER_HOME", ".")) / "FreeSurferColorLUT.txt"
    parser.add_argument("--lut", type=Path, default=default_lut if default_lut.is_file() else None)
    parser.add_argument("--out", type=Path, default=Path("vera-viewer.html"))
    args = parser.parse_args()
    fig, groups = build(args.segmentation, args.electrodes, args.space, read_lut(args.lut))
    write_viewer(fig, groups, args.out)
    print(f"{groups.count('Surfaces')} surfaces, {groups.count('Electrodes')} electrodes: {args.out}")


if __name__ == "__main__":
    main()
