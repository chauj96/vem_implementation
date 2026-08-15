import sys
import pyvista as pv
import numpy as np
from scipy.io import savemat

"""
Read a 3D polyhedral VTU mesh and reconstruct the geometric data required
by the MFD/VEM discretization.

The VTU file provides the mesh vertices and polyhedral connectivity. From
this information, the reader constructs a unique global face structure and
computes the geometric quantities needed by the discretization, including
face areas, face centroids, unit normals, cell volumes and cell centroids.

Face geometry is reconstructed by triangulating each polygonal face about
its vertex mean. Cell volumes and centroids are then computed from the
oriented boundary-face triangulations.

As a consistency check, the reconstructed geometry is verified against the
discrete geometric identity

    C^T N = |E| I,

where the rows of C contain the vectors from the cell centroid to the face
centroids and the rows of N contain the outward face area vectors. This
identity is required for exact discrete divergence of linear fields and is
a fundamental geometric consistency condition for the MFD/VEM formulation.
"""

input_file = sys.argv[1]
output_file = sys.argv[2]

mesh = pv.read(input_file)

face_dict = {}
face_struct = []
cell_struct = []

# construct global face_struct and cell connectivity
for cid in range(mesh.n_cells):

    cell = mesh.get_cell(cid)
    cell_center_guess = np.mean(cell.points, axis=0)

    cell_faces = []
    cell_orientations = []
    cell_face_normals = []

    for face in cell.faces:

        point_ids = tuple(face.point_ids)
        face_key = tuple(sorted(point_ids))

        verts = np.asarray(face.points, dtype=float)
        face_vertex_mean = np.mean(verts, axis=0)

        area_total = 0.0
        centroid_accum = np.zeros(3)
        nhat_total = np.zeros(3)
        triangles = []

        for k in range(len(verts)):

            a = face_vertex_mean
            b = verts[k]
            c = verts[(k+1) % len(verts)]

            nhat = np.cross(b-a, c-a)
            Atri = 0.5 * np.linalg.norm(nhat)

            if Atri < 1e-14:
                continue

            Ctri = (a+b+c)/3.0

            area_total += Atri
            centroid_accum += Atri * Ctri
            nhat_total += nhat

            triangles.append((a, b, c, nhat))

        if area_total < 1e-14:

            print("\nDegenerate face")
            print("cell =", cid)
            print("point_ids =", point_ids)
            print("nverts =", len(verts))
            print("verts =\n", verts)

            rank = np.linalg.matrix_rank(verts - verts[0])

            print("rank =", rank)
            print("area_total =", area_total)

            continue

        face_center = centroid_accum / area_total
        face_normal = nhat_total / np.linalg.norm(nhat_total)
        face_area = area_total

        # store global face once
        if face_key not in face_dict:

            fid = len(face_struct)
            face_dict[face_key] = fid

            face_struct.append({
                "verts": point_ids,
                "center": face_center,
                "normal": face_normal,
                "area": face_area,
                "triangles": triangles,
                "cells": [cid]
            })

        else:

            fid = face_dict[face_key]

            if cid not in face_struct[fid]["cells"]:
                face_struct[fid]["cells"].append(cid)

        # orientation for this cell
        nf = face_struct[fid]["normal"]
        Cf = face_struct[fid]["center"]

        signf = np.sign(np.dot(Cf - cell_center_guess, nf))

        if signf == 0:
            raise RuntimeError(f"Zero orientation in cell {cid}, face {fid}")

        cell_faces.append(fid)
        cell_orientations.append(signf)
        cell_face_normals.append(signf * nf)

    cell_struct.append({
        "center": None,
        "volume": None,
        "faces": cell_faces,
        "faces_orientation": cell_orientations,
        "face_normals": np.asarray(cell_face_normals)
    })

# reconstruct cell volume and centroid
for cid in range(len(cell_struct)):

    V_total = 0.0
    M = np.zeros(3)

    face_ids = cell_struct[cid]["faces"]
    signs = cell_struct[cid]["faces_orientation"]

    for fid, signf in zip(face_ids, signs):

        for a, b, c, nhat in face_struct[fid]["triangles"]:

            nhat_out = signf * nhat

            # Volume: V = 1/6 sum a · nhat
            V_total += np.dot(a, nhat_out) / 6.0

            # Centroid numerator
            for d in range(3):

                M[d] += nhat_out[d] * (
                    (a[d] + b[d])**2
                    + (b[d] + c[d])**2
                    + (c[d] + a[d])**2
                ) / 48.0

    if abs(V_total) < 1e-14:
        raise RuntimeError(f"Zero cell volume at cell {cid}")

    cell_struct[cid]["volume"] = abs(V_total)
    cell_struct[cid]["center"] = M / V_total

# CTN CHECK
# print("\n=== CTN CHECK ===")

errs = []

for cid in range(len(cell_struct)):

    Cc = cell_struct[cid]["center"]
    Vc = cell_struct[cid]["volume"]

    face_ids = cell_struct[cid]["faces"]
    signs = cell_struct[cid]["faces_orientation"]

    C = np.zeros((len(face_ids), 3))
    N = np.zeros((len(face_ids), 3))

    for k, (fid, signf) in enumerate(zip(face_ids, signs)):

        Cf = face_struct[fid]["center"]
        nf = face_struct[fid]["normal"]
        Af = face_struct[fid]["area"]

        C[k, :] = Cf - Cc
        N[k, :] = Af * signf * nf

    CTN = C.T @ N
    err = np.linalg.norm(CTN - Vc*np.eye(3), ord="fro") / max(abs(Vc), 1e-14)

    errs.append(err)

tol = 1e-10
max_err = np.max(errs)

print(f"Discrete divergence identity error min  = {np.min(errs):.6e}")
print(f"Discrete divergence identity error mean = {np.mean(errs):.6e}")
print(f"Discrete divergence identity error max  = {max_err:.6e}")

if max_err < tol:
    print(f"Geometric consistency check: PASS (tol = {tol:.0e})")
    print("Mesh geometry is consistent with the MFD/VEM discrete divergence identity.")
else:
    print(f"Geometric consistency check: WARNING (tol = {tol:.0e})")
    print("Mesh geometry does not satisfy the MFD/VEM discrete divergence identity within tolerance.")

nCells = len(cell_struct)
nFaces = len(face_struct)

print("\n3D mesh model info:")
print(f"  {mesh.n_points} vertices")
print(f"  {nCells} cells")
print(f"  {nFaces} faces")

# MATLAB STYLE EXPORT!

V3 = np.asarray(mesh.points, dtype=np.float64)

# face_struct
face_dtype = [
    ("cells", "O"),
    ("verts", "O"),
    ("center", "O"),
    ("normal", "O"),
    ("area", "O")
]

face_struct_matlab = np.empty((len(face_struct), 1), dtype=face_dtype)

for f, fs in enumerate(face_struct):

    face_struct_matlab[f, 0]["cells"] = \
        np.asarray([c + 1 for c in fs["cells"]], dtype=np.float64).reshape(1, -1)

    face_struct_matlab[f, 0]["verts"] = \
        np.asarray([v + 1 for v in fs["verts"]], dtype=np.float64).reshape(1, -1)

    face_struct_matlab[f, 0]["center"] = \
        np.asarray(fs["center"], dtype=np.float64).reshape(3, 1)

    face_struct_matlab[f, 0]["normal"] = \
        np.asarray(fs["normal"], dtype=np.float64).reshape(3, 1)

    face_struct_matlab[f, 0]["area"] = \
        np.asarray([[fs["area"]]], dtype=np.float64)


# cell_struct
cell_dtype = [
    ("center", "O"),
    ("volume", "O"),
    ("faces", "O"),
    ("faces_orientation", "O"),
    ("face_normals", "O")
]

cell_struct_matlab = np.empty((len(cell_struct), 1), dtype=cell_dtype)

for c, cs in enumerate(cell_struct):

    cell_struct_matlab[c, 0]["center"] = \
        np.asarray(cs["center"], dtype=np.float64).reshape(3, 1)

    cell_struct_matlab[c, 0]["volume"] = \
        np.asarray([[cs["volume"]]], dtype=np.float64)

    cell_struct_matlab[c, 0]["faces"] = \
        np.asarray([f + 1 for f in cs["faces"]], dtype=np.float64).reshape(1, -1)

    cell_struct_matlab[c, 0]["faces_orientation"] = \
        np.asarray(cs["faces_orientation"], dtype=np.float64).reshape(1, -1)

    cell_struct_matlab[c, 0]["face_normals"] = \
        np.asarray(cs["face_normals"], dtype=np.float64)


# cells3D
cells3D_matlab = np.empty((len(cell_struct), 1), dtype=object)

for c, cs in enumerate(cell_struct):

    vertex_ids = []

    for fid in cs["faces"]:
        vertex_ids.extend(face_struct[fid]["verts"])

    vertex_ids = np.unique(vertex_ids)

    cells3D_matlab[c, 0] = \
        np.asarray(vertex_ids + 1, dtype=np.float64).reshape(1, -1)


# save
savemat(
    output_file,
    {
        "V3": V3,
        "face_struct": face_struct_matlab,
        "cell_struct": cell_struct_matlab,
        "cells3D": cells3D_matlab
    }
)