# VEM Implementation

MATLAB implementation of a mixed and hybridized Virtual Element Method (VEM) for three-dimensional linear elasticity on general polyhedral meshes. The code provides several options for **(1) the VEM formulation**, **(2) mesh geometry** and **(3) linear solver** as described below.

## Formulations
#### Mixed VEM

The mixed formulation assembles and solves the global saddle-point system for stress and displacement:

```math
\begin{bmatrix}
K & B^T \\
B & 0
\end{bmatrix}
\begin{bmatrix}
\sigma \\
u
\end{bmatrix}
=
\begin{bmatrix}
g \\
-f
\end{bmatrix}.
```

#### Hybridized VEM

The hybridized formulation treats element stress unknowns independently and enforces traction continuity across interior faces through Lagrange multipliers:

```math
\begin{bmatrix}
\widehat{K} & \widehat{B}^T & -C^T \\
\widehat{B} & 0 & 0 \\
C & 0 & 0
\end{bmatrix}
\begin{bmatrix}
\widehat{\sigma} \\
u \\
\lambda
\end{bmatrix}
=
\begin{bmatrix}
g \\
-f \\
0
\end{bmatrix}.
```

The element-local stress and displacement unknowns are eliminated, reducing the global problem to

$$
H\lambda = h,
\qquad
H \in \mathbb{R}^{6N_f^{\mathrm{int}} \times 6N_f^{\mathrm{int}}}.
$$

Once $\lambda$ is obtained, $(\sigma_E,u_E)$ are recovered independently on each element.

The formulation can be selected in `main.m` using

```matlab
formOption = 'mixed';
formOption = 'hybrid';
```

## Mesh Geometry

Three 3D mesh models are provided by default with mesh data available in both MATLAB `.mat` and `.vtu` formats:

1. **Unit cube**  
   A structured Cartesian mesh generated directly using MRST.

2. **Two-fault model**  
   A faulted polyhedral mesh.

3. **Polyhedral Voronoi model**  
   A general polyhedral Voronoi mesh.

The provided `.mat` files contain preprocessed mesh geometry and can be used directly in MATLAB. The `.vtu` files can be used to reconstruct the geometry through the VTU reader.

### Reading Custom VTU Meshes

The VTU reader is also supported for users who want to run the VEM implementation on their own 3D polyhedral meshes. The reader uses **PyVista** to extract the polyhedral topology and reconstruct the geometric quantities required by the VEM formulation including

- global faces and cell-to-face connectivity
- face areas, centroids, and normals
- cell volumes and centroids
- local face orientations

The reconstructed geometry is additionally checked against the discrete divergence consistency identity which provides a numerical consistency check for the geometric quantities used by the MFD/VEM discretization.

### Python Environment Setup

So the VTU reader requires Python with `PyVista`, `NumPy`, and `SciPy`. A Conda environment file is provided with the repository.

Clone the repository and enter the project directory:

```bash
git clone https://github.com/chauj96/vem_implementation.git
cd vem_implementation
```

Create the provided Conda environment:

```bash
conda env create -f environment.yml
```

Confirm that the `vem` environment was created:

```bash
conda env list
```

You should see an environment named `vem`.

Activate the environment:

```bash
conda activate vem
```

Verify that the required Python packages are available:

```bash
python -c "import pyvista, numpy, scipy; print('VEM environment ready')"
```

If the environment is configured correctly, this should print

```text
VEM environment ready
```

The reader automatically searches for the `vem` Conda environment and uses it to execute the Python mesh reader. 
Once the environment is ready, a custom VTU mesh can be loaded directly in `main.m` using

```matlab
filename = 'meshes/your_mesh.vtu';
[cell_struct, face_struct, V3, cells3D] = readVTU(filename);
```

The mesh geometry is then automatically reconstructed and the user is notified whether the mesh is suitable for use with the VEM discretization.

## Linear Solver and Running the Code

Once the mesh and VEM formulation are selected, the corresponding linear solver can be specified. The available solver options for each formulation are summarized below.

| Formulation | Solver | Description |
|---|---|---|
| **Mixed VEM** | **Direct** | Direct solution of the global saddle-point system. |
|  | **Iterative** | GMRES with a block Schur complement preconditioner. |
| **Hybridized VEM** | **PCG** | Automatically used for the reduced multiplier system $H\lambda=h$, followed by element-local recovery. |

When the options are fully configured, run

```matlab
main
```
