# VEM Implementation

MATLAB implementation of a mixed Virtual Element Method (VEM) for three-dimensional linear elasticity on general polyhedral meshes. The code supports simulations on general 3D polyhedral meshes. Mesh geometry can either be loaded from preprocessed MATLAB `.mat` files or reconstructed directly from `.vtu` files.

## Available Meshes

Three example mesh models are currently included:

1. **Unit cube**  
   A structured Cartesian mesh generated using MRST.

2. **Two-fault model**  
   A faulted polyhedral mesh. Preprocessed mesh geometry is provided as a MATLAB `.mat` file.

3. **Polyhedral Voronoi model**  
   A general polyhedral Voronoi mesh. Preprocessed mesh geometry is also provided as a MATLAB `.mat` file.

The two-fault and polyhedral Voronoi examples can therefore be run directly from the provided `.mat` files without using the Python-based VTU reader.

## Reading Custom VTU Meshes

The code also provides a general VTU mesh reader for users who want to run the VEM implementation on their own polyhedral meshes.

The VTU reader uses **PyVista** to extract the polyhedral topology and reconstructs the geometric quantities required by the MFD/VEM formulation, including

- global faces and cell-to-face connectivity,
- face areas, centroids, and normals,
- cell volumes and centroids,
- local face orientations.

The reconstructed geometry is additionally checked against the discrete divergence consistency identity which provides a numerical consistency check for the geometric quantities used by the MFD/VEM discretization.

### Python Environment Setup

The VTU reader requires Python with `PyVista`, `NumPy`, and `SciPy`. A Conda environment file is provided with the repository.

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

The MATLAB VTU reader automatically searches for the `vem` Conda environment and uses it to execute the Python mesh reader.

## Running the Code

Open the repository in MATLAB and run

```matlab
main
```

The mesh and solver options (direct or iterative) can be selected near the beginning of `main.m`.

For a custom VTU mesh, specify the mesh file and load it using

```matlab
filename = 'meshes/your_mesh.vtu';
[cell_struct, face_struct, V3, cells3D] = readVTU(filename);
```

The mesh geometry will be reconstructed and checked before being used by the VEM implementation.
