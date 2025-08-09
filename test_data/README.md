# Test Data

This directory contains small test datasets for unit testing.

## Files

- `test_system.psf`: Topology file (61.57 MB)
- `test_trajectory.xtc`: Trajectory with 8 frames (15.66 MB)

## Generation

These files were generated from the full trajectory using `create_test_data.py`.
- Original PSF: step5_assembly.psf
- Original XTC: md_wrapped.xtc
- Original trajectory: 88799 frames
- Test trajectory: frames 60000 to 60140 (step 20)

## Usage

```python
import MDAnalysis as mda
u = mda.Universe('test_data/test_system.psf', 'test_data/test_trajectory.xtc')
print(f'Loaded {len(u.atoms)} atoms, {len(u.trajectory)} frames')
```
