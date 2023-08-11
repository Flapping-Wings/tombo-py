This test suite compares the calculations from the original MATLAB implementation with those from `tombo-py`.
**This is admittedly not great design, but it is what we are using for now.**

The data from the MATLAB implementation is in the `matlab_data` folder.
This data was saved using MATLAB's `save()` function.
The data files named after functions were saved after their respective function calls,
and the `loop_data` files were taken at the end of the velocity calculations in the main simulation loop
(**BEFORE** the call to `tbshedB`).
