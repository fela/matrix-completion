Matrix Completion
=================

This code implements a matrix completion and decomposition algorithm, based on [this][1] paper.

[1]:http://arxiv.org/abs/1102.4807

The main method you will want to use is `matrix_decomposition()` which decomposes an optionally incomplete matrix, into two components, one low rank, and one sparse. It takes the following arguments:

* `Y`: the matrix to decompose. It is supposed that **`Y` is equal to `TH + GA + W`**, where `TH` is an approximately low rank matrix, `GA` is a sparse "spiky" matrix, and `W` is noise. A full matrix needs to be given, but some parts can be ignored, as specified by the mask.
* `Mask`: a mask can be given to treat the matrix as incomplete. It must be of the same dimentions as `Y`. It must have `0` or `False` in the positions where `Y` is incomplete, `1` or `True` elsewhere. Defaults to `None`.
* `lambda_d`: regularization paramater for the low rank `TH` matrix.
* `mu_d`: regularization parameter for the sparse `GA` matrix. Use higher values if no spikes are expected.
* `alpha`: parameter that limits the maximum element of the low rank TH matrix. Bigger matrices will need bigger alpha` values.

The method will **return** `TH` and `GA`.

Additionally the code contains some methods and classes to generate synthetic matrices and to test the decomposition in various ways.

