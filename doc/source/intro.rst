Introduction
------------

PySUNDIALS is a python package providing python bindings for the SUNDIALS suite
of solvers. It is being developed by the triple 'J' group at Stellenbosch
University, South Africa. While python bindings for SUNDIALS will hopefully be
generally useful in the computational scientific community, they are being
developed with the specific aim of providing a robust underlying numerical
solver capable of implementing models conforming to the Systems Biology Markup
Language specification (version 2), including triggers, events, and delays.  As
such the development process is partially driven by the continuing parallel
development of PySCeS.

This documentation serves to introduce and act as a reference for PySUNDIALS. As
such it focuses on the differences in behaviour between PySUNDIALS and SUNDIALS.
If you wish to know more about the underlying mathematical considerations, or
wish a more in depth discussion of how to go about using SUNDIALS in general,
the SUNDIALS documentation is the place to start. We assume a basic knowledge of
initial value problems and their computational solutions on the part of the
reader.
