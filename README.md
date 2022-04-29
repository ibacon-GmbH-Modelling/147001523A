Folder with BYOM and DEBtox2019 code used for the analysis of C. Dubia data.
See copyright notices in subfolders.

The original MATLAB implementation was expanded with the use of a C++ code to
speed up the solution of the ODEs.

This implementation uses the boost::odeint libraries, available at the following 
address:

[https://www.boost.org/](https://www.boost.org/)

(the version used in this project is v1.77.0)

The C++ code has to be compiled within MATLAB using the MEX C++ APIs
linking the external boost libraries:

```
>> mex COMPFLAGS='$COMPFLAGS -std=c++11' test_derivatives.cpp -I<path to boost libraries>
```

The original MATLAB code can still be run by substituting the file
`call_deri.m` with `call_deri_old.m`
in the DEBtox_2019_v45a folder.