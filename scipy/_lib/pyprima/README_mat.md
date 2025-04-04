This is the README file for using PRIMA under MATLAB.

## Prerequisites

PRIMA supports MATLAB R2018a and later releases. To use PRIMA, you need first
set up the [MEX](https://www.mathworks.com/help/matlab/ref/mex.html) of your MATLAB.
**The setup of MEX is a pure MATLAB usage problem and it has nothing to do with PRIMA.**

To see whether your MEX is ready, run the following code in MATLAB:

```matlab
mex('-setup', '-v', 'fortran'); mex('-v', fullfile(matlabroot, 'extern', 'examples', 'refbook', 'timestwo.F'));
```

If this completes successfully, then your MEX is ready. Otherwise, it is not, and
you may try the [`setup_mex` package](https://github.com/equipez/setup_mex) at
```bash
https://github.com/equipez/setup_mex
```
It will help you to set MEX up on Windows or macOS (the setup of MEX is trivial on Linux).
In case `setup_mex` does not work, you need to consult a local MATLAB expert or the technical support of
MathWorks about "[how to set up MEX](https://www.mathworks.com/help/matlab/ref/mex.html)", which is
**not** part of PRIMA.


## Installation

Download and decompress the source code package of PRIMA,
or clone it from [GitHub](https://github.com/libprima/prima) or [Gitee](https://gitee.com/zaikunzhang/prima).
You will obtain a folder containing setup.m. Place this folder at the location where you
want PRIMA to be installed. In MATLAB, change the directory to this folder,
and execute the following command:

```matlab
setup
```

If this command runs successfully, PRIMA is installed. You may execute the
following command in MATLAB to verify the installation:

```matlab
testprima
```


## Usage

PRIMA provides a MATLAB function `prima`, which can solve general
constrained or unconstrained optimization problems without using derivatives.

The `prima` function can automatically identify the type of your problem
and then call one of Powell's solvers, namely COBYLA, UOBYQA, NEWUOA, BOBYQA,
and LINCOA. The user can also specify the solver by setting the `solver` field
of the options passed to `prima`.

The `prima` function is designed to be compatible with the `fmincon`
function available in the Optimization Toolbox of MATLAB. You can call `prima`
in exactly the same way as calling `fmincon`.

For detailed syntax of `prima`, use the standard `help` command of MATLAB:

```matlab
help prima
```



## Uninstall

PRIMA can be uninstalled using the setup.m script by executing the following
command in MATLAB:

```matlab
setup uninstall
```

## Contact

In case of problems, [open a GitHub issue](https://github.com/libprima/prima/issues) or [contact
Zaikun Zhang](https://www.zhangzk.net).
