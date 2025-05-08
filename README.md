<div align="center">

# QSLD

## Quantum Simulation Library in D

</div>

### Introduction

QSLD is a quantum simulation library, mainly for my own benefit in learning more about quantum computing and the technical details behind what is happening at a mathematical level and also to learn the D programming language. The library should be used at your own risk knowing that it is not optimized to the maximum extent nor is it intended to be super production ready. If you would like to use a quantum simulation library which is production ready and optimized please use Qiskit instead of QSLD. I will provide ample documentation for anyone who would like to use this library for recreational purposes. To access it, refer to the wiki (not available yet). This library embraces the idea of simplicity and therefore you will have to build the entire library from scratch into a dynamically linked library file to then link with your main D program. Instructions for doing all of this are below.

### Dependencies For Building

- `dmd` or `gdc` or `ldc2`
- `rdmd`
- `make` (optional)

### Building

Before building, clone the project and cd into the root.

If you are using dmd and do not want to compile any examples in the `examples/` directory, You may use the following command:

```console
$ ./build.d 
```
If you are using another compiler (i.e. gdc or ldc2), you must specify the compiler name like so:

```console
$ ./build.d -c gdc
```
If you would like to compile an example, you must specify an example name to compile:

```console
$ ./build.d -e example_name
```
you can also do multiple in a row:

```console
$ ./build.d -e example_name1 -e example_name2
```

### Linking

To link the `libqsld.a` file with your main D program, either put it in a global directory for libraries like `/usr/local/lib` or `/usr/lib` and use the following command:

```console
$ <compiler-name> -L="-lqsld" -of=<bin-name> <d-filename>.d
```

or keep it in the root of the project and use the following:

```console
$ <compiler-name> -L-L/path/to/project/libsqld.a -L="-lqsld" -of=<bin-name> <d-filename>.d

```

Linking with any of the examples in `examples/` will be done **automatically** by `build.d`.

### Cleanup

If you would like to cleanup all of the compiled files generated after compilation, you can do:

```console
$ ./build.d -d
```

**Note** this will delete `libqsld.a` and any object or binary files related to the examples.
