# boxcode2d

This repository contains an order-independent version of
the original Ethridge-Greengard volume FMM for the Poisson
equation.

## Compiling

Running

```
make test
```

in the top-level directory should compile the static library and
run a test.

## Usage

There is not much documentation for now. See test/test_pbox2drouts.f
for a test that runs the FMM.

## Generating additional tables (expert users)

To generate a new table, go to the gen directory. Run

```
make -f gentabs_pbox2d.make
```

And enter the order and table type you'd like to make.
This will generate a file with a name like pbox2dtab*.f in
the gen/ folder. Move this file to the tab/ folder. Then,
edit src/pbox2dreftabs.f to add this option to the list.
Note that high order tables take a long time to generate.
You can select new compilation options by editing a make.inc
file in that directory.

## License

The boxcode2d library is available under an Apache 2.0 license, see
LICENSE for details. This applies to the contents of the src folder
and the compiled libraries. 

The adaptive integration code used to generate tables 
and run unit tests (gen/dcuhre.f) is subject to the
ACM digital library licensing (see gen/ACM_LICENSE), which
is more restrictive than the Apache 2.0 License. That is why
dcuhre.f is not included in the library. Files in the
test folder may have their own copyright information. 

## TODO

- Fast generation of tables?
- Automate incorporating the tables better (generate pbox2dreftabs.f?)
- Helmholtz/modified Helmholtz codes tables/FMMs
- Hessian tables for Poisson