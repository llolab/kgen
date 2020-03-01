# kgen

The *kgen* software provides functionality to efficiently compute inverses of genetic similarity matrices for partially observed multi-phenotype data. The source, manual and build instructions for kgen are provided in the Supplementary Material for the paper *Kinship solutions for partially observed multi-phenotype data* (published in the *Journal of Computational Biology*), and in this repository.

## Files

* **`supmat.pdf`** - The Supplementary Material for the paper *Kinship solutions for partially observed multi-phenotype data*. This material includes additional information about the experiments, and also the *kgen* manual.
* **`kgen`** - A static *ELF* binary for the *kgen* software built with `gcc 9.1.0` and a `Xeon E5-2683` on `CentOS 7.6.1810`.
* **`Makefile`** - Build instructions for the *kgen* software.
* **`kgen.c`** and **`kgen.h`** - Source code for the *kgen* software. This source is released under the BSD 2-clause license (provided below), and this source makes use of the *Intel MKL* or *netlib* packages.
* **`libqrupdate.a`** - A static library for the libqrupdate package (by Jaroslav Hajek). This package is released under the GPLv2 license, and source code is available at [sourceforge](https://sourceforge.net/projects/qrupdate/).

## Changelog

* February 29th 2020 - *kgen v1*. Initial release. The flag `-t` is currently unimplemented.

## License

```
kgen v1. Copyright (c) 2020. Lloyd T. Elliott.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the
following conditions are met:

1. Redistributions of source code must retain the above
copyright notice, this  list of conditions and the
following disclaimer.

2. Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials
provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE   OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```
