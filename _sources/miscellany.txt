Miscellany
==========

Motivation
----------
Primer3 is a very mature and widely used primer design platform; however, 
integrating it into more complex pipelines has historically required obtuse
subprocess wrappers and inefficient file IO. We hope that this Python API will
allow for cleaner and more efficient integration of the Primer3 core 
functionalities into your projects and applications. 

Derivative work
---------------
The included ``libprimer3`` source is a dervative of the Primer3 v2.3.7 
library that includes a number of optimizations and improvements to make it 
more amenable to Python C API / Cython bindings. Every effort has been made
to insure that the underlying algorithms and calculations deterministically 
product the same results as the vanilla 2.3.7 Primer3 library.

Future work
-----------
We are happy with the current state of **Primer3-py** and have no intention of
making any major changes to the API or underlying code base in the near future.
We will provide support for bug fixes and compatibility improvements, and will
consider incorporating any changes included in major releases of the main 
Primer3 project.

Contributions
-------------
We welcome any contributions that improve the stability, compatibility or 
test coverage (shoutout to pramasoul (https://github.com/pramasoul) for his 
help with the testing framework). The best way to interface with the 
project and dev team is through Github 
(https://github.com/libnano/primer3-py).

Licensing
---------
Citations should reference the `lastest Primer3 paper 
<http://nar.oxfordjournals.org/content/early/2012/06/21/nar.gks596>`_::

  Untergasser, Andreas, et al. "Primer3—new capabilities and interfaces." 
  Nucleic acids research 40.15 (2012): e115-e115.
  doi: 10.1093/nar/gks596

All project code, including the derivative Primer3 library, is licensed 
under GPLv2. The included Python and Python C API bindings are 
Copyright (c) 2014 Ben Pruitt, Nick Conway; Wyss Institute for 
Biologically Inspired Engineering.

See LICENSE for full GPLv2 license.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
