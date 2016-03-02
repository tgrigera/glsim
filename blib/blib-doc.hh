/*
 * blib-doc.hh -- Documentation for the basic library
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015
 * by Tomas S. Grigera.
 * 
 * glsim is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, with the additional requirements of
 * attribution and nonmisrepresentation. You can use either version 3, or
 * (at your option) any later version.
 * 
 * Additional terms under GNU GPL version 3 section 7:
 * 
 * When you redistribute this software, you are required to preserve its
 * author attributions. If you distribute verbatim copies, you must not
 * alter the AUTHORS file or attributions inserted in the source files,
 * and you must not change the software's name. If you distribute a
 * modified copy, then you must give clear notice that your work is
 * different from but based on glsim. You must distribute it under a
 * different name, but include a prominent notice specifying that "(your
 * package) is based on glsim version x.x", and provide a pointer to the
 * glsim distribution.
 *
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */


/**\defgroup Blib Basic library
 
The main library is `glsim` proper.  The utility programs and the
lattice and off-lattice libraries depend on the modules included here.
It includes class hierarchies to represent the main simulation
abstractions (Parameters, Configuration, Environment and Simulation),
as well as utilities for logging and error handling and interfaces for
basic object such as random number generators and self-describing disk
files.

### Error handling and logging

This includes a simple log stream, the exceptions thrown by the
library and a class to allow printing of source context and stack
backtrace when catching the exceptions.  Documented \link Error
here\endlink.

### Interfaces to other basic libraries

Convenient objetc-orinted interfaces are provided for \link Random
random numbers\endlink (from the GSL library) and for \link HDF HDF5
files\endlink.

### Main simulation abstractions

This includes de base of the four main hierarchies of the library,
representing \link Parameter parameters\endlink, the \link
Configuration configuration\endlink, the \link Environment
environent\endlink, and the \link Simulation simulation\endlink.


*****************************************************************

@{ \defgroup Error Error handling, logging and debugging aid
    \ingroup Blib
@{

\defgroup Scontext Source context and backtrace

These classes aid in reporting source context (position in a source
file) and backtrace information.  This information is included in the
\link Exceptions exceptions\endlink thrown by `glsim`, and can be
printed when catching the exception to aid debugging.

\defgroup Exceptions Exceptions thrown by the library

Here we define the base exceptions for glsim.  We create our own
`Logic_error` and `Runtime_error` exceptions that inherit
from the standard `logic_error` and `runtime_error` exceptions,
adding a `Source_context` argument.  These classes add the source
context description to the exception description argument so that it
can be displayed through the `what()` method, plus another method
that gives access to the Backtrace object stored in
Source_context.

\defgroup Logging Logging

We implement a simple log stream class.  The idea is that `main()`
will initialize a global `glsim::logs` object (of type
`glsim::logger`) with the desired verbosity level.  Code in other
functions or classes can then write to the log stream indicating the
loglevel of the message as

     logs(glsim::info) << "Information message\n";

To initialize the log stream, the main module must set the verbosity
level and actual output stream calling
`glsim::logs.set_stream(stream,loglevel)`.

*NOT IMPLEMENTED:* A second output stream can be set with
set_additional_stream (same syntax), possibly with a different
verbosity level.  This will result in some messages being copied to
both streams (depending on the levels).

@}

\defgroup Interfaces Interfaces to other basic libraries

@{

\defgroup Random Random numbers

Random numbers are not needed in all simulations, but they are so
frequently used that we provide an interface integrated into `glsim`'s
conventions.  `glsim` uses the Gnu Scientific Library's generators.
We provide a C++ wrapper around the GSL interface, which is not hard
since the GSL design is object-oriented.  We use a class to hold one
of the available generators, and another set that represent random
number distributions.  These provide the actual numbers to be used in
the simulation, but rely on (pseudo)random numbers provided by
Random_number_generator.

The design of the classes allows to define a random distribution
object (descended from Distribution_base) without providing a
generator.  This is allows the programmer to write a class or function
that uses random numbers without initializing its own private
generator (which in general is not desired as, among other things it
would introduce the need for several seeds) or requiring to receive a
reference or pointer to a generator as an argument.  In this way the
user of the class still gets to decide which generator will be used,
while keeping the interface simpler.  This is achieved with
a scope mechanism as that of the Environment classes.  Only one
generator per scope is allowed (though it is fine to define multiple
distributions sharing scope and thus generator).

Saving and loading are supported, but loading into an object
initialized to use a GSL generator different than the one it is
initialized to is not correctly handled.  This should be considered a
bug.

The scope cannot be changed through loading; it is not written
but taken and left unchanged from creation time.  An exception is when
deserializing through pointers, where the scope cannot be provided
beforehand.  This maybe counterintuitive and can perhaps be considered
a bug; but at present I see no cheap solution.



\defgroup HDF In/out through HDF5 files

A set of classes that represents HDF5 library objects (files, groups,
and datatypes) so that they may be more easily created and
manipulated.  There is also a class (HDF_record_file) that makes it
easy to create and read a file with a simple record structure (a
header plus a table-like structure).

\defgroup FILE File manipulation

\defgroup FFT Interface to 3rd party Fast Fourier Transforms

@}


\defgroup Abstractions Main simulation abstractions

@{

Here are the base objects of our main simulation abstractions:

 - Parameters
 - Configuration
 - Environment
 - Simulation


\defgroup Parameters Parameters

Parameters are represented by one or more objects of a type derived
from \link glsim::Parameters\endlink.  The base class is able to parse
a file with a `parameter=value` structure (a `.ini` file) through the
`program_options` component of the `Boost` library.  The library user
will derive from glsim::Parameters, defining the actual parameters in
the derived constructor by calling
Parameters::parameter_file_options().  Once all parameters are
defined, the parameter file is parsed by calling
Parameters::parse().

Since parameters are tied to particular simulation components (e.g. an
observable might use a parameter to set the observation frequency), to
preserve modularity it is convenient to scatter the definition of the
parameters across several source locations.  However, parsing needs to
be done only once all parameters are defined (it is far easier to use
Boost this way).  One can do this by defining a single object from a
type that is derived from Parameters, but it turns out that it is much
more convenient to have different objects (descending from Parameters
but of different type), defined in places scattered all over the code,
but which share the `.ini` file.  Our solution is to specify a \link
Scopes scope\endlink name (a string) when defining each Parameters
object.  All Parameters objects defined within a scope will
consolidate the parameter information in one `boost::program_options`
object, and the parsing function will read all of them from a single
file.  Different scopes should correspond to different files.  If this
complexity is not needed, scoping can be ignored and a default scope
will be used.


### The command line

To parse the command line, one can create an object derived from
[[ParametersCL]].  Only one such object should exist (but it can share
scope with other [[Parameters]]).  In this case two parser functions
must be called, one for the command line and one for the file (though
the second call can be automated, see [[ParametersCL]] below).  But
note that [[ParametersCL]] descends from [[Parameters]], thus it too
consolidates all parameter definitions. As result, if a
[[ParametersCL]] object exists, the parameters defined in other
[[Parameters]] objects in the same scope can also be specified as long
options in the command line, overriding the file values.




\defgroup Environment Environment

The concept of environment is one of the major abstractions of the
library, and this module contains the classes that represent it.

By our definition, there should be an Environment class to match any
Simulation class.  The Environment class corresponding to a given
Simulation must hold _all_ the data needed to perform the simulation.
The Parameters objects are meant to be used only to initialize the
environment, and not for storing data.  The environment is endowed
with methods to allow saving to disk and restoring, while Parameters
classes have no such capability.  This is important, as the
environment can evolve during the simulation (see the method step()).

The parts of the environment essential for the simulation should be
represented by an object of a type derived from SimEnvironment below.
However, observables (which see), which will require a part of the
environment to store their internal data, are designed to be used as
``plug-ins'', and thus it is inconvenient to require that all the
environment be in a single object.  Rather we allow to define several
objects, belonging to several hierarchies (descending from
Environment), which will consolidate with each other through a
mechanism similar to the singleton construction, so that saving and
restoring can be done through a single object.  To allow for several
independent environments to coexist (as is useful for instance in
parallel simulations), environments will be consolidated only if they
belong to the same \link Scope scope\endlink labeled with a string and
defined at the time the object is created.

The intended use of the hierarchy is as folows: for every simulation
to be created, there will be an object of a type derived from
SimEnvironment.  Each of these would be in a different scope.
Plugin-like object such as observables will create private
environments (derived from Environment), in one of the scopes used
when creating the `SimEnvironment`s.  The simulation then will deal
with the `SimEnvironment`s only, and will call methods to initialize,
load, or save automatically all the consolidated environments in the
scope (in effect working conceptually as a single environment).

See the class documentation for how to correctly derive from
Environment and how the saving and restoring of environments can be
controlled.


\defgroup Simulation Simulation
\brief The base simulation classes

Currently there is only one class here (Simulation), which must be the
base for all simulations.  The user will inherit from Simulation,
defining the pure virtual step(), which carries out the actual
simulation step.  Then, after proper initialization of configuration
and environment (descended from Configuration and SimEnvironment), the
simulation object is created and the simulation is started by calling
Simulation::run() (note that Simulation wile _not_ attemtp to
initialize configuration or environment, see the prepare() function
for a way to perform this inizialization more easily).

Define step-based and target-based simulations; explain how both can
be dealt with here.  env.run_completed must be set by child sim or
specific environment.

Define run and stage.

IMPORTANT: step() must decide when the sim is finished (set
run_completed), and compute time and tolerance when applicable.


\defgroup Observable Observable
	   \brief Observables

These are notes for an observable.

For things like magnetization etc which are normally computed from
within the simulation step, an observable hierarchy is of questionable
value.  Rather we shall provide (but not now), files that will be
monitored for consistent checkpointing.

The objects sketched here are useful for the ``plug-in'' situation,
where the observation is completely independent from the simulation.
It might even be argued that this should be independent from the sim
and observe the configuration alone. We shall see.


For checkpointing I need

 - a FILE-castable object that remembers last position
 - a way to distinguish whether is second+ stage or not

 - if second stage, move to end and don't write header
 - if first stage, open file (delete if existing), write header.
 - need and environment for that.
 - need an aditional ``INIT'' method for observables? ---> es mejor
 que no, y que se arregle step.  Para optimizar, en el futuro,
 existira el step_is_noop() (en BaseEnv), y un first_step (para que la
 primera vez llame a una función y luego llame a otra).



@}


@}

*/
