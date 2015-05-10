/*
 * blib.hh -- basic library description
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
 

This chapter describes some classes that implement basic functionality
and help.  The next chapters describe the classes that represent the
main glsim abstractions.


@{

\defgroup Error Error handling and debugging aid

\subsection{Source context and backtrace}

We first define classes to aid in reporting source context (position
in a source file) and backtrace information.  This information is
 included in the exceptions defined below, and can be printed when
catching the exception to aid debugging


\defgroup Exceptions Exceptions for glsim (including source context)

Here we define the base exceptions for glsim.  We create our own
`Logic_error` and `Runtime_error` exceptions that inherit
from the standard `logic_error` and `runtime_error` exceptions,
adding a `Source_context` argument.  These classes add the source
context description to the exception description argument (so that it
can be displayed through the `what()` method, plus another method
that gives access to the Backtrace object stored in
Source_context.

\defgroup Logging Logging

This providesa simple log stream class.  The idea is that `main()` will
initialize a global `glsim::logs` object (of type `glsim::logger`
with the desired verbosity level.  Applicaiton code then writes to the
log stream indicating the loglevel of the message as 

     logs(glsim::info) << "Information message\n";

To initialize the log stream, the main module must (see
example in chapter XXXXREF) set the verbosity level and actual
output stream calling `glsim::logs.set_stream(stream,loglevel)`.

\textbf{NOT IMPLEMENTED} A second output stream can be set with
[[set_additional_stream]] (same syntax), possibly with a different
verbosity level.  This will result in some messages being copied to
both streams (depending on the levels).


\defgroup Random Random numbers

Random numbers are not needed in all simulations, but they are so
frequently used that we provide an interface integrated into \glsim's
conventions.  \glsim\ uses the GSL's generators.  These classes are
basically a C++ wrapper around the GSL interface, which is not hard since
the GSL design is object-oriented.  We use a class to hold one of the
available generators, and another set that represent random number
distributions.  These provide the actual numbers to be used in the
simulation, but rely on (pseudo)random numbers provided by
[[random_number_generator]].

The design of the classes allows to define a random distribution
object (descended from [[distribution_base]]) without providing a
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


\defgroup Simulation Main simulation objects


@ \chapter{Parameters}

We provide here two classes, which correspond to our parameters
abstraction.  The library user will derive from here to define
parameters as needed as shown below.  The actual parsing is done
through Boost::program\_options.  One or more [[Parameters]]
descendants will be used to define the parameters to be read from a
parameter ([[.ini]]) file.  Typically, parameters will be defined at
several places, but parsing needs to be done only once all parameters
are defined (it is far easier to use Boost this way).  One can do this
by defining a single object from a type that is derived from
[[Parameters]].  It turns out however that it is much more convenient
to have different objects (descending from [[Parameters]] but of
different type), defined in places scattered all over the code, but
which share the [[.ini]] file.  Our solution is to define each
[[Parameters]] object with a \emph{scope} (designated by a string).
All [[Parameters]] objects defined within a scope will consolidate the
parameters information in one [[boost::program_options]] object, and
the parsing function will read all of them from a single file.
Different scopes should correspond to different files, if this
complexity is not needed, scoping can be ignored and a default scope
will be used.

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


@ \section{Parameters}

To read parameters, the user declares a class inherited from
[[Parameters]].  The constructor of the derived class must declare
the parameters to be read by calling [[parm_file_options]], which is
an object of type [[options_description]] from
[[boost::program_options]].  See the example in [[test]] and the Boost
documentation for the declaration syntax.  Optionally, a scope can be
given when the object is created.  It is legal to define two objects
of the same class with different scopes.

To actually read the parameters from the file, one must call
[[Parameters::parse(char*)]] passing it a file name.  The parsing
should be done only once (in the simulation, from an [[Environent]]
object, which see).  After that, [[Parameters::count]] and
[[Parameters::value]] may be called to load parameters as desired.

The backend for parameter reading is Boost::program\_options.  Though
we allow that the parameter definition be scattered all over, the
definitios are actually collected in a single static object.  All the
parameters must be defined by the time the parser is called.  As a
result, the library user \emph{must not declary any global
  [[Parameter]] object,} or the static member objects may fail to be
properly initialized.



@}

*/

/*

\include{environment}
\include{configuration}
\include{simulation}

*/


/*

@ The implementation is rather simple.  [[operator()]] simply returns
a reference to the actual stream that will do the output.  This stream
can be
\begin{enumerate}
\item an actual output stream such as [[std::cerr]] or a
  [[std::ofstream]] linked to a file,
\item a tee stream to implement writing to two sinks with just one
  insertion operator, or
\item a null stream that discards everything sent to it.
\end{enumerate}

The constructor initializes the logstream at the highest level and
outputing to [[std::cout]] (this should be changed by [[main()]]), and
sets up the null stream.  The easiest way to do the last appears to be
the trick suggested in StackOverflow
(http://stackoverflow.com/questions/6240950/platform-independent-dev-null-in-c/6240980\#6240980):
define a [[std::ostream]] initialized with a null pointer.  This
creates a stream without a [[streambuf]] buffer, so that the stream is
left in an error state and never outputs anything (I also believe that
it skips all formatting, thus calls to the insertion operator should
be rather fast).


@ To set the actual stream and loglevel we must store pointers to the
null stream or to the actual stream as appopriate in the [[logstream]]
array.  To allow simultaneous output to two streams we will have to
add an aditional ``tee'' stream (TODO).

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% How to build a tee stream
%
%
% From http://wordaligned.org/articles/cpp-streambufs (by Thomas
% Guest).
%

<<not used now>>=
#include <streambuf>

class teebuf: public std::streambuf
{
public:
    // Construct a streambuf which tees output to both input
    // streambufs.
    teebuf(std::streambuf * sb1, std::streambuf * sb2)
        : sb1(sb1)
        , sb2(sb2)
    {
    }
private:
    // This tee buffer has no buffer. So every character "overflows"
    // and can be put directly into the teed buffers.
    virtual int overflow(int c)
    {
        if (c == EOF)
        {
            return !EOF;
        }
        else
        {
            int const r1 = sb1->sputc(c);
            int const r2 = sb2->sputc(c);
            return r1 == EOF || r2 == EOF ? EOF : c;
        }
    }
    
    // Sync both teed buffers.
    virtual int sync()
    {
        int const r1 = sb1->pubsync();
        int const r2 = sb2->pubsync();
        return r1 == 0 && r2 == 0 ? 0 : -1;
    }   
private:
    std::streambuf * sb1;
    std::streambuf * sb2;
};

// Helper class

class teestream : public std::ostream
{
public:
    // Construct an ostream which tees output to the supplied
    // ostreams.
    teestream(std::ostream & o1, std::ostream & o2);
private:
    teebuf tbuf;
};

teestream::teestream(std::ostream & o1, std::ostream & o2)
  : std::ostream(&tbuf)
  , tbuf(o1.rdbuf(), o2.rdbuf())
{
}

// Usage

#include <fstream>
#include <iostream>
#include <teestream>

int main()
{
    std::ofstream log("hello-world.log");
    teestream tee(std::cout, log);
    tee << "Hello, world!\n";
    return 0;
}

// More generic version

template <typename char_type,
          typename traits = std::char_traits<char_type> >
class basic_teebuf:
    public std::basic_streambuf<char_type, traits>
{
public:
    typedef typename traits::int_type int_type;
    
    basic_teebuf(std::basic_streambuf<char_type, traits> * sb1,
                 std::basic_streambuf<char_type, traits> * sb2)
      : sb1(sb1)
      , sb2(sb2)
    {
    }
    
private:    
    virtual int sync()
    {
        int const r1 = sb1->pubsync();
        int const r2 = sb2->pubsync();
        return r1 == 0 && r2 == 0 ? 0 : -1;
    }
    
    virtual int_type overflow(int_type c)
    {
        int_type const eof = traits::eof();
        
        if (traits::eq_int_type(c, eof))
        {
            return traits::not_eof(c);
        }
        else
        {
            char_type const ch = traits::to_char_type(c);
            int_type const r1 = sb1->sputc(ch);
            int_type const r2 = sb2->sputc(ch);
            
            return
                traits::eq_int_type(r1, eof) ||
                traits::eq_int_type(r2, eof) ? eof : c;
        }
    }
    
private:
    std::basic_streambuf<char_type, traits> * sb1;
    std::basic_streambuf<char_type, traits> * sb2;
};

typedef basic_teebuf<char> teebuf;

*/
