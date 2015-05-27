/*
 * log.cc -- a simple logging stream
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

#include "exception.hh"
#define _LOGGING_CC_
#include "log.hh"

namespace glsim{

/**
   The main program should call this function to choose the stream
   that the log messages will be sent to, and at which verbosity
   level.  Message with level lower than the level set here are
   ignored.
 */
logger& logger::set_stream(std::ostream& str,loglevel level)
{
  // To set the actual stream and loglevel we must store pointers to the
  // null stream or to the actual stream as appopriate in the logstream
  // array.  To allow simultaneous output to two streams we will have to
  // add an aditional ``tee'' stream (TODO).  stream1.str=&str;

  stream1.level=level;
  stream1.str=&str;
  for (int i=0; i<level; i++)
    logstream[i]=&nullstream;
  for (int i=level; i<=error; i++)
    logstream[i]=stream1.str;

  return *this;
}

logger& logger::set_additional_stream(std::ostream&,loglevel)
{
  throw Unimplemented("set_additional_stream",HERE);

  return *this;
}

} /* namespace */

/*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% How to build a tee stream
%
%
% From http://wordaligned.org/articles/cpp-streambufs (by Thomas
% Guest).
%

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
