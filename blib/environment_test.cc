/*
 * environment_test.cc -- test loading and saving of environments
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

#include "environment.hh"

using namespace glsim;

class floating_pars : public Parameters {
public:
  floating_pars(const char* scope=Parameters::default_scope);
} ;

floating_pars::floating_pars(const char* scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("float.alpha",po::value<double>())
    ("float.beta",po::value<double>())
    ;
}

class floating_environment : public Environment {
public:
  floating_environment(const char* scope=Parameters::default_scope);
  void step() {}

  double alpha,beta;

protected:
  void init_local(),warm_init_local();
  void step_local() {}

private:
  floating_pars par;

  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  friend class BaseEnvironment;
  virtual void vserial(iarchive_t &ar) {ar >> *this;}
  virtual void vserial(oarchive_t &ar) {ar << *this;}
} ;

inline floating_environment::floating_environment(const char* scope) :
  Environment(scope),
  par(scope),
  alpha(0), beta(0)
{}

template <typename Archive>
void floating_environment::serialize(Archive &ar,const unsigned int version)
{
  ar & boost::serialization::base_object<Environment>(*this);
  ar & alpha;
  ar & beta;
}

void floating_environment::init_local()
{
  Environment::init_local();
  floating_environment::warm_init_local();
}

void floating_environment::warm_init_local()
{
  Environment::warm_init_local();
  alpha=par.value("float.alpha").as<double>();
  beta=par.value("float.beta").as<double>();
}

void print_environment(glsim::SimEnvironment &env,floating_environment &fl)
{
  std::cout
    << "   config file ini/fin = " << env.configuration_file_ini << '/'
    << env.configuration_file_fin << '\n'
    << "   env file ini/fin    = " << env.environment_file_ini << '/'
    << env.environment_file_fin << '\n'
    << "   title               = " << env.title << '\n'
    << "   log_interval        = " << env.log_interval << '\n'
    << "   run_completed       = " << env.run_completed << '\n'
    << "   [float]\n"
    << "   alpha               = " << fl.alpha << '\n'
    << "   beta                = " << fl.beta << '\n'
    ;
}


int main(int argc, char *argv[])
{
  int rcode=1;

  SimulationCL CL("environmment_test","(C) TSG","Myscope");
  SimEnvironment main_env("Myscope");
  floating_environment fl("Myscope");

  try {

    CL.parse_command_line(argc,argv);
    main_env.init();
    fl.alpha*=1.2;
    main_env.run_completed=true;
    print_environment(main_env,fl);
    main_env.save();

    main_env.log_interval=10000;
    main_env.run_completed=false;
    fl.alpha=0;
    fl.beta=0;
    std::cout << "\nSaved and modified:\n";
    print_environment(main_env,fl);

    main_env.environment_file_ini=main_env.environment_file_fin;
    main_env.load();
    std::cout << "\nReloaded:\n";
    print_environment(main_env,fl);

    std::cout << "\nNew envs:\n";
    SimEnvironment newSE;
    newSE.environment_file_ini=main_env.environment_file_fin;
    floating_environment newfl;
    newSE.load();
    print_environment(newSE,newfl);

    rcode=0;

  } catch (const glsim::Early_stop &e) {
  } catch (const glsim::Usage_error &e) {
    CL.show_usage();
  } catch (const glsim::Runtime_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
  } catch (const glsim::Logic_error& e) {
    std::cerr << e.what() << '\n';
    std::cerr << e.backtrace();
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << '\n';
  }
  return rcode;
}
