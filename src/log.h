/*  Copyright (C) 2018-2023 by Diego Conti, diego.conti@unipi.it      
                                                                     
    This file is part of DEMONbLAST
	                                                                     
    DEMONbLAST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DEMONbLAST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DEMONbLAST.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef LOG_H
#define LOG_H

#ifdef NDEBUG
#include <ostream>
class Log {
public:
  template<typename StringType> void* counter(StringType&& counter_description, int step=1) const {return nullptr;}
};

template<typename T>
Log& operator<<(Log& log, T&&) {
  return log;
}
inline Log& operator<<(Log& log, std::ostream&(*)(std::ostream&)) {return log;}

#else
#include <fstream>
#include <string>


class Counter {
  friend class Log;
  std::string description;
  int count=0;
  int step;
  template<typename StringType> Counter(StringType&& counter_description, int print_step) : description{std::forward<StringType>(counter_description)},step{print_step} {}
};

class Log {
  std::ofstream os;
public:
  Log();
  template<typename T> void print(T&& t) {os<<std::forward<T>(t);}
  void print(Counter& counter);
  template<typename StringType> Counter counter(StringType&& counter_description, int step=1) const {return {std::forward<StringType>(counter_description),step};}
};
template<typename T>
Log& operator<<(Log& log, T&& t) {
  log.print(std::forward<T>(t));
  return log;
}
inline Log& operator<<(Log& log, std::ostream&(*f)(std::ostream&)) {
  log.print(f);
  return log;
}
#endif

extern thread_local Log nice_log;

#endif
