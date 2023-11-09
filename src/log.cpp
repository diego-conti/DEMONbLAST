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
#include "log.h"
#include <thread>
#include <sstream>
thread_local Log nice_log;

#ifndef NDEBUG
std::string current_thread_id() {
    std::stringstream s;
    s<<std::this_thread::get_id();
    return s.str();    
}

Log::Log() : os{"log/"+current_thread_id()} {}

void Log::print(Counter& counter) {
  if (++counter.count%counter.step==0)
    os<<counter.description<<" "<<++counter.count<<std::endl;
  else os<<'.';
}

#endif
