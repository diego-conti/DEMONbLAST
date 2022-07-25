/*  Copyright (C) 2018 by Diego Conti, diego.conti@unimib.it      
                                                                     
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
#ifndef OPTIONS_H
#define OPTIONS_H

template<typename Option>
struct Options {
  using underlying=std::underlying_type_t<Option>;
  Option options=Option::dflt;
  void log() const {nice_log<<"options = "<<static_cast<underlying>(options)<<endl;}
  void set(Option option) {as_underlying_type()|=static_cast<underlying>(option);}
  void clear(Option option) {set(option); as_underlying_type()^=static_cast<underlying>(option);}
  bool has(Option option) const {return static_cast<underlying>(options)&static_cast<underlying>(option);}
  bool operator==(Options other_options) const {return static_cast<underlying>(options)==static_cast<underlying>(other_options.options);}
  bool operator==(Option other_options) const {return static_cast<underlying>(options)==static_cast<underlying>(other_options);}
  template<typename Opt>
  bool operator!=(Opt other_options) const {return !(*this==other_options);}
private:
  underlying& as_underlying_type() {
  	return *reinterpret_cast<underlying*>(&options);
  }
};

#endif
