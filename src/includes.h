#ifndef INCLUDES_H
#define INCLUDES_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <sstream>
#include <filesystem>

#include <vector>
#include <list>
#include <map>
#include <utility>
#include <set>
#include <initializer_list>

#include <string>

#include <memory>

#include <algorithm>
#include <numeric>
#include <limits>

#include <stdexcept>
#include <cassert>

#include <boost/logic/tribool.hpp>

#include <optional>
#include <wedge/wedge.h>

using std::ostream;
using std::istream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::stringstream;

using std::vector;
using std::list;
using std::map;
using std::pair;
using std::make_pair;
using std::set;
using std::initializer_list;

using std::string;
using namespace std::string_literals;

using std::unique_ptr;
using std::make_unique;

using std::numeric_limits;
using std::swap;
using std::move;

using std::optional;
using std::nullopt;

using boost::logic::tribool;
using boost::logic::indeterminate;

using namespace Wedge;
#endif
