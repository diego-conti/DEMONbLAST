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
#ifndef TASKRUNNER_H
#define TASKRUNNER_H

#include <future>
using std::future;

class TaskWithFileOutput {
public:
  TaskWithFileOutput(string filename) : filename_{filename} {}
  virtual future<void> run_and_write_to_file()=0;
protected:
  ofstream open_file() const {
    return ofstream{filename_,std::ofstream::out | std::ofstream::trunc};
   }
private:
  string filename_;
};

template<typename Closure, typename Args>
class TaskWithFileOutputImpl : public TaskWithFileOutput {
  Closure closure_;
  Args args_;
public:
  TaskWithFileOutputImpl(Closure&& closure, string filename, Args&& args) : TaskWithFileOutput{filename}, closure_{std::forward<Closure>(closure)}, args_{std::forward<Args>(args)} {}
 
  future<void> run_and_write_to_file() override {
  		using namespace std;
      return async(launch::async,
        [this] () {
            stringstream output;
            closure_(args_,output);
            if (!output.str().empty()) open_file()<<output.str();
        });
  }
};



template<typename Closure, typename Args>
unique_ptr<TaskWithFileOutput> make_task_with_file_output(Closure&& closure,string filename, Args&& args) {
  return unique_ptr<TaskWithFileOutput>{new TaskWithFileOutputImpl<Closure,Args>{std::forward<Closure>(closure),filename,std::forward<Args>(args)}};
}

class TaskRunner {
  vector<unique_ptr<TaskWithFileOutput>> tasks;
public:
  template<typename Closure, typename FunctionMappingArgsToFileName, typename ListOfArgs>
  TaskRunner(Closure&& closure,  FunctionMappingArgsToFileName&& function_mapping_args_to_filename, const ListOfArgs& args) {
      for (auto some_args: args) {
        auto filename = function_mapping_args_to_filename(some_args);
        tasks.push_back(make_task_with_file_output(std::forward<Closure>(closure),filename, std::move(some_args)));
      }
  }
  void run_and_write_to_file() {
    vector<future<void>> handles;
    for (auto& task: tasks)
      handles.push_back(task->run_and_write_to_file());
    for (auto& handle : handles) {
       handle.get();      
    }
  }
};

#endif
