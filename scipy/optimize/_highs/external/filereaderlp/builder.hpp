#ifndef __READERLP_BUILDER_HPP__
#define __READERLP_BUILDER_HPP__

#include <map>
#include <memory>
#include <string>

#include "model.hpp"

struct Builder { 
   std::map<std::string, std::shared_ptr<Variable>> variables; 

   Model model;

   std::shared_ptr<Variable> getvarbyname(std::string name) {
      if (variables.count(name) == 0) {
         variables[name] = std::shared_ptr<Variable>(new Variable(name));
         model.variables.push_back(variables[name]);
      }
      return variables[name];
   }
};

#endif
