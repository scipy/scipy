/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "HighsModelBuilder.h"

#include <math.h>

#include "lp_data/HConst.h"

HighsModelBuilder::~HighsModelBuilder() {
  while (this->variables.size() > 0) {
    HighsVar* variable;
    variable = this->variables.front();
    this->variables.pop_front();

    // find all coefficients corresponding to the variable
    VarConsCoefsMap::iterator it =
        this->variableConstraintCoefficientMap.find(variable);
    if (it != this->variableConstraintCoefficientMap.end()) {
      std::list<HighsLinearConsCoef*>* coefficients = it->second;

      while (coefficients->size() > 0) {
        HighsLinearConsCoef* coef = coefficients->front();
        coefficients->pop_front();
        // remove coefficient from constraint

        CoefConsMap::iterator iter = this->coefficientConstraintMap.find(coef);
        assert(iter != this->coefficientConstraintMap.end());
        HighsLinearCons* constraint = iter->second;
        VarConsCoefMap::iterator iterator =
            constraint->linearCoefs.find(variable);
        assert(iterator != constraint->linearCoefs.end());
        constraint->linearCoefs.erase(iterator);
        this->coefficientConstraintMap.erase(iter);
        delete coef;
      }
      VarConsMap::iterator iter = this->variableConstraintMap.find(variable);
      if (iter != variableConstraintMap.end()) {
        std::list<HighsLinearCons*>* conslist = iter->second;
        conslist->clear();
        this->variableConstraintMap.erase(iter);
        delete conslist;
      }

      this->variableConstraintCoefficientMap.erase(it);
      delete coefficients;
    }

    delete variable;
  }

  while (this->linearConstraints.size() > 0) {
    HighsLinearCons* constraint;
    constraint = this->linearConstraints.front();
    this->linearConstraints.pop_front();

    delete constraint;
  }
}

HighsVar::HighsVar(const char* name, double lo, double hi, double obj,
                   HighsVarType type) {
  // create a copy of the name
  if (name != NULL) {
    int namelen = strlen(name);
    this->name = new char[namelen + 1];
    strcpy(this->name, name);
  } else {
    this->name = NULL;
  }

  // copy all remaining data
  this->lowerBound = fmax(-HIGHS_CONST_INF, lo);
  this->upperBound = fmin(HIGHS_CONST_INF, hi);
  this->obj = obj;
  this->type = type;
}

HighsVar::~HighsVar() {
  if (this->name != NULL) {
    delete[] this->name;
  }
}

HighsCons::HighsCons(const char* name, double lo, double hi) {
  // create a copy of the name
  if (name != NULL) {
    int namelen = strlen(name);
    this->name = new char[namelen + 1];
    strcpy(this->name, name);
  } else {
    this->name = NULL;
  }

  // copy all remaining data
  this->lowerBound = lo;
  this->upperBound = hi;
}

HighsCons::~HighsCons() {
  if (this->name != NULL) {
    delete[] this->name;
  }
}

HighsLinearCons::HighsLinearCons(const char* name, double lo, double hi)
    : HighsCons(name, lo, hi) {}

HighsLinearCons::~HighsLinearCons() {}

HighsLinearConsCoef::HighsLinearConsCoef(HighsVar* var, double coef) {
  this->var = var;
  this->coef = coef;
}

HighsLinearConsCoef::~HighsLinearConsCoef() {}

void HighsModelBuilder::HighsCreateVar(const char* name, double lo, double hi,
                                       double obj, HighsVarType type,
                                       HighsVar** var) {
  if (name != NULL) {
    // make sure name is available
    VarMap::iterator it = this->variableMap.find(name);
    if (it != this->variableMap.end()) {
      // name already in use
      // TODO: Error Message
      return;
    }
  }

  // create the new variable and add it to the model
  *var = new HighsVar(name, lo, hi, obj, type);
  this->variables.push_back(*var);
  if (name != NULL) {
    this->variableMap.insert(VarMap::value_type((*var)->name, *var));
  }
}

void HighsModelBuilder::HighsCreateVar(const char* name, HighsVar** var) {
  this->HighsCreateVar(name, 0.0, HIGHS_CONST_INF, 0.0, HighsVarType::CONT,
                       var);
}

void HighsModelBuilder::HighsGetOrCreateVarByName(const char* name,
                                                  HighsVar** var) {
  this->HighsGetVarByName(name, var);
  if (*var == NULL) {
    this->HighsCreateVar(name, var);
  }
}

void HighsModelBuilder::HighsCreateVar(HighsVar** var) {
  this->HighsCreateVar(NULL, var);
}

void HighsModelBuilder::HighsGetVarByName(const char* name, HighsVar** var) {
  VarMap::iterator it = this->variableMap.find(name);
  if (it != this->variableMap.end()) {
    *var = it->second;
  } else {
    // variable not found
    // TODO: Error Message
    *var = NULL;
  }
}

void HighsModelBuilder::HighsRemoveVar(HighsVar* var) {
  // check that variable is no longer used in any constraints
  // TODO

  // remove variable from map
  VarMap::iterator it = this->variableMap.find(var->name);
  if (it == this->variableMap.end()) {
    // variable no longer in Model?
    // TODO: Error Message
    return;
  }
  this->variableMap.erase(var->name);

  // remove variable from list
  // TODO
  return;
}

void HighsModelBuilder::HighsCreateLinearCons(const char* name, double lo,
                                              double hi,
                                              HighsLinearCons** cons) {
  if (name != NULL) {
    // make sure name is available
    ConsMap::iterator it = this->constraintMap.find(name);
    if (it != this->constraintMap.end()) {
      // name already in use
      // TODO: Error Message
      return;
    }
  }

  // create the new constraint and add it to the model
  *cons = new HighsLinearCons(name, lo, hi);
  this->linearConstraints.push_back(*cons);
  if (name != NULL) {
    this->constraintMap.insert(ConsMap::value_type((*cons)->name, *cons));
  }
}

void HighsModelBuilder::HighsCreateLinearCons(const char* name,
                                              HighsLinearCons** cons) {
  this->HighsCreateLinearCons(name, -HIGHS_CONST_INF, HIGHS_CONST_INF, cons);
}

void HighsModelBuilder::HighsCreateLinearCons(HighsLinearCons** cons) {
  this->HighsCreateLinearCons(NULL, cons);
}

void HighsModelBuilder::HighsGetLinearConsByName() {
  // TODO
}

void HighsModelBuilder::HighsDestroyLinearCons() {
  // TODO
}

void HighsModelBuilder::HighsCreateLinearConsCoef(
    HighsVar* var, double coef, HighsLinearConsCoef** consCoef) {
  *consCoef = new HighsLinearConsCoef(var, coef);
  VarConsCoefsMap::iterator it =
      this->variableConstraintCoefficientMap.find(var);
  if (it != this->variableConstraintCoefficientMap.end()) {
    it->second->push_back(*consCoef);
    ;
  } else {
    std::list<HighsLinearConsCoef*>* coefList =
        new std::list<HighsLinearConsCoef*>;
    coefList->push_back(*consCoef);
    this->variableConstraintCoefficientMap.insert(
        VarConsCoefsMap::value_type(var, coefList));
  }
}

int HighsModelBuilder::getNumberOfVariables() { return this->variables.size(); }

void HighsModelBuilder::HighsAddLinearConsCoefToCons(
    HighsLinearCons* cons, HighsLinearConsCoef* coef) {
  VarConsCoefMap::iterator it = cons->linearCoefs.find(coef->var);
  if (it != cons->linearCoefs.end()) {
    // constraint already has a coefficient for this variable
  } else {
    coefficientConstraintMap.insert(CoefConsMap::value_type(coef, cons));
    cons->linearCoefs.insert(VarConsCoefMap::value_type(coef->var, coef));
    VarConsMap::iterator it = this->variableConstraintMap.find(coef->var);
    if (it != this->variableConstraintMap.end()) {
      it->second->push_back(cons);
    } else {
      std::list<HighsLinearCons*>* consList = new std::list<HighsLinearCons*>;
      consList->push_back(cons);
      this->variableConstraintMap.insert(
          VarConsMap::value_type(coef->var, consList));
    }
  }
}

void HighsModelBuilder::HighsBuildTechnicalModel(HighsLp* lp) {
  lp->numCol_ = this->variables.size();
  lp->numRow_ = this->linearConstraints.size();

  lp->sense_ = this->objSense;

  // determine order of variables
  HighsVar** variables = new HighsVar*[lp->numCol_];
  for (int i = 0; i < lp->numCol_; i++) {
    HighsVar* front = this->variables.front();
    this->variables.pop_front();
    this->variables.push_back(front);
    variables[i] = front;
    lp->colCost_.push_back(front->obj);
    lp->colLower_.push_back(front->lowerBound);
    lp->colUpper_.push_back(front->upperBound);
  }

  // determine order of constraints
  HighsLinearCons** constraints = new HighsLinearCons*[lp->numRow_];
  for (int i = 0; i < lp->numRow_; i++) {
    HighsLinearCons* front = this->linearConstraints.front();
    this->linearConstraints.pop_front();
    this->linearConstraints.push_back(front);
    constraints[i] = front;
    lp->rowLower_.push_back(front->lowerBound);
    lp->rowUpper_.push_back(front->upperBound);
  }

  // handle constraints
  lp->Astart_.clear();
  lp->Astart_.push_back(0);
  for (int var = 0; var < lp->numCol_; var++) {
    VarConsCoefsMap::iterator iter =
        this->variableConstraintCoefficientMap.find(variables[var]);
    if (iter != this->variableConstraintCoefficientMap.end()) {
      std::list<HighsLinearConsCoef*>* coefs = iter->second;
      int numberOfCoefficients = coefs->size();

      lp->Astart_.push_back(lp->Astart_[var] + numberOfCoefficients);

      for (int coef = 0; coef < numberOfCoefficients; coef++) {
        HighsLinearConsCoef* front = coefs->front();
        coefs->pop_front();
        coefs->push_back(front);
        lp->Avalue_.push_back(front->coef);
        CoefConsMap::iterator it = this->coefficientConstraintMap.find(front);
        if (it != this->coefficientConstraintMap.end()) {
          // find index of constraint
          HighsCons* currentCons = it->second;
          for (int cons = 0; cons < lp->numRow_; cons++) {
            if (constraints[cons] == currentCons) {
              lp->Aindex_.push_back(cons);
              break;
            }
          }
        } else {
          // ERROR
        }
      }
    }
  }

  delete[] variables;
  delete[] constraints;
}
