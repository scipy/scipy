/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderLp.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "io/FilereaderLp.h"

#include <cstdarg>

#include "lp_data/HConst.h"
#include "util/stringutil.h"

FilereaderLp::FilereaderLp() {
  this->isFileBufferFullyRead = true;
  this->readingPosition = NULL;
  this->file = NULL;
  this->linelength = 0;
  this->status = LP_FILEREADER_STATUS::SUCCESS;
}

void emptyTokenQueue(std::list<LpToken*>& list) {
  while (list.size() > 0) {
    LpToken* token = list.front();
    list.pop_front();
    delete token;
  }
}

FilereaderLp::~FilereaderLp() {
  emptyTokenQueue(this->tokenQueue);
  emptyTokenQueue(this->objectiveSection);
  emptyTokenQueue(this->constraintSection);
  emptyTokenQueue(this->boundsSection);
  emptyTokenQueue(this->binSection);
  emptyTokenQueue(this->generalSection);
  emptyTokenQueue(this->semiSection);
  emptyTokenQueue(this->sosSection);
}

FilereaderRetcode FilereaderLp::readModelFromFile(const HighsOptions& options,
                                                  HighsLp& model) {
  HighsModelBuilder m;
  const char* filename = options.model_file.c_str();
  this->readModelFromFile(filename, m);
  m.HighsBuildTechnicalModel(&model);
  return FilereaderRetcode::OK;
}

FilereaderRetcode FilereaderLp::readModelFromFile(const char* filename,
                                                  HighsModelBuilder& model) {
  this->file = fopen(filename, "r");
  if (file == NULL) {
    return FilereaderRetcode::FILENOTFOUND;
  }

  this->tokenizeInput();
  if (this->status != LP_FILEREADER_STATUS::ERROR) this->splitTokens();
  if (this->status != LP_FILEREADER_STATUS::ERROR)
    this->handleObjectiveSection(model);
  if (this->status != LP_FILEREADER_STATUS::ERROR)
    this->handleConstraintSection(model);
  if (this->status != LP_FILEREADER_STATUS::ERROR)
    this->handleBoundsSection(model);
  if (this->status != LP_FILEREADER_STATUS::ERROR)
    this->handleBinarySection(model);
  if (this->status != LP_FILEREADER_STATUS::ERROR)
    this->handleGeneralSection(model);
  if (this->status != LP_FILEREADER_STATUS::ERROR)
    this->handleSemiSection(model);
  if (this->status != LP_FILEREADER_STATUS::ERROR) {
    FilereaderRetcode filereader_return_code = this->handleSosSection(model);
    if (filereader_return_code != FilereaderRetcode::OK)
      return FilereaderRetcode::PARSERERROR;
  }
  assert(this->tokenQueue.size() == 0);

  fclose(file);
  if (this->status != LP_FILEREADER_STATUS::ERROR) {
    return FilereaderRetcode::OK;
  } else {
    return FilereaderRetcode::PARSERERROR;
  }
}

void FilereaderLp::handleBinarySection(HighsModelBuilder& model) {
  if (this->binSection.size() == 0) {
    return;
  }

  LpToken* token;
  token = this->binSection.front();
  this->binSection.pop_front();
  assert(token->type == LpTokenType::SECTIONKEYWORD);
  assert(((LpTokenSectionKeyword*)token)->section == LpSectionKeyword::BIN);
  delete token;

  while (this->binSection.size() > 0) {
    LpToken* token = this->binSection.front();
    assert(token->type == LpTokenType::VARIDENTIFIER);

    HighsVar* variable;
    model.HighsGetOrCreateVarByName(((LpTokenVarIdentifier*)token)->value,
                                    &variable);
    // update bounds if necessary
    if (variable->lowerBound == 0.0 &&
        variable->upperBound == HIGHS_CONST_INF) {
      variable->upperBound = 1.0;
    }
    variable->type = HighsVarType::BIN;

    this->binSection.pop_front();
    delete token;
  }

  assert(this->binSection.size() == 0);
}

void FilereaderLp::handleGeneralSection(HighsModelBuilder& model) {
  if (this->generalSection.size() == 0) {
    return;
  }

  LpToken* token;
  token = this->generalSection.front();
  this->generalSection.pop_front();
  assert(token->type == LpTokenType::SECTIONKEYWORD);
  assert(((LpTokenSectionKeyword*)token)->section == LpSectionKeyword::GEN);
  delete token;

  while (this->generalSection.size() > 0) {
    token = this->generalSection.front();
    assert(token->type == LpTokenType::VARIDENTIFIER);

    HighsVar* variable;
    model.HighsGetOrCreateVarByName(((LpTokenVarIdentifier*)token)->value,
                                    &variable);
    variable->type = HighsVarType::GEN;

    this->generalSection.pop_front();
    delete token;
  }

  assert(this->generalSection.size() == 0);
}

void FilereaderLp::handleSemiSection(HighsModelBuilder& model) {
  if (this->semiSection.size() == 0) {
    return;
  }

  LpToken* token;
  token = this->semiSection.front();
  this->semiSection.pop_front();
  assert(token->type == LpTokenType::SECTIONKEYWORD);
  assert(((LpTokenSectionKeyword*)token)->section == LpSectionKeyword::SEMI);
  delete token;

  while (this->semiSection.size() > 0) {
    token = this->semiSection.front();
    assert(token->type == LpTokenType::VARIDENTIFIER);

    HighsVar* variable;
    model.HighsGetOrCreateVarByName(((LpTokenVarIdentifier*)token)->value,
                                    &variable);
    variable->type = HighsVarType::SEMI;

    this->semiSection.pop_front();
    delete token;
  }

  assert(this->semiSection.size() == 0);
}

FilereaderRetcode FilereaderLp::handleSosSection(HighsModelBuilder& model) {
#ifdef HiGHSDEV
  printf("SoS section is not currenlty supported by the .lp filereader.\n");
#endif

  if (this->sosSection.size() == 0) {
    return FilereaderRetcode::OK;
  }

  LpToken* token;
  token = this->sosSection.front();
  this->sosSection.pop_front();
  assert(token->type == LpTokenType::SECTIONKEYWORD);
  assert(((LpTokenSectionKeyword*)token)->section == LpSectionKeyword::SOS);
  delete token;

  while (this->sosSection.size() > 0) {
    LpToken* token = this->sosSection.front();

    this->sosSection.pop_front();
    delete token;
  }
  return FilereaderRetcode::NOT_IMPLEMENTED;
}

void FilereaderLp::handleBoundsSection(HighsModelBuilder& model) {
  if (this->boundsSection.size() == 0) {
    return;
  }

  LpToken* token;
  token = this->boundsSection.front();
  this->boundsSection.pop_front();
  assert(token->type == LpTokenType::SECTIONKEYWORD);
  assert(((LpTokenSectionKeyword*)token)->section == LpSectionKeyword::BOUNDS);
  delete token;

  while (this->boundsSection.size() > 1) {
    LpToken* current = this->boundsSection.front();
    this->boundsSection.pop_front();
    LpToken* next = this->boundsSection.front();
    this->boundsSection.pop_front();
    // LpToken* nextnext = this->boundsSection.front();
    // cases: c < x < c, c < x, x < c, x free,
    if (current->type == LpTokenType::VARIDENTIFIER &&
        next->type == LpTokenType::FREE) {
      HighsVar* variable;
      model.HighsGetOrCreateVarByName(((LpTokenVarIdentifier*)current)->value,
                                      &variable);
      variable->lowerBound = -HIGHS_CONST_INF;
      variable->upperBound = HIGHS_CONST_INF;
      delete current;
      delete next;
    } else if (current->type == LpTokenType::CONSTANT) {
      assert(next->type == LpTokenType::COMPARISON);
      assert(this->boundsSection.size() > 0);
      LpToken* nextnext = this->boundsSection.front();
      this->boundsSection.pop_front();
      assert(nextnext->type == LpTokenType::VARIDENTIFIER);

      assert(((LpTokenComparison*)next)->comparison ==
             LpComparisonIndicator::LEQ);
      HighsVar* variable;
      model.HighsGetOrCreateVarByName(((LpTokenVarIdentifier*)nextnext)->value,
                                      &variable);
      variable->lowerBound = ((LpTokenConstant*)current)->value;

      delete current;
      delete next;

      if (this->boundsSection.size() > 0) {
        LpToken* nextnextnext = this->boundsSection.front();
        if (nextnextnext->type == LpTokenType::COMPARISON) {
          this->boundsSection.push_front(nextnext);
        } else {
          delete nextnext;
        }
      } else {
        delete nextnext;
      }
    } else if (current->type == LpTokenType::VARIDENTIFIER) {
      assert(next->type == LpTokenType::COMPARISON);
      assert(this->boundsSection.size() > 0);
      LpToken* nextnext = this->boundsSection.front();
      this->boundsSection.pop_front();
      assert(nextnext->type == LpTokenType::CONSTANT);

      assert(((LpTokenComparison*)next)->comparison ==
             LpComparisonIndicator::LEQ);
      HighsVar* variable;
      model.HighsGetOrCreateVarByName(((LpTokenVarIdentifier*)current)->value,
                                      &variable);
      variable->upperBound = ((LpTokenConstant*)nextnext)->value;

      delete current;
      delete nextnext;
      delete next;
    } else {
      HighsLogMessage(stdout, HighsMessageType::ERROR,
                      "Error when parsing bounds section.\n");
      this->status = LP_FILEREADER_STATUS::ERROR;
      delete current;
      delete next;
      return;
    }
  }
}

void FilereaderLp::handleConstraintSection(HighsModelBuilder& model) {
  assert(this->constraintSection.size() > 0);

  LpToken* token;
  token = this->constraintSection.front();
  this->constraintSection.pop_front();
  assert(token->type == LpTokenType::SECTIONKEYWORD);
  assert(((LpTokenSectionKeyword*)token)->section == LpSectionKeyword::CON);
  delete token;

  while (this->constraintSection.size() > 0) {
    LpToken* current = this->constraintSection.front();
    // handle next constraint
    HighsLinearCons* constraint;

    // create constraint (with name if available)
    if (current->type == LpTokenType::CONSIDENTIFIER) {
      model.HighsCreateLinearCons(((LpTokenConsIdentifier*)current)->value,
                                  &constraint);
      this->constraintSection.pop_front();
      delete current;
    } else {
      model.HighsCreateLinearCons(&constraint);
    }

    current = this->constraintSection.front();
    while (current->type != LpTokenType::COMPARISON) {
      this->constraintSection.pop_front();
      LpToken* next = this->constraintSection.front();
      if (next->type == LpTokenType::COMPARISON) {
        next = NULL;
      }

      if (current->type == LpTokenType::VARIDENTIFIER &&
          (next == NULL || next->type == LpTokenType::CONSTANT)) {
        // variable with implicit coefficient
        HighsVar* variable;
        HighsLinearConsCoef* coefficient;
        model.HighsGetOrCreateVarByName(((LpTokenVarIdentifier*)current)->value,
                                        &variable);
        model.HighsCreateLinearConsCoef(variable, 1.0, &coefficient);
        model.HighsAddLinearConsCoefToCons(constraint, coefficient);
        delete current;
      } else if (current->type == LpTokenType::CONSTANT &&
                 next->type == LpTokenType::VARIDENTIFIER) {
        // variable with explicit coefficient
        HighsVar* variable;
        HighsLinearConsCoef* coefficient;
        model.HighsGetOrCreateVarByName(((LpTokenVarIdentifier*)next)->value,
                                        &variable);
        model.HighsCreateLinearConsCoef(
            variable, ((LpTokenConstant*)current)->value, &coefficient);
        model.HighsAddLinearConsCoefToCons(constraint, coefficient);
        delete current;
        this->constraintSection.pop_front();
        delete next;
      } else {
        // error
        HighsLogMessage(stdout, HighsMessageType::ERROR,
                        "Error when parsing constraint section\n");
        this->status = LP_FILEREADER_STATUS::ERROR;
        delete current;
        return;
      }
      current = this->constraintSection.front();
    }
    assert(current->type == LpTokenType::COMPARISON);
    this->constraintSection.pop_front();
    assert(this->constraintSection.size() > 0);
    LpToken* next = this->constraintSection.front();
    assert(next->type == LpTokenType::CONSTANT);
    this->constraintSection.pop_front();

    switch (((LpTokenComparison*)current)->comparison) {
      case LpComparisonIndicator::LEQ:
        constraint->upperBound = ((LpTokenConstant*)next)->value;
        break;
      case LpComparisonIndicator::EQ:
        constraint->lowerBound = ((LpTokenConstant*)next)->value;
        constraint->upperBound = ((LpTokenConstant*)next)->value;
        break;
      case LpComparisonIndicator::GEQ:
        constraint->lowerBound = ((LpTokenConstant*)next)->value;
        break;
      default:
        // L or G: ignore for now, maybe error?
        break;
    }
    delete current;
    delete next;
  }
}

void FilereaderLp::handleObjectiveSection(HighsModelBuilder& model) {
  assert(this->objectiveSection.size() > 0);

  LpToken* token;
  // handle objective sense
  token = this->objectiveSection.front();
  this->objectiveSection.pop_front();
  assert(token->type == LpTokenType::SECTIONKEYWORD);
  assert(((LpTokenSectionKeyword*)token)->section == LpSectionKeyword::OBJ);

  if (((LpTokenObjectiveSectionKeyword*)token)->objectiveType !=
      LpObjectiveSectionKeywordType::MIN) {
    assert(((LpTokenObjectiveSectionKeyword*)token)->objectiveType ==
           LpObjectiveSectionKeywordType::MAX);
    model.objSense = ObjSense::MAXIMIZE;
  }
  delete token;

  if (this->objectiveSection.size() == 0) {
    return;
  }

  // handle objective name
  token = this->objectiveSection.front();
  if (token->type == LpTokenType::CONSIDENTIFIER) {
    this->objectiveSection.pop_front();
    delete token;
  }

  // 3 cases: constant+varidentifier, constant, varidentifier, later: quadratic
  // terms
  while (this->objectiveSection.size() > 0) {
    LpToken* current = this->objectiveSection.front();
    this->objectiveSection.pop_front();
    LpToken* next = NULL;
    if (this->objectiveSection.size() > 0) {
      next = this->objectiveSection.front();
    }

    if (current->type == LpTokenType::CONSTANT &&
        (next == NULL || next->type == LpTokenType::CONSTANT)) {
      // standalone constanst aka objective offset
      model.objOffset = ((LpTokenConstant*)current)->value;
      delete current;
    } else if (current->type == LpTokenType::CONSTANT &&
               next->type == LpTokenType::VARIDENTIFIER) {
      // variable with constant
      this->objectiveSection.pop_front();
      HighsVar* variable;
      model.HighsGetOrCreateVarByName(((LpTokenVarIdentifier*)next)->value,
                                      &variable);
      variable->obj = ((LpTokenConstant*)current)->value;
      delete current;
      delete next;
    } else if (current->type == LpTokenType::VARIDENTIFIER) {
      // variable with implied constant
      HighsVar* variable;
      model.HighsGetOrCreateVarByName(((LpTokenVarIdentifier*)current)->value,
                                      &variable);
      variable->obj = 1;
      delete current;
    } else {
      // error
      HighsLogMessage(stdout, HighsMessageType::ERROR,
                      "Error when parsing objective section.\n");
      this->status = LP_FILEREADER_STATUS::ERROR;
      delete current;
      return;
    }
  }
}

void FilereaderLp::splitTokens() {
  std::list<LpToken*>* dest = NULL;
  while (this->tokenQueue.size() > 0) {
    LpToken* token = this->tokenQueue.front();
    assert(token->type == LpTokenType::SECTIONKEYWORD);
    LpTokenSectionKeyword* sectionKeyword = (LpTokenSectionKeyword*)token;
    switch (sectionKeyword->section) {
      case LpSectionKeyword::OBJ:
        dest = &this->objectiveSection;
        break;
      case LpSectionKeyword::CON:
        dest = &this->constraintSection;
        break;
      case LpSectionKeyword::BOUNDS:
        dest = &this->boundsSection;
        break;
      case LpSectionKeyword::BIN:
        dest = &this->binSection;
        break;
      case LpSectionKeyword::GEN:
        dest = &this->generalSection;
        break;
      case LpSectionKeyword::SEMI:
        dest = &this->semiSection;
        break;
      case LpSectionKeyword::SOS:
        dest = &this->sosSection;
        break;
      case LpSectionKeyword::END:
        this->tokenQueue.pop_front();
        delete token;
        return;
      case LpSectionKeyword::NONE:
        // error
        this->status = LP_FILEREADER_STATUS::ERROR;
        HighsLogMessage(stdout, HighsMessageType::ERROR,
                        "Error when splitting tokens.\n");
        return;
    }
    do {
      this->tokenQueue.pop_front();
      dest->push_back(token);
      token = this->tokenQueue.front();
    } while (token != NULL && token->type != LpTokenType::SECTIONKEYWORD);
  }
}

FilereaderRetcode FilereaderLp::tokenizeInput() {
  // add extra new line at beginning (to identify initial LpSectionKeyword
  // properly)
  LpToken* newToken = new LpToken(LpTokenType::LINEEND);
  this->tokenQueue.push_back(newToken);
  bool cont;
  do {
    cont = this->readNextToken();
  } while (cont);

  return FilereaderRetcode::OK;
}

LpSectionKeyword FilereaderLp::tryParseLongSectionKeyword(const char* str,
                                                          int* characters) {
  char s1[LP_MAX_NAME_LENGTH];
  char s2[LP_MAX_NAME_LENGTH];
  char s3[LP_MAX_LINE_LENGTH];

  int nread = sscanf(str, "%s %s%n", s1, s2, characters);
  if (nread == 2) {
    sprintf(s3, "%s %s", s1, s2);
    char* s4 = strClone(s3);
    strToLower(s4);
    if (strcmp(s4, LP_KEYWORD_ST[0]) == 0) {
      return LpSectionKeyword::CON;
    }
    if (strcmp(s4, LP_KEYWORD_ST[1]) == 0) {
      return LpSectionKeyword::CON;
    }
  }

  nread = sscanf(str, "%s%n", s1, characters);
  if (nread == 1) {
    if (strcmp(s1, LP_KEYWORD_SEMI[0]) == 0) {
      return LpSectionKeyword::SEMI;
    }
  }

  return LpSectionKeyword::NONE;
}

bool FilereaderLp::readNextToken() {
  LpToken* previousToken = this->tokenQueue.back();
  bool previousTokenWasLineEnd = false;
  if (previousToken->type == LpTokenType::LINEEND) {
    previousTokenWasLineEnd = true;
  }

  // check if we need to read a new line of the file
  if (this->isFileBufferFullyRead) {
    char* eof = fgets(this->fileBuffer, BUFFERSIZE, this->file);

    // check if the file ended
    if (eof == NULL) {
      if (previousTokenWasLineEnd) {
        this->tokenQueue.pop_back();
        delete previousToken;
      }
      return false;
    } else {
      this->isFileBufferFullyRead = false;
      this->readingPosition = this->fileBuffer;
    }
  }

  // check if a comment starts at the current reading position
  if (*this->readingPosition == '\\') {
    if (!previousTokenWasLineEnd) {
      LpToken* newToken = new LpToken(LpTokenType::LINEEND);
      this->tokenQueue.push_back(newToken);
    }
    this->isFileBufferFullyRead = true;
    return true;
  }

  // check if a keyword containing a space/hyphen starts the the current reading
  // position
  int charactersConsumed;
  LpSectionKeyword longKeyword = this->tryParseLongSectionKeyword(
      this->readingPosition, &charactersConsumed);
  if (previousTokenWasLineEnd && longKeyword != LpSectionKeyword::NONE) {
    LpTokenSectionKeyword* newToken = new LpTokenSectionKeyword(longKeyword);
    this->tokenQueue.pop_back();
    delete previousToken;
    this->tokenQueue.push_back(newToken);
    this->readingPosition += charactersConsumed;
    return true;
  }

  // check if a (standard) constant starts at the current reading position
  int nread = sscanf(this->readingPosition, "%lf%n", &this->constantBuffer,
                     &charactersConsumed);
  if (nread == 1) {
    if (this->constantBuffer >= HIGHS_CONST_INF) {
      this->constantBuffer = HIGHS_CONST_INF;
    }
    int multiplier = 1;
    if (previousToken->type == LpTokenType::SIGN) {
      this->tokenQueue.pop_back();
      multiplier = ((LpTokenSign*)previousToken)->sign;
    }
    this->readingPosition += charactersConsumed;
    LpTokenConstant* newToken =
        new LpTokenConstant(this->constantBuffer * multiplier);
    if (previousTokenWasLineEnd) {
      this->tokenQueue.pop_back();
      delete previousToken;
    }
    this->tokenQueue.push_back(newToken);
    return true;
  }

  // check if a non-standard constant starts at the current reading position
  // TODO

  // read string, check if it is a keyword, a variable name, a constraint name,
  // 'free', or infinity (constant)
  nread = sscanf(this->readingPosition, "%[^][\t\n:+<>= -]%n",
                 this->stringBuffer, &charactersConsumed);
  if (nread == 1) {
    // check if it is a section keyword
    LpSectionKeyword keyword = this->tryParseSectionKeyword(this->stringBuffer);
    if (previousTokenWasLineEnd && keyword != LpSectionKeyword::NONE) {
      LpTokenSectionKeyword* newToken;
      if (keyword == LpSectionKeyword::OBJ) {
        newToken = new LpTokenObjectiveSectionKeyword(
            this->parseObjectiveSectionKeyword(this->stringBuffer));
      } else {
        newToken = new LpTokenSectionKeyword(keyword);
      }
      this->tokenQueue.pop_back();
      delete previousToken;

      this->tokenQueue.push_back(newToken);
      this->readingPosition += charactersConsumed;

      return true;
    }

    // check if it is infinity
    bool isInfinity =
        this->isKeyword(this->stringBuffer, LP_KEYWORD_INF, LP_KEYWORD_INF_N);
    if (isInfinity) {
      int multiplier = 1;
      if (!previousTokenWasLineEnd &&
          previousToken->type == LpTokenType::SIGN) {
        this->tokenQueue.pop_back();
        multiplier = ((LpTokenSign*)previousToken)->sign;
      }
      LpTokenConstant* newToken =
          new LpTokenConstant(HIGHS_CONST_INF * multiplier);
      if (previousTokenWasLineEnd) {
        this->tokenQueue.pop_back();
        delete previousToken;
      }
      this->tokenQueue.push_back(newToken);
      this->readingPosition += charactersConsumed;
      return true;
    }

    // check if it is 'free'
    bool isFree =
        this->isKeyword(this->stringBuffer, LP_KEYWORD_FREE, LP_KEYWORD_FREE_N);
    if (isFree) {
      LpToken* newToken = new LpToken(LpTokenType::FREE);
      if (previousTokenWasLineEnd) {
        // should not happen
        this->tokenQueue.pop_back();
        delete previousToken;
        HighsLogMessage(stdout, HighsMessageType::ERROR,
                        "Error when parsing file.\n");
        this->status = LP_FILEREADER_STATUS::ERROR;
      }
      this->tokenQueue.push_back(newToken);
      this->readingPosition += charactersConsumed;
      return true;
    }

    // check if it is a constraint name
    // TODO: check if name is allowed
    if (previousTokenWasLineEnd &&
        *(this->readingPosition + charactersConsumed) == ':') {
      LpTokenConsIdentifier* newToken =
          new LpTokenConsIdentifier(this->stringBuffer);
      if (previousTokenWasLineEnd) {
        this->tokenQueue.pop_back();
        delete previousToken;
      }
      this->tokenQueue.push_back(newToken);
      this->readingPosition += charactersConsumed + 1;
      return true;
    }

    // check if it is a variable name
    // TODO: check if name is allowed
    if (!previousTokenWasLineEnd && previousToken->type == LpTokenType::SIGN) {
      this->tokenQueue.pop_back();
      LpTokenConstant* newToken =
          new LpTokenConstant(((LpTokenSign*)previousToken)->sign);
      this->tokenQueue.push_back(newToken);
    }
    LpTokenVarIdentifier* newToken =
        new LpTokenVarIdentifier(this->stringBuffer);
    if (previousTokenWasLineEnd) {
      this->tokenQueue.pop_back();
      delete previousToken;
    }
    this->tokenQueue.push_back(newToken);
    this->readingPosition += charactersConsumed;
    return true;
  }

  // read single character, check if it is a lineend, whitespace (tab or space),
  // (partial) comparison, colon (should not happen), sign, or bracket
  if (*this->readingPosition == '\0') {
    HighsLogMessage(stdout, HighsMessageType::ERROR,
                    "NULL character read. Should not have happened.\n");
    this->isFileBufferFullyRead = true;
    this->status = LP_FILEREADER_STATUS::ERROR;
    return false;
  }

  char symbol;
  nread = sscanf(this->readingPosition, "%c", &symbol);
  if (nread == 1) {
    LpToken* newToken;
    switch (symbol) {
      case '\n':
        if (!previousTokenWasLineEnd) {
          newToken = new LpToken(LpTokenType::LINEEND);
          this->tokenQueue.push_back(newToken);
        }
        this->isFileBufferFullyRead = true;
        return true;

      case ':':
        HighsLogMessage(stdout, HighsMessageType::ERROR,
                        "COLON character read. Should not have happened.\n");
        this->readingPosition += 1;
        this->status = LP_FILEREADER_STATUS::ERROR;
        return false;

      case ' ':
      case '\t':
        this->readingPosition += 1;
        return previousToken;

      case '+':
        newToken = new LpTokenSign(1);
        if (previousTokenWasLineEnd) {
          this->tokenQueue.pop_back();
          delete previousToken;
        }
        this->tokenQueue.push_back(newToken);
        this->readingPosition += 1;
        return true;

      case '-':
        newToken = new LpTokenSign(-1);
        if (previousTokenWasLineEnd) {
          this->tokenQueue.pop_back();
          delete previousToken;
        }
        this->tokenQueue.push_back(newToken);
        this->readingPosition += 1;
        return true;

      case '[':
        newToken = new LpToken(LpTokenType::BRACKETOPEN);
        if (previousTokenWasLineEnd) {
          this->tokenQueue.pop_back();
          delete previousToken;
        }
        this->tokenQueue.push_back(newToken);
        this->readingPosition += 1;
        return true;

      case ']':
        newToken = new LpToken(LpTokenType::BRACKETCLOSE);
        if (previousTokenWasLineEnd) {
          this->tokenQueue.pop_back();
          delete previousToken;
        }
        this->tokenQueue.push_back(newToken);
        this->readingPosition += 1;
        return true;

      case '<':
        if (previousTokenWasLineEnd) {
          this->tokenQueue.pop_back();
          delete previousToken;
        }
        if (!previousTokenWasLineEnd &&
            previousToken->type == LpTokenType::COMPARISON) {
          ((LpTokenComparison*)previousToken)
              ->upgrade(LpComparisonIndicator::L);
          this->readingPosition += 1;

          return true;
        } else {
          newToken = new LpTokenComparison(LpComparisonIndicator::L);
          this->readingPosition += 1;
          this->tokenQueue.push_back(newToken);
          return true;
        }

      case '>':
        if (previousTokenWasLineEnd) {
          this->tokenQueue.pop_back();
          delete previousToken;
        }
        if (!previousTokenWasLineEnd &&
            previousToken->type == LpTokenType::COMPARISON) {
          ((LpTokenComparison*)previousToken)
              ->upgrade(LpComparisonIndicator::G);
          this->readingPosition += 1;

          return true;
        } else {
          newToken = new LpTokenComparison(LpComparisonIndicator::G);
          this->readingPosition += 1;
          this->tokenQueue.push_back(newToken);
          return true;
        }

      case '=':
        if (previousTokenWasLineEnd) {
          this->tokenQueue.pop_back();
          delete previousToken;
        }
        if (!previousTokenWasLineEnd &&
            previousToken->type == LpTokenType::COMPARISON) {
          ((LpTokenComparison*)previousToken)
              ->upgrade(LpComparisonIndicator::EQ);
          this->readingPosition += 1;

          return true;
        } else {
          newToken = new LpTokenComparison(LpComparisonIndicator::EQ);
          this->readingPosition += 1;
          this->tokenQueue.push_back(newToken);
          return true;
        }
      default:
        HighsLogMessage(stdout, HighsMessageType::ERROR, "Unknown symbol\n");
        return false;
    }
  }
  return false;
}

bool FilereaderLp::isKeyword(const char* orig, const char* const* keywords,
                             const int nkeywords) {
  char* str = strClone(orig);
  strToLower(str);
  int i;
  for (i = 0; i < nkeywords; i++) {
    if (strcmp(str, keywords[i]) == 0) {
      delete[] str;
      return true;
    }
  }
  delete[] str;
  return false;
}

LpObjectiveSectionKeywordType FilereaderLp::parseObjectiveSectionKeyword(
    const char* str) {
  // min?
  if (isKeyword(str, LP_KEYWORD_MIN, LP_KEYWORD_MIN_N)) {
    return LpObjectiveSectionKeywordType::MIN;
  }

  // max?
  if (isKeyword(str, LP_KEYWORD_MAX, LP_KEYWORD_MAX_N)) {
    return LpObjectiveSectionKeywordType::MAX;
  }

  return LpObjectiveSectionKeywordType::NONE;
}

LpSectionKeyword FilereaderLp::tryParseSectionKeyword(const char* str) {
  // min?
  if (isKeyword(str, LP_KEYWORD_MIN, LP_KEYWORD_MIN_N)) {
    return LpSectionKeyword::OBJ;
  }

  // max?
  if (isKeyword(str, LP_KEYWORD_MAX, LP_KEYWORD_MAX_N)) {
    return LpSectionKeyword::OBJ;
  }

  // st?
  if (isKeyword(str, LP_KEYWORD_ST, LP_KEYWORD_ST_N)) {
    return LpSectionKeyword::CON;
  }

  // bounds?
  if (isKeyword(str, LP_KEYWORD_BOUNDS, LP_KEYWORD_BOUNDS_N)) {
    return LpSectionKeyword::BOUNDS;
  }

  // gen?
  if (isKeyword(str, LP_KEYWORD_GEN, LP_KEYWORD_GEN_N)) {
    return LpSectionKeyword::GEN;
  }

  // bin?
  if (isKeyword(str, LP_KEYWORD_BIN, LP_KEYWORD_BIN_N)) {
    return LpSectionKeyword::BIN;
  }

  // semi?
  if (isKeyword(str, LP_KEYWORD_SEMI, LP_KEYWORD_SEMI_N)) {
    return LpSectionKeyword::SEMI;
  }

  // sos?
  if (isKeyword(str, LP_KEYWORD_SOS, LP_KEYWORD_SOS_N)) {
    return LpSectionKeyword::SOS;
  }

  // end?
  if (isKeyword(str, LP_KEYWORD_END, LP_KEYWORD_END_N)) {
    return LpSectionKeyword::END;
  }

  return LpSectionKeyword::NONE;
}

void FilereaderLp::writeToFile(const char* format, ...) {
  va_list argptr;
  va_start(argptr, format);
  int tokenlength = vsprintf(this->stringBuffer, format, argptr);
  if (this->linelength + tokenlength >= LP_MAX_LINE_LENGTH) {
    fprintf(this->file, "\n");
    fprintf(this->file, "%s", this->stringBuffer);
    this->linelength = tokenlength;
  } else {
    fprintf(this->file, "%s", this->stringBuffer);
    this->linelength += tokenlength;
  }
}

void FilereaderLp::writeToFileLineend() {
  fprintf(this->file, "\n");
  this->linelength = 0;
}

HighsStatus FilereaderLp::writeModelToFile(const HighsOptions& options,
                                           const char* filename,
                                           HighsLp& model) {
  this->file = fopen(filename, "w");

  // write comment at the start of the file
  this->writeToFile("\\ %s", LP_COMMENT_FILESTART);
  this->writeToFileLineend();

  // write objective
  this->writeToFile("%s", model.sense_ == ObjSense::MINIMIZE
                              ? LP_KEYWORD_MIN[0]
                              : LP_KEYWORD_MAX[0]);
  this->writeToFileLineend();
  this->writeToFile(" obj: ");
  for (int i = 0; i < model.numCol_; i++) {
    this->writeToFile("%+g x%d ", model.colCost_[i], (i + 1));
  }
  this->writeToFileLineend();

  // write constraint section, lower & upper bounds are one constraint each
  this->writeToFile("%s", LP_KEYWORD_ST[2]);
  this->writeToFileLineend();
  for (int row = 0; row < model.numRow_; row++) {
    if (model.rowLower_[row] == model.rowUpper_[row]) {
      // equality constraint
      this->writeToFile(" con%d: ", row + 1);
      for (int var = 0; var < model.numCol_; var++) {
        for (int idx = model.Astart_[var]; idx < model.Astart_[var + 1];
             idx++) {
          if (model.Aindex_[idx] == row) {
            this->writeToFile("%+g x%d ", model.Avalue_[idx], var + 1);
          }
        }
      }
      this->writeToFile("= %+g", model.rowLower_[row]);
      this->writeToFileLineend();
    } else {
      if (model.rowLower_[row] >= -10E10) {
        // has a lower bounds
        this->writeToFile(" con%dlo: ", row + 1);
        for (int var = 0; var < model.numCol_; var++) {
          for (int idx = model.Astart_[var]; idx < model.Astart_[var + 1];
               idx++) {
            if (model.Aindex_[idx] == row) {
              this->writeToFile("%+g x%d ", model.Avalue_[idx], var + 1);
            }
          }
        }
        this->writeToFile(">= %+g", model.rowLower_[row]);
        this->writeToFileLineend();
      } else if (model.rowUpper_[row] <= 10E10) {
        // has an upper bounds
        this->writeToFile(" con%dup: ", row + 1);
        for (int var = 0; var < model.numCol_; var++) {
          for (int idx = model.Astart_[var]; idx < model.Astart_[var + 1];
               idx++) {
            if (model.Aindex_[idx] == row) {
              this->writeToFile("%+g x%d ", model.Avalue_[idx], var + 1);
            }
          }
        }
        this->writeToFile("<= %+g", model.rowLower_[row]);
        this->writeToFileLineend();
      } else {
        // constraint has infinite lower & upper bounds so not a proper
        // constraint, does not get written
      }
    }
  }

  // write bounds section
  this->writeToFile("%s", LP_KEYWORD_BOUNDS[0]);
  this->writeToFileLineend();
  for (int i = 0; i < model.numCol_; i++) {
    // if both lower/upper bound are +/-infinite: [name] free
    if (model.colLower_[i] > -HIGHS_CONST_INF &&
        model.colUpper_[i] < HIGHS_CONST_INF) {
      this->writeToFile(" %+g <= x%d <= %+g", model.colLower_[i], i + 1,
                        model.colUpper_[i]);
      this->writeToFileLineend();
    } else if (model.colLower_[i] <= -HIGHS_CONST_INF &&
               model.colUpper_[i] < HIGHS_CONST_INF) {
      this->writeToFile(" -inf <= x%d <= %+g", i + 1, model.colUpper_[i]);
      this->writeToFileLineend();

    } else if (model.colLower_[i] > -HIGHS_CONST_INF &&
               model.colUpper_[i] >= HIGHS_CONST_INF) {
      this->writeToFile(" %+g <= x%d <= +inf", model.colLower_[i], i + 1);
      this->writeToFileLineend();
    } else {
      this->writeToFile(" x%d %s", i + 1, LP_KEYWORD_FREE[0]);
      this->writeToFileLineend();
    }
  }

  // write binary section
  this->writeToFile("%s", LP_KEYWORD_BIN[0]);
  this->writeToFileLineend();

  // write general section
  this->writeToFile("%s", LP_KEYWORD_GEN[0]);
  this->writeToFileLineend();

  // write semi section
  this->writeToFile("%s", LP_KEYWORD_SEMI[1]);
  this->writeToFileLineend();

  // write end
  this->writeToFile("%s", LP_KEYWORD_END[0]);
  this->writeToFileLineend();

  fclose(this->file);
  return HighsStatus::OK;
}
