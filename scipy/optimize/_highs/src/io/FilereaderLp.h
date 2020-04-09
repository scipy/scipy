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

#ifndef IO_FILEREADER_LP_H_
#define IO_FILEREADER_LP_H_

#include <list>

#include "io/Filereader.h"
#include "io/HighsIO.h"

#define BUFFERSIZE 561
#define LP_MAX_LINE_LENGTH 560
#define LP_MAX_NAME_LENGTH 255

enum class LP_FILEREADER_STATUS { NONE, SUCCESS, ERROR };

#define LP_COMMENT_FILESTART ("File written by Highs .lp filereader")

const char* const LP_KEYWORD_MIN[] = {"minimize", "min", "minimum"};
const char* const LP_KEYWORD_MAX[] = {"maximize", "max", "maximum"};
const char* const LP_KEYWORD_ST[] = {"subject to", "such that", "st", "s.t."};
const char* const LP_KEYWORD_BOUNDS[] = {"bounds", "bound"};
const char* const LP_KEYWORD_INF[] = {"infinity", "inf"};
const char* const LP_KEYWORD_FREE[] = {"free"};
const char* const LP_KEYWORD_GEN[] = {"general", "generals", "gen"};
const char* const LP_KEYWORD_BIN[] = {"binary", "binaries", "bin"};
const char* const LP_KEYWORD_SEMI[] = {"semi-continuous", "semi", "semis"};
const char* const LP_KEYWORD_SOS[] = {"sos"};
const char* const LP_KEYWORD_END[] = {"end"};

const int LP_KEYWORD_MIN_N = 3;
const int LP_KEYWORD_MAX_N = 3;
const int LP_KEYWORD_ST_N = 4;
const int LP_KEYWORD_BOUNDS_N = 2;
const int LP_KEYWORD_INF_N = 2;
const int LP_KEYWORD_FREE_N = 1;
const int LP_KEYWORD_GEN_N = 3;
const int LP_KEYWORD_BIN_N = 3;
const int LP_KEYWORD_SEMI_N = 3;
const int LP_KEYWORD_SOS_N = 1;
const int LP_KEYWORD_END_N = 1;

enum class LpSectionKeyword {
  NONE,
  OBJ,
  CON,
  BOUNDS,
  GEN,
  BIN,
  SEMI,
  SOS,
  END
};

enum class LpObjectiveSectionKeywordType { NONE, MIN, MAX };

enum class LpComparisonIndicator { LEQ, L, EQ, G, GEQ };

enum LpTokenType {
  NONE,
  VARIDENTIFIER,
  CONSIDENTIFIER,
  SECTIONKEYWORD,
  FREE,
  CONSTANT,
  SIGN,
  COLON,
  BRACKETOPEN,
  BRACKETCLOSE,
  COMPARISON,
  LINEEND,
  FILEEND
};

const char* const LpTokenTypeString[] = {
    "NONE",        "VARIDENTIFIER", "CONSIDENTIFIER", "SECTIONKEYWORD",
    "FREE",        "CONSTANT",      "SIGN",           "COLON",
    "BRACKETOPEN", "BRACKETCLOSE",  "COMPARISON",     "LINEEND",
    "FILEEND"};

const char* const LpSectionKeywordString[] = {
    "NONE", "OBJ", "ST", "BOUNDS", "GEN", "BIN", "SEMI", "SOS", "END"};

const char* const LpObjectiveSectionKeywordString[] = {"NONE", "MIN", "MAX"};

class LpToken {
 public:
  LpTokenType type;
  virtual void print() {
    HighsLogMessage(stdout, HighsMessageType::INFO, "%s ",
                    LpTokenTypeString[type]);
  }

  virtual ~LpToken() { ; }

  LpToken(LpTokenType t) { this->type = t; }
};

class LpTokenVarIdentifier : public LpToken {
 public:
  char* value;

  LpTokenVarIdentifier(char* v) : LpToken(LpTokenType::VARIDENTIFIER) {
    int len;
    len = strlen(v);
    this->value = new char[len + 1];
    strcpy(this->value, v);
  }

  ~LpTokenVarIdentifier() { delete[] this->value; }
};

class LpTokenConsIdentifier : public LpToken {
 public:
  char* value;

  LpTokenConsIdentifier(char* v) : LpToken(LpTokenType::CONSIDENTIFIER) {
    int len;
    len = strlen(v);
    this->value = new char[len + 1];
    strcpy(this->value, v);
  }

  ~LpTokenConsIdentifier() { delete[] this->value; }
};

class LpTokenSectionKeyword : public LpToken {
 public:
  LpSectionKeyword section;
  LpTokenSectionKeyword(LpSectionKeyword k)
      : LpToken(LpTokenType::SECTIONKEYWORD) {
    this->type = LpTokenType::SECTIONKEYWORD;
    this->section = k;
  }
};

class LpTokenObjectiveSectionKeyword : public LpTokenSectionKeyword {
 public:
  LpObjectiveSectionKeywordType objectiveType =
      LpObjectiveSectionKeywordType::MIN;
  LpTokenObjectiveSectionKeyword()
      : LpTokenSectionKeyword(LpSectionKeyword::OBJ) {}
  LpTokenObjectiveSectionKeyword(LpObjectiveSectionKeywordType t)
      : LpTokenSectionKeyword(LpSectionKeyword::OBJ) {
    this->objectiveType = t;
  }
};

class LpTokenConstant : public LpToken {
 public:
  double value;
  bool isPositiveInfinity;
  bool isNegativeInfinity;

  LpTokenConstant(double v) : LpToken(LpTokenType::CONSTANT) {
    this->value = v;
  }
};

class LpTokenSign : public LpToken {
 public:
  int sign;
  LpTokenSign(int s) : LpToken(LpTokenType::SIGN) { this->sign = s; }
};

class LpTokenComparison : public LpToken {
 public:
  LpComparisonIndicator comparison;

  LpTokenComparison(LpComparisonIndicator c)
      : LpToken(LpTokenType::COMPARISON) {
    this->comparison = c;
  }

  void upgrade(LpComparisonIndicator c) {
    switch (this->comparison) {
      case LpComparisonIndicator::G:
        if (c == LpComparisonIndicator::EQ) {
          this->comparison = LpComparisonIndicator::GEQ;
        } else {
          // error
          HighsLogMessage(stdout, HighsMessageType::ERROR,
                          "Invalid comparison indicator.\n");
        }
        break;
      case LpComparisonIndicator::L:
        if (c == LpComparisonIndicator::EQ) {
          this->comparison = LpComparisonIndicator::LEQ;
        } else {
          // error
          HighsLogMessage(stdout, HighsMessageType::ERROR,
                          "Invalid comparison indicator.\n");
        }
        break;
      case LpComparisonIndicator::EQ:
        if (c == LpComparisonIndicator::EQ) {
          // do nothing
        } else if (c == LpComparisonIndicator::G) {
          this->comparison = LpComparisonIndicator::GEQ;
        } else if (c == LpComparisonIndicator::L) {
          this->comparison = LpComparisonIndicator::LEQ;
        } else {
          // error
          HighsLogMessage(stdout, HighsMessageType::ERROR,
                          "Invalid comparison indicator.\n");
        }
        break;
      default:
        // error
        HighsLogMessage(stdout, HighsMessageType::ERROR,
                        "Invalid comparison indicator.\n");
    }
  }
};

class FilereaderLp : public Filereader {
 public:
  FilereaderRetcode readModelFromFile(const HighsOptions& options,
                                      HighsLp& model);
  FilereaderRetcode readModelFromFile(const char* filename,
                                      HighsModelBuilder& model);
  HighsStatus writeModelToFile(const HighsOptions& options,
                               const char* filename, HighsLp& model);
  FilereaderLp();
  ~FilereaderLp();

 private:
  // list of all tokens after initial scan of file
  std::list<LpToken*> tokenQueue;

  // tokens split according to their section
  std::list<LpToken*> objectiveSection;
  std::list<LpToken*> constraintSection;
  std::list<LpToken*> boundsSection;
  std::list<LpToken*> binSection;
  std::list<LpToken*> generalSection;
  std::list<LpToken*> sosSection;
  std::list<LpToken*> semiSection;

  FILE* file;
  char fileBuffer[BUFFERSIZE];
  char stringBuffer[BUFFERSIZE];
  char* readingPosition;
  bool isFileBufferFullyRead;
  double constantBuffer;

  // functions to read files
  bool isKeyword(const char* str, const char* const* keywords,
                 const int nkeywords);
  LpSectionKeyword tryParseLongSectionKeyword(const char* str, int* characters);
  LpSectionKeyword tryParseSectionKeyword(const char* str);
  LpObjectiveSectionKeywordType parseObjectiveSectionKeyword(const char* str);

  FilereaderRetcode tokenizeInput();
  bool readNextToken();
  void splitTokens();

  void handleObjectiveSection(HighsModelBuilder& model);
  void handleConstraintSection(HighsModelBuilder& model);
  void handleBoundsSection(HighsModelBuilder& model);
  void handleBinarySection(HighsModelBuilder& model);
  void handleGeneralSection(HighsModelBuilder& model);
  void handleSemiSection(HighsModelBuilder& model);
  FilereaderRetcode handleSosSection(HighsModelBuilder& model);

  LP_FILEREADER_STATUS status;

  // functions to write files
  int linelength;
  void writeToFile(const char* format, ...);
  void writeToFileLineend();
};

#endif
