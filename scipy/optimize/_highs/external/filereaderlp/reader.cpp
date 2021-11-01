#include "reader.hpp"

#include "builder.hpp"

#include <cstdio>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

void inline lpassert(bool condition) {
   if (!condition) {
      throw std::invalid_argument("File not existant or illegal file format.");
   }
}

#define LP_MAX_NAME_LENGTH 255
#define LP_MAX_LINE_LENGTH 560

const std::string LP_KEYWORD_MIN[] = {"minimize", "min", "minimum"};
const std::string LP_KEYWORD_MAX[] = {"maximize", "max", "maximum"};
const std::string LP_KEYWORD_ST[] = {"subject to", "such that", "st", "s.t."};
const std::string LP_KEYWORD_BOUNDS[] = {"bounds", "bound"};
const std::string LP_KEYWORD_INF[] = {"infinity", "inf"};
const std::string LP_KEYWORD_FREE[] = {"free"};
const std::string LP_KEYWORD_GEN[] = {"general", "generals", "gen"};
const std::string LP_KEYWORD_BIN[] = {"binary", "binaries", "bin"};
const std::string LP_KEYWORD_SEMI[] = {"semi-continuous", "semi", "semis"};
const std::string LP_KEYWORD_SOS[] = {"sos"};
const std::string LP_KEYWORD_END[] = {"end"};

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

enum class RawTokenType {
   NONE,
   STR,
   CONS,
   LESS,
   GREATER,
   EQUAL,
   COLON,
   LNEND,
   FLEND,
   BRKOP,
   BRKCL,
   PLUS,
   MINUS,
   HAT,
   SLASH,
   ASTERISK
};

struct RawToken {
   RawTokenType type;
   inline bool istype(RawTokenType t) {
      return this->type == t;
   }
   RawToken(RawTokenType t) : type(t) {} ;
};

struct RawStringToken : RawToken {
   std::string value;
   RawStringToken(std::string v) : RawToken(RawTokenType::STR), value(v) {};
};

struct RawConstantToken : RawToken {
   double value;
   RawConstantToken(double v) : RawToken(RawTokenType::CONS), value(v) {};
};

enum class ProcessedTokenType {
   NONE,
   SECID,
   VARID,
   CONID,
   CONST,
   FREE,
   BRKOP,
   BRKCL,
   COMP,
   LNEND,
   SLASH,
   ASTERISK,
   HAT
};

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

enum class LpComparisonType { LEQ, L, EQ, G, GEQ };

struct ProcessedToken {
   ProcessedTokenType type;
   ProcessedToken(ProcessedTokenType t) : type(t) {};
};

struct ProcessedTokenSectionKeyword : ProcessedToken {
   LpSectionKeyword keyword;
   ProcessedTokenSectionKeyword(LpSectionKeyword k) : ProcessedToken(ProcessedTokenType::SECID), keyword(k) {};
};

struct ProcessedTokenObjectiveSectionKeyword : ProcessedTokenSectionKeyword {
   LpObjectiveSectionKeywordType objsense;
   ProcessedTokenObjectiveSectionKeyword(LpObjectiveSectionKeywordType os) : ProcessedTokenSectionKeyword(LpSectionKeyword::OBJ), objsense(os) {};
};

struct ProcessedConsIdToken : ProcessedToken {
   std::string name;
   ProcessedConsIdToken(std::string n) : ProcessedToken(ProcessedTokenType::CONID), name(n) {};
};

struct ProcessedVarIdToken : ProcessedToken {
   std::string name;
   ProcessedVarIdToken(std::string n) : ProcessedToken(ProcessedTokenType::VARID), name(n) {};
};

struct ProcessedConstantToken : ProcessedToken {
   double value;
   ProcessedConstantToken(double v) : ProcessedToken(ProcessedTokenType::CONST), value(v) {};
};

struct ProcessedComparisonToken : ProcessedToken {
   LpComparisonType dir;
   ProcessedComparisonToken(LpComparisonType d) : ProcessedToken(ProcessedTokenType::COMP), dir(d) {};
};

class Reader {
private:
   FILE* file;
   std::vector<std::unique_ptr<RawToken>> rawtokens;
   std::vector<std::unique_ptr<ProcessedToken>> processedtokens;
   std::map<LpSectionKeyword, std::vector<std::unique_ptr<ProcessedToken>>> sectiontokens;
   
   char linebuffer[LP_MAX_LINE_LENGTH+1];
   bool linebufferrefill;
   char* linebufferpos;

   Builder builder;

   void tokenize();
   void readnexttoken(bool& done);
   void processtokens();
   void splittokens();
   void processsections();
   void processnonesec();
   void processobjsec();
   void processconsec();
   void processboundssec();
   void processbinsec();
   void processgensec();
   void processsemisec();
   void processsossec();
   void processendsec();
   void parseexpression(std::vector<std::unique_ptr<ProcessedToken>>& tokens, std::shared_ptr<Expression> expr, unsigned int& i);

public:
   Reader(std::string filename) : file(fopen(filename.c_str(), "r")) {
      lpassert(file != NULL);
   };

   ~Reader() {
      fclose(file);
   }

   Model read();
};

Model readinstance(std::string filename) {
   Reader reader(filename);
   return reader.read();
}

bool isstrequalnocase(const std::string str1, const std::string str2) {
   unsigned int len = str1.size();
    if (str2.size() != len)
        return false;
    for (unsigned int i = 0; i < len; ++i)
        if (tolower(str1[i]) != tolower(str2[i]))
            return false;
    return true;
}

bool iskeyword(const std::string str, const std::string* keywords, const int nkeywords) {
   for (int i=0; i<nkeywords; i++) {
      if (isstrequalnocase(str, keywords[i])) {
         return true;
      }
   }
   return false;
}

LpObjectiveSectionKeywordType parseobjectivesectionkeyword(const std::string str) {
   if (iskeyword(str, LP_KEYWORD_MIN, LP_KEYWORD_MIN_N)) {
      return LpObjectiveSectionKeywordType::MIN;
   }

   if (iskeyword(str, LP_KEYWORD_MAX, LP_KEYWORD_MAX_N)) {
      return LpObjectiveSectionKeywordType::MAX;
   }
   
   return LpObjectiveSectionKeywordType::NONE;
}

LpSectionKeyword parsesectionkeyword(const std::string& str) {
   if (parseobjectivesectionkeyword(str) != LpObjectiveSectionKeywordType::NONE) {
      return LpSectionKeyword::OBJ;
   }

   if (iskeyword(str, LP_KEYWORD_ST, LP_KEYWORD_ST_N)) {
      return LpSectionKeyword::CON;
   }

   if (iskeyword(str, LP_KEYWORD_BOUNDS, LP_KEYWORD_BOUNDS_N)) {
      return LpSectionKeyword::BOUNDS;
   }

   if (iskeyword(str, LP_KEYWORD_BIN, LP_KEYWORD_BIN_N)) {
      return LpSectionKeyword::BIN;
   }

   if (iskeyword(str, LP_KEYWORD_GEN, LP_KEYWORD_GEN_N)) {
      return LpSectionKeyword::GEN;
   }

   if (iskeyword(str, LP_KEYWORD_SEMI, LP_KEYWORD_SEMI_N)) {
      return LpSectionKeyword::SEMI;
   }

   if (iskeyword(str, LP_KEYWORD_SOS, LP_KEYWORD_SOS_N)) {
      return LpSectionKeyword::SOS;
   }

   if (iskeyword(str, LP_KEYWORD_END, LP_KEYWORD_END_N)) {
      return LpSectionKeyword::END;
   }

   return LpSectionKeyword::NONE;
}

Model Reader::read() {
   tokenize();
   processtokens();
   splittokens();
   processsections();

   return builder.model;
}

void Reader::processnonesec() {
   lpassert(sectiontokens[LpSectionKeyword::NONE].empty());
}

void Reader::parseexpression(std::vector<std::unique_ptr<ProcessedToken>>& tokens, std::shared_ptr<Expression> expr, unsigned int& i) {
   if (tokens.size() - i >= 1 && tokens[0]->type == ProcessedTokenType::CONID) {
      expr->name = ((ProcessedConsIdToken*)tokens[0].get())->name;
      i++;
   }

   while (i<tokens.size()) {
      // const var
      if (tokens.size() - i >= 2
      && tokens[i]->type == ProcessedTokenType::CONST
      && tokens[i+1]->type == ProcessedTokenType::VARID) {
         std::string name = ((ProcessedVarIdToken*)tokens[i+1].get())->name;
         
         std::shared_ptr<LinTerm> linterm = std::shared_ptr<LinTerm>(new LinTerm());
         linterm->coef = ((ProcessedConstantToken*)tokens[i].get())->value;
         linterm->var = builder.getvarbyname(name);
         expr->linterms.push_back(linterm);

         i += 2;
         continue;
      }

      // const
      if (tokens.size() - i  >= 1 && tokens[i]->type == ProcessedTokenType::CONST) {
         expr->offset = ((ProcessedConstantToken*)tokens[i].get())->value;
         i++;
         continue;
      }
      
      // var
      if (tokens.size() - i  >= 1 && tokens[i]->type == ProcessedTokenType::VARID) {
         std::string name = ((ProcessedVarIdToken*)tokens[i].get())->name;
         
         std::shared_ptr<LinTerm> linterm = std::shared_ptr<LinTerm>(new LinTerm());
         linterm->coef = 1.0;
         linterm->var = builder.getvarbyname(name);
         expr->linterms.push_back(linterm);

         i++;
         continue;
      }

      // quadratic expression
      if (tokens.size() - i >= 2 && tokens[i]->type == ProcessedTokenType::BRKOP) {
         i++;
         while (i < tokens.size() && tokens[i]->type != ProcessedTokenType::BRKCL) {
            // const var hat const
            if (tokens.size() - i >= 4
            && tokens[i]->type == ProcessedTokenType::CONST
            && tokens[i+1]->type == ProcessedTokenType::VARID
            && tokens[i+2]->type == ProcessedTokenType::HAT
            && tokens[i+3]->type == ProcessedTokenType::CONST) {
               std::string name = ((ProcessedVarIdToken*)tokens[i+1].get())->name;

               lpassert (((ProcessedConstantToken*)tokens[i+3].get())->value == 2.0);

               std::shared_ptr<QuadTerm> quadterm = std::shared_ptr<QuadTerm>(new QuadTerm());
               quadterm->coef = ((ProcessedConstantToken*)tokens[i].get())->value;
               quadterm->var1 = builder.getvarbyname(name);
               quadterm->var2 = builder.getvarbyname(name);
               expr->quadterms.push_back(quadterm);

               i += 4;
               continue;
            }

            // var hat const
            if (tokens.size() - i >= 3
            && tokens[i]->type == ProcessedTokenType::VARID
            && tokens[i+1]->type == ProcessedTokenType::HAT
            && tokens[i+2]->type == ProcessedTokenType::CONST) {
               std::string name = ((ProcessedVarIdToken*)tokens[i].get())->name;

               lpassert (((ProcessedConstantToken*)tokens[i+2].get())->value == 2.0);

               std::shared_ptr<QuadTerm> quadterm = std::shared_ptr<QuadTerm>(new QuadTerm());
               quadterm->coef = 1.0;
               quadterm->var1 = builder.getvarbyname(name);
               quadterm->var2 = builder.getvarbyname(name);
               expr->quadterms.push_back(quadterm);

               i += 3;
               continue;
            }

            // const var asterisk var
            if (tokens.size() - i >= 4
            && tokens[i]->type == ProcessedTokenType::CONST
            && tokens[i+1]->type == ProcessedTokenType::VARID
            && tokens[i+2]->type == ProcessedTokenType::ASTERISK
            && tokens[i+3]->type == ProcessedTokenType::VARID) {
               std::string name1 = ((ProcessedVarIdToken*)tokens[i+1].get())->name;
               std::string name2 = ((ProcessedVarIdToken*)tokens[i+3].get())->name;

               std::shared_ptr<QuadTerm> quadterm = std::shared_ptr<QuadTerm>(new QuadTerm());
               quadterm->coef = ((ProcessedConstantToken*)tokens[i].get())->value;
               quadterm->var1 = builder.getvarbyname(name1);
               quadterm->var2 = builder.getvarbyname(name2);
               expr->quadterms.push_back(quadterm);

               i += 4;
               continue;
            }

            // var asterisk var
            if (tokens.size() - i >= 3
            && tokens[i]->type == ProcessedTokenType::VARID
            && tokens[i+1]->type == ProcessedTokenType::ASTERISK
            && tokens[i+2]->type == ProcessedTokenType::VARID) {
               std::string name1 = ((ProcessedVarIdToken*)tokens[i].get())->name;
               std::string name2 = ((ProcessedVarIdToken*)tokens[i+2].get())->name;

               std::shared_ptr<QuadTerm> quadterm = std::shared_ptr<QuadTerm>(new QuadTerm());
               quadterm->coef = 1.0;
               quadterm->var1 = builder.getvarbyname(name1);
               quadterm->var2 = builder.getvarbyname(name2);
               expr->quadterms.push_back(quadterm);

               i += 3;
               continue;
            }
         }
         lpassert(tokens.size() - i >= 3);
         lpassert(tokens[i]->type == ProcessedTokenType::BRKCL);
         lpassert(tokens[i+1]->type == ProcessedTokenType::SLASH);
         lpassert(tokens[i+2]->type == ProcessedTokenType::CONST);
         lpassert(((ProcessedConstantToken*)tokens[i+2].get())->value == 2.0);
         i += 3;
         continue;
      }

      break;
   }
}

void Reader::processobjsec() {
   builder.model.objective = std::shared_ptr<Expression>(new Expression);
   unsigned int i = 0;   
   parseexpression(sectiontokens[LpSectionKeyword::OBJ], builder.model.objective, i);
   lpassert(i == sectiontokens[LpSectionKeyword::OBJ].size());
}

void Reader::processconsec() {
   unsigned int i=0;
   while (i<sectiontokens[LpSectionKeyword::CON].size()) {
      std::shared_ptr<Constraint> con = std::shared_ptr<Constraint>(new Constraint);
      parseexpression(sectiontokens[LpSectionKeyword::CON], con->expr, i);
      lpassert(sectiontokens[LpSectionKeyword::CON].size() - i >= 2);
	  lpassert(sectiontokens[LpSectionKeyword::CON][i]->type == ProcessedTokenType::COMP);
      lpassert(sectiontokens[LpSectionKeyword::CON][i+1]->type == ProcessedTokenType::CONST);
      double value = ((ProcessedConstantToken*)sectiontokens[LpSectionKeyword::CON][i+1].get())->value;
      switch (((ProcessedComparisonToken*)sectiontokens[LpSectionKeyword::CON][i].get())->dir) {
         case LpComparisonType::EQ:
            con->lowerbound = con->upperbound = value;
            break;
         case LpComparisonType::LEQ:
            con->upperbound = value;
            break;
         case LpComparisonType::GEQ:
            con->lowerbound = value;
            break;
         default:
            lpassert(false);
      }
      i += 2;
      builder.model.constraints.push_back(con);
   }
}

void Reader::processboundssec() {
   unsigned int i=0;
   while (i<sectiontokens[LpSectionKeyword::BOUNDS].size()) {
      // VAR free
      if (sectiontokens[LpSectionKeyword::BOUNDS].size() - i >= 2
         && sectiontokens[LpSectionKeyword::BOUNDS][i]->type == ProcessedTokenType::VARID
         && sectiontokens[LpSectionKeyword::BOUNDS][i+1]->type == ProcessedTokenType::FREE) {
         std::string name = ((ProcessedVarIdToken*)sectiontokens[LpSectionKeyword::BOUNDS][i].get())->name;
         std::shared_ptr<Variable> var = builder.getvarbyname(name);
         var->lowerbound = -std::numeric_limits<double>::infinity(); 
         var->upperbound = std::numeric_limits<double>::infinity();
         i += 2;
		 continue;
      }

	  // CONST COMP VAR COMP CONST
	  if (sectiontokens[LpSectionKeyword::BOUNDS].size() - i >= 5
		  && sectiontokens[LpSectionKeyword::BOUNDS][i]->type == ProcessedTokenType::CONST
		  && sectiontokens[LpSectionKeyword::BOUNDS][i + 1]->type == ProcessedTokenType::COMP
		  && sectiontokens[LpSectionKeyword::BOUNDS][i + 2]->type == ProcessedTokenType::VARID
		  && sectiontokens[LpSectionKeyword::BOUNDS][i + 3]->type == ProcessedTokenType::COMP
		  && sectiontokens[LpSectionKeyword::BOUNDS][i + 4]->type == ProcessedTokenType::CONST) {
		  lpassert(((ProcessedComparisonToken*)sectiontokens[LpSectionKeyword::BOUNDS][i + 1].get())->dir == LpComparisonType::LEQ);
		  lpassert(((ProcessedComparisonToken*)sectiontokens[LpSectionKeyword::BOUNDS][i + 3].get())->dir == LpComparisonType::LEQ);

		  double lb = ((ProcessedConstantToken*)sectiontokens[LpSectionKeyword::BOUNDS][i].get())->value;
		  double ub = ((ProcessedConstantToken*)sectiontokens[LpSectionKeyword::BOUNDS][i + 4].get())->value;

		  std::string name = ((ProcessedVarIdToken*)sectiontokens[LpSectionKeyword::BOUNDS][i + 2].get())->name;
		  std::shared_ptr<Variable> var = builder.getvarbyname(name);

		  var->lowerbound = lb;
		  var->upperbound = ub;

		  i += 5;
		  continue;
	  }

      // CONST COMP VAR
      if (sectiontokens[LpSectionKeyword::BOUNDS].size() - i >= 3
      && sectiontokens[LpSectionKeyword::BOUNDS][i]->type == ProcessedTokenType::CONST
      && sectiontokens[LpSectionKeyword::BOUNDS][i+1]->type == ProcessedTokenType::COMP
      && sectiontokens[LpSectionKeyword::BOUNDS][i+2]->type == ProcessedTokenType::VARID) {
         double value = ((ProcessedConstantToken*)sectiontokens[LpSectionKeyword::BOUNDS][i].get())->value;
         std::string name = ((ProcessedVarIdToken*)sectiontokens[LpSectionKeyword::BOUNDS][i+2].get())->name;
         std::shared_ptr<Variable> var = builder.getvarbyname(name);
         LpComparisonType dir = ((ProcessedComparisonToken*)sectiontokens[LpSectionKeyword::BOUNDS][i+1].get())->dir;

         lpassert(dir != LpComparisonType::L && dir != LpComparisonType::G);

         switch (dir) {
            case LpComparisonType::LEQ:
               var->lowerbound = value;
               break;
            case LpComparisonType::GEQ:
               var->upperbound = value;
               break;
            case LpComparisonType::EQ:
               var->lowerbound = var->upperbound = value;
               break;
            default:
               lpassert(false);
         }
         i += 3;
         continue;
      }

      // VAR COMP CONST
      if (sectiontokens[LpSectionKeyword::BOUNDS].size() -i >= 3
      && sectiontokens[LpSectionKeyword::BOUNDS][i]->type == ProcessedTokenType::VARID
      && sectiontokens[LpSectionKeyword::BOUNDS][i+1]->type == ProcessedTokenType::COMP
      && sectiontokens[LpSectionKeyword::BOUNDS][i+2]->type == ProcessedTokenType::CONST) {
         double value = ((ProcessedConstantToken*)sectiontokens[LpSectionKeyword::BOUNDS][i+2].get())->value;
         std::string name = ((ProcessedVarIdToken*)sectiontokens[LpSectionKeyword::BOUNDS][i].get())->name;
         std::shared_ptr<Variable> var = builder.getvarbyname(name);
         LpComparisonType dir = ((ProcessedComparisonToken*)sectiontokens[LpSectionKeyword::BOUNDS][i+1].get())->dir;

         lpassert(dir != LpComparisonType::L && dir != LpComparisonType::G);

         switch (dir) {
            case LpComparisonType::LEQ:
               var->upperbound = value;
               break;
            case LpComparisonType::GEQ:
               var->lowerbound = value;
               break;
            case LpComparisonType::EQ:
               var->lowerbound = var->upperbound = value;
               break;
            default:
               lpassert(false);
         }
         i += 3;
         continue;
      }
      
	  lpassert(false);
   }
}

void Reader::processbinsec() {
   for (unsigned int i=0; i<sectiontokens[LpSectionKeyword::BIN].size(); i++) {
      lpassert(sectiontokens[LpSectionKeyword::BIN][i]->type == ProcessedTokenType::VARID);
      std::string name = ((ProcessedVarIdToken*)sectiontokens[LpSectionKeyword::BIN][i].get())->name;
      std::shared_ptr<Variable> var = builder.getvarbyname(name);
      var->type = VariableType::BINARY;
   }
}

void Reader::processgensec() {
   for (unsigned int i=0; i<sectiontokens[LpSectionKeyword::GEN].size(); i++) {
      lpassert(sectiontokens[LpSectionKeyword::GEN][i]->type == ProcessedTokenType::VARID);
      std::string name = ((ProcessedVarIdToken*)sectiontokens[LpSectionKeyword::GEN][i].get())->name;
      std::shared_ptr<Variable> var = builder.getvarbyname(name);
      var->type = VariableType::GENERAL;
   }
}

void Reader::processsemisec() {
   for (unsigned int i=0; i<sectiontokens[LpSectionKeyword::SEMI].size(); i++) {
      lpassert(sectiontokens[LpSectionKeyword::SEMI][i]->type == ProcessedTokenType::VARID);
      std::string name = ((ProcessedVarIdToken*)sectiontokens[LpSectionKeyword::SEMI][i].get())->name;
      std::shared_ptr<Variable> var = builder.getvarbyname(name);
      var->type = VariableType::SEMICONTINUOUS;
   }
}

void Reader::processsossec() {
   // TODO
   lpassert(sectiontokens[LpSectionKeyword::SOS].empty());
}

void Reader::processendsec() {
   lpassert(sectiontokens[LpSectionKeyword::END].empty());
}

void Reader::processsections() {
   processnonesec();
   processobjsec();
   processconsec();
   processboundssec();
   processgensec();
   processbinsec();
   processsemisec();
   processsossec();
   processendsec();
}

void Reader::splittokens() {
   LpSectionKeyword currentsection = LpSectionKeyword::NONE;
   
   for (unsigned int i=0; i < processedtokens.size(); ++i) {
      if (processedtokens[i]->type == ProcessedTokenType::SECID) {
         currentsection = ((ProcessedTokenSectionKeyword*)processedtokens[i].get())->keyword;
         
         if (currentsection == LpSectionKeyword::OBJ) {
            switch(((ProcessedTokenObjectiveSectionKeyword*)processedtokens[i].get())->objsense) {
               case LpObjectiveSectionKeywordType::MIN:
                  builder.model.sense = ObjectiveSense::MIN;
                  break;
               case LpObjectiveSectionKeywordType::MAX:
                  builder.model.sense = ObjectiveSense::MAX;
                  break;
               default:
                  lpassert(false);
            }
         }

         // make sure this section did not yet occur
         lpassert(sectiontokens[currentsection].empty());
      } else {
         sectiontokens[currentsection].push_back(std::move(processedtokens[i]));
      }
   }
}

void Reader::processtokens() {
   unsigned int i = 0;
   
   while (i < this->rawtokens.size()) {
      fflush(stdout);

      // long section keyword semi-continuous
      if (rawtokens.size() - i >= 3 && rawtokens[i]->istype(RawTokenType::STR) && rawtokens[i+1]->istype(RawTokenType::MINUS) && rawtokens[i+2]->istype(RawTokenType::STR)) {
         std::string temp = ((RawStringToken*)rawtokens[i].get())->value + "-" + ((RawStringToken*)rawtokens[i+2].get())->value;
         LpSectionKeyword keyword = parsesectionkeyword(temp);
         if (keyword != LpSectionKeyword::NONE) {
            processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedTokenSectionKeyword(keyword)));
            i += 3;
            continue;
         }
      }

      // long section keyword subject to/such that
      if (rawtokens.size() - i >= 2 && rawtokens[i]->istype(RawTokenType::STR) && rawtokens[i+1]->istype(RawTokenType::STR)) {
         std::string temp = ((RawStringToken*)rawtokens[i].get())->value + " " + ((RawStringToken*)rawtokens[i+1].get())->value;
         LpSectionKeyword keyword = parsesectionkeyword(temp);
         if (keyword != LpSectionKeyword::NONE) {
            processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedTokenSectionKeyword(keyword)));
            i += 2;
            continue;
         }
      }

      // other section keyword
      if (rawtokens[i]->istype(RawTokenType::STR)) {
         LpSectionKeyword keyword = parsesectionkeyword(((RawStringToken*)rawtokens[i].get())->value);
         if (keyword != LpSectionKeyword::NONE) {
            if (keyword == LpSectionKeyword::OBJ) {
               LpObjectiveSectionKeywordType kw = parseobjectivesectionkeyword(((RawStringToken*)rawtokens[i].get())->value);
               processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedTokenObjectiveSectionKeyword(kw)));
            } else {
               processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedTokenSectionKeyword(keyword)));
            }
            i++;
            continue;
         }
      }

      // constraint identifier?
      if (rawtokens.size() - i >= 2 && rawtokens[i]->istype(RawTokenType::STR) && rawtokens[i+1]->istype(RawTokenType::COLON)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedConsIdToken(((RawStringToken*)rawtokens[i].get())->value)));
         i += 2;
         continue;
      }

      // check if free
      if (rawtokens[i]->istype(RawTokenType::STR) && iskeyword(((RawStringToken*)rawtokens[i].get())->value, LP_KEYWORD_FREE, LP_KEYWORD_FREE_N)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedToken(ProcessedTokenType::FREE)));
         i++;
         continue;
      }

      // check if infinty
      if (rawtokens[i]->istype(RawTokenType::STR) && iskeyword(((RawStringToken*)rawtokens[i].get())->value, LP_KEYWORD_INF, LP_KEYWORD_INF_N)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedConstantToken(std::numeric_limits<double>::infinity())));
         i++;
         continue;
      }

      // assume var identifier
      if (rawtokens[i]->istype(RawTokenType::STR)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedVarIdToken(((RawStringToken*)rawtokens[i].get())->value)));
         i++;
         continue;
      }

      // + Constant
      if (rawtokens.size() - i >= 2 && rawtokens[i]->istype(RawTokenType::PLUS) && rawtokens[i+1]->istype(RawTokenType::CONS)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedConstantToken(((RawConstantToken*)rawtokens[i+1].get())->value)));
         i += 2;
         continue;
      }

      // - constant
      if (rawtokens.size() - i >= 2 && rawtokens[i]->istype(RawTokenType::MINUS) && rawtokens[i+1]->istype(RawTokenType::CONS)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedConstantToken(-((RawConstantToken*)rawtokens[i+1].get())->value)));
         i += 2;
         continue;
      }

      // +
      if (rawtokens[i]->istype(RawTokenType::PLUS)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedConstantToken(1.0)));
         i++;
         continue;
      }

      // -
      if (rawtokens[i]->istype(RawTokenType::MINUS)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedConstantToken(-1.0)));
         i++;
         continue;
      }

      // constant
      if (rawtokens[i]->istype(RawTokenType::CONS)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedConstantToken(((RawConstantToken*)rawtokens[i].get())->value)));
         i++;
         continue;
      }

      // [
      if (rawtokens[i]->istype(RawTokenType::BRKOP)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedToken(ProcessedTokenType::BRKOP)));
         i++;
         continue;
      }

      // ]
      if (rawtokens[i]->istype(RawTokenType::BRKCL)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedToken(ProcessedTokenType::BRKCL)));
         i++;
         continue;
      }

      // /
      if (rawtokens[i]->istype(RawTokenType::SLASH)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedToken(ProcessedTokenType::SLASH)));
         i++;
         continue;
      }

      // *
      if (rawtokens[i]->istype(RawTokenType::ASTERISK)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedToken(ProcessedTokenType::ASTERISK)));
         i++;
         continue;
      }

      // ^
      if (rawtokens[i]->istype(RawTokenType::HAT)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedToken(ProcessedTokenType::HAT)));
         i++;
         continue;
      }

      // <=
      if (rawtokens.size() - i >= 2 && rawtokens[i]->istype(RawTokenType::LESS) && rawtokens[i+1]->istype(RawTokenType::EQUAL)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedComparisonToken(LpComparisonType::LEQ)));
         i += 2;
         continue;
      }

      // <
      if (rawtokens[i]->istype(RawTokenType::LESS)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedComparisonToken(LpComparisonType::L)));
         i++;
         continue;
      }

      // >=
      if (rawtokens.size() - i >= 2 && rawtokens[i]->istype(RawTokenType::GREATER) && rawtokens[i+1]->istype(RawTokenType::EQUAL)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedComparisonToken(LpComparisonType::GEQ)));
         i += 2;
         continue;
      }

      // >
      if (rawtokens[i]->istype(RawTokenType::GREATER)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedComparisonToken(LpComparisonType::G)));
         i++;
         continue;
      }

      // =
      if (rawtokens[i]->istype(RawTokenType::EQUAL)) {
         processedtokens.push_back(std::unique_ptr<ProcessedToken>(new ProcessedComparisonToken(LpComparisonType::EQ)));
         i++;
         continue;
      }

      // FILEEND
      if (rawtokens[i]->istype(RawTokenType::FLEND)) {
         i++;
         continue;
      }

      // catch all unknown symbols
      lpassert(false);
      break;
   }
}

// reads the entire file and separates 
void Reader::tokenize() {
   this->linebufferrefill = true;
   bool done = false;
   while(true) {
      this->readnexttoken(done);
      if (this->rawtokens.size() >= 1 && this->rawtokens.back()->type == RawTokenType::FLEND) {
         break;
      }
   }
}

void Reader::readnexttoken(bool& done) {
   done = false;
   if (this->linebufferrefill) {
      char* eof = fgets(this->linebuffer, LP_MAX_LINE_LENGTH+1, this->file);
      this->linebufferpos = this->linebuffer;
      this->linebufferrefill = false;

      // fgets returns NULL if end of file reached (EOF following a \n)
      if (eof == NULL) {
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::FLEND)));
         done = true;
         return;
      }
   }

   // check single character tokens
   char nextchar = *this->linebufferpos;

   switch (nextchar) {
      // check for comment
      case '\\':
         this->linebufferrefill = true;
         return;
      
      // check for bracket opening
      case '[':
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::BRKOP)));
         this->linebufferpos++;
         return;

      // check for bracket closing
      case ']':
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::BRKCL)));
         this->linebufferpos++;
         return;

      // check for less sign
      case '<':
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::LESS)));
         this->linebufferpos++;
         return;

      // check for greater sign
      case '>':
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::GREATER)));
         this->linebufferpos++;
         return;

      // check for equal sign
      case '=':
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::EQUAL)));
         this->linebufferpos++;
         return;
      
      // check for colon
      case ':':
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::COLON)));
         this->linebufferpos++;
         return;

      // check for plus
      case '+':
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::PLUS)));
         this->linebufferpos++;
         return;

      // check for hat
      case '^':
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::HAT)));
         this->linebufferpos++;
         return;

      // check for hat
      case '/':
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::SLASH)));
         this->linebufferpos++;
         return;

      // check for asterisk
      case '*':
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::ASTERISK)));
         this->linebufferpos++;
         return;
      
      // check for minus
      case '-':
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::MINUS)));
         this->linebufferpos++;
         return;

      // check for whitespace
      case ' ':
      case '\t':
         this->linebufferpos++;
         return;

      // check for line end
      case '\n':
         this->linebufferrefill = true;
         return;

      // check for file end (EOF at end of some line)
      case '\0': 
         this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawToken(RawTokenType::FLEND)));
         done = true;
         return;
   }
   
   // check for double value
   double constant;
   int ncharconsumed;
   int nread = sscanf(this->linebufferpos, "%lf%n", &constant, &ncharconsumed);
   if (nread == 1) {
      this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawConstantToken(constant)));
      this->linebufferpos += ncharconsumed;
      return;
   }

   // assume it's an (section/variable/constraint) idenifier
   char stringbuffer[LP_MAX_NAME_LENGTH+1];
   nread = sscanf(this->linebufferpos, "%[^][\t\n\\:+<>^= /-]%n",
                 stringbuffer, &ncharconsumed);
   if (nread == 1) {
      this->rawtokens.push_back(std::unique_ptr<RawToken>(new RawStringToken(stringbuffer)));
      this->linebufferpos += ncharconsumed;
      return;
   }
   
   lpassert(false);
}
