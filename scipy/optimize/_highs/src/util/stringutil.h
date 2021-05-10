/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <ctype.h>

#include <cstring>
#include <string>

void strRemoveWhitespace(char* str);
char* strClone(const char* str);
int strIsWhitespace(const char* str);
void strToLower(char* str);
void strTrim(char* str);

std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");
std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");
std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ");

bool is_empty(std::string& str, const std::string& chars = "\t\n\v\f\r ");
bool is_empty(char c, const std::string& chars = "\t\n\v\f\r ");
bool is_end(std::string& str, int end,
            const std::string& chars = "\t\n\v\f\r ");

// todo: replace with pair of references rather than string ret value to avoid
// copy and also using function below. or do it properly with iterators.
std::string first_word(std::string& str, int start);
int first_word_end(std::string& str, int start);

#endif
