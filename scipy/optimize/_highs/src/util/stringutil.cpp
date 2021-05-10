/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "util/stringutil.h"

void strRemoveWhitespace(char* str) {
  char* dest = str;
  do
    while (isspace(*str)) str++;
  while ((*dest++ = *str++));
}

char* strClone(const char* str) {
  size_t n;
  n = strlen(str);

  char* cpy = new char[n + 1];
  strcpy(cpy, str);
  return cpy;
}

int strIsWhitespace(const char* str) {
  while (*str != '\0') {
    if (!isspace((unsigned char)*str)) {
      return 0;
    }
    str++;
  }
  return 1;
}

void strToLower(char* str) {
  int i;
  for (i = 0; str[i] != '\0'; i++) {
    str[i] = (char)tolower(str[i]);
  }
}

void strTrim(char* str) {
  int i;
  int begin = 0;
  int end = strlen(str) - 1;

  while (isspace((unsigned char)str[begin])) begin++;

  while ((end >= begin) && isspace((unsigned char)str[end])) end--;

  // Shift all characters back to the start of the string array.
  for (i = begin; i <= end; i++) str[i - begin] = str[i];

  str[i - begin] = '\0';  // Null terminate string.
}

std::string& ltrim(std::string& str, const std::string& chars) {
  str.erase(0, str.find_first_not_of(chars));
  return str;
}

std::string& rtrim(std::string& str, const std::string& chars) {
  str.erase(str.find_last_not_of(chars) + 1);
  return str;
}

std::string& trim(std::string& str, const std::string& chars) {
  return ltrim(rtrim(str, chars), chars);
}

bool is_empty(char c, const std::string& chars) {
  int pos = chars.find_first_of(c);
  if (pos == -1 || pos == (int)chars.size()) return false;
  return true;
}

bool is_empty(std::string& str, const std::string& chars) {
  int pos = str.find_first_not_of(chars);
  if (pos == -1 || pos == (int)str.size()) return true;
  return false;
}

bool is_end(std::string& str, int end, const std::string& chars) {
  int pos = str.find_first_not_of(chars, end);
  if (pos == -1 || pos == (int)str.size()) return true;
  return false;
}

int first_word_end(std::string& str, int start) {
  const std::string chars = "\t\n\v\f\r ";
  int next_word_start = str.find_first_not_of(chars, start);
  int next_word_end = str.find_first_of(chars, next_word_start);
  if (next_word_end < 0 || next_word_end > (int)str.size()) return str.size();
  return next_word_end;
}

std::string first_word(std::string& str, int start) {
  const std::string chars = "\t\n\v\f\r ";
  int next_word_start = str.find_first_not_of(chars, start);
  int next_word_end = str.find_first_of(chars, next_word_start);
  return str.substr(next_word_start, next_word_end - next_word_start);
}
