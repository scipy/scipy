/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <ctype.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <string>
#include <vector>

static const std::string kDefaultNonChars = "\t\n\v\f\r ";
static const char* kLegalIntegerCharacters = "+-0123456789eE";
static const char* kLegalDoubleCharacters = ".+-0123456789eE";

/*
void strRemoveWhitespace(char* str);
char* strClone(const char* str);
int strIsWhitespace(const char* str);
void strToLower(char* str);
void strTrim(char* str);
*/
// std::string& str_tolower(std::string s);

void tolower(std::string& str);
void toupper(std::string& str);

std::string& ltrim(std::string& str,
                   const std::string& chars = kDefaultNonChars);
std::string& rtrim(std::string& str,
                   const std::string& chars = kDefaultNonChars);
std::string& trim(std::string& str,
                  const std::string& chars = kDefaultNonChars);

template <class Iterator>
std::string joinString(Iterator begin, Iterator end,
                       const std::string& separator) {
  std::string joined;
  for (Iterator it = begin; it != end; ++it) {
    if (!joined.empty()) joined += separator;
    joined += *it;
  }
  return joined;
}

// simple class for rendering aligned text table
template <size_t NumCols>
class HighsTextTable {
 public:
  typedef std::array<std::string, NumCols> Row;

  explicit HighsTextTable(const Row& headers) : headers_(headers) {
    resetWidths();
  }

  void addRow(const Row& row) {
    rows_.push_back(row);
    for (size_t i = 0; i < NumCols; ++i) {
      widths_[i] = std::max(widths_[i], row[i].size());
    }
  }

  void clear() {
    rows_.clear();
    resetWidths();
  }

  std::string render(const std::string& separator = "  ") const {
    std::string rendered_table;
    const size_t row_capacity = getRowCapacity(separator);
    rendered_table.reserve((rows_.size() + 2) * (row_capacity + 1));

    renderRow(rendered_table, headers_, separator);
    rendered_table += '\n';

    // render separator row
    rendered_table.append(widths_[0], '-');

    for (size_t i = 1; i < NumCols; ++i) {
      rendered_table += separator;
      rendered_table.append(widths_[i], '-');
    }

    // render rows
    for (size_t row_index = 0; row_index < rows_.size(); ++row_index) {
      rendered_table += '\n';
      renderRow(rendered_table, rows_[row_index], separator);
    }

    return rendered_table;
  }

  size_t numColumns() const { return NumCols; }
  size_t numRows() const { return rows_.size(); }

  template <size_t Column>
  void sortRows() {
    static_assert(Column < NumCols, "Column out of range");
    std::sort(rows_.begin(), rows_.end(), [](const Row& lhs, const Row& rhs) {
      return lhs[Column] < rhs[Column];
    });
  }

 private:
  void resetWidths() {
    for (size_t i = 0; i < NumCols; ++i) {
      widths_[i] = headers_[i].size();
    }
  }

  size_t getRowCapacity(const std::string& separator) const {
    size_t row_capacity = 0;
    for (size_t i = 0; i < NumCols; ++i) {
      row_capacity += widths_[i];
    }
    if (NumCols > 1) row_capacity += (NumCols - 1) * separator.size();

    return row_capacity;
  }

  void renderRow(std::string& out, const Row& row,
                 const std::string& separator) const {
    for (size_t i = 0; i < NumCols; ++i) {
      if (i) out += separator;
      out += row[i];
      out.append(widths_[i] - row[i].size(), ' ');
    }
  }

  Row headers_;
  std::array<size_t, NumCols> widths_;
  std::vector<Row> rows_;
};

bool is_empty(std::string& str, const std::string& chars = kDefaultNonChars);
bool is_empty(char c, const std::string& chars = kDefaultNonChars);
bool is_end(std::string& str, size_t end,
            const std::string& chars = kDefaultNonChars);

// todo: replace with pair of references rather than string ret value to avoid
// copy and also using function below. or do it properly with iterators.
std::string first_word(std::string& str, size_t start);
size_t first_word_end(std::string& str, size_t start);

#if defined(_WIN32) || defined(_WIN64)
#include <string.h>  // _strnicmp
#else
#include <strings.h>  // strncasecmp
#endif

// check if substring is contained in str, case insensitive
inline bool is_substring_case(const char* str, const char* substring) {
  if (!str || !substring) return false;

  const size_t substring_len = strlen(substring);
  if (substring_len == 0) return true;

  const size_t str_len = strlen(str);
  if (substring_len > str_len) return false;

  const char* str_last_possible = str + str_len - substring_len;
  for (; str <= str_last_possible; ++str) {
#if defined(_WIN32) || defined(_WIN64)
    if (_strnicmp(str, substring, substring_len) == 0) return true;
#else
    if (strncasecmp(str, substring, substring_len) == 0) return true;
#endif
  }
  return false;
}

#endif
