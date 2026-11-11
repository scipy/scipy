/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/LpReader.cpp
 * @brief Reader for the CPLEX LP file format
 *
 * The reader is a hand-written recursive-descent parser for the following
 * grammar, which is based on the CPLEX LP format and the LP reader in
 * MathOptInterface.jl. Keywords are case-insensitive. `|` separates
 * alternatives, `[x]` means x is optional and `{x}` means x is repeated zero
 * or more times.
 *
 *   <file>          := { <section> } "end" | <empty>
 *   <section>       := <objective> | <constraints> | <bounds> | <integers>
 *                    | <binaries> | <semis>
 *   <objective>     := ("minimize" | "maximize") [<name>] <expression>
 *   <constraints>   := "subject to" { <constraint> }
 *   <bounds>        := "bounds" { <bound> }
 *   <integers>      := "general" { <identifier> }
 *   <binaries>      := "binary" { <identifier> }
 *   <semis>         := "semi-continuous" { <identifier> }
 *
 *   <constraint>    := [<name>] <expression> <comparison> <number> <newline>
 *   <bound>         := <identifier> "free" <newline>
 *                    | <identifier> <comparison> <number> <newline>
 *                    | <number> <comparison> <identifier>
 *                      [<comparison> <number>] <newline>
 *
 *   <name>          := (<identifier> | <number>) ":"
 *   <expression>    := [<signs>] <term> { <signs> <term> }
 *   <term>          := <constant> [["*"] <identifier>] | <identifier>
 *                    | <quadratic>
 *   <quadratic>     := "[" [[<signs>] <quad-term> { <signs> <quad-term> }] "]"
 *                      ["/" "2"]
 *   <quad-term>     := [<constant> ["*"]] <identifier>
 *                      ("^" "2" | "*" <identifier>)
 *   <number>        := [<signs>] (<constant> | "inf" | "infinity")
 *   <signs>         := ("+" | "-") { "+" | "-" }
 *   <comparison>    := "<" | "<=" | "=<" | ">" | ">=" | "=>" | "=" | "=="
 *
 * Some details are handled by the lexer and not the grammar:
 *
 *  - Everything from a backslash to the end of the line is a comment.
 *  - An identifier is only a section keyword if it is the first token on its
 *    line and it isn't followed by ":", so `end: x >= 1` is a constraint
 *    named "end". The keywords that are more than one token long are
 *    "subject to", "such that" and "semi-continuous", and all of their
 *    tokens must be on the same line.
 *  - Identifiers can contain any printable character except whitespace and
 *    the operators \ : + - * / ^ < > = [ ]. They cannot start with a digit
 *    or "/", or with "." followed by a digit. This is more permissive than
 *    CPLEX, so that names like "...100" written by HiGHS can be read.
 *  - <newline> means the next token must be on a new line, or be the end of
 *    the file.
 *  - The quadratic part of the objective is 0.5 x'Qx. Without the "/ 2", the
 *    bracketed terms are doubled.
 *
 * Some things are syntactically valid, but are not supported by HiGHS, and
 * are reported as errors: quadratic constraints, SOS constraints, and
 * indicator constraints.
 *
 * For backwards compatibility, a constant on the left-hand side of a
 * constraint is ignored with a warning, rather than being moved to the
 * right-hand side.
 *
 * The file is read in chunks, and the model is built as the file is parsed,
 * so the whole file is never held in memory. Diagnostics record the line
 * that they refer to, and the line is read from the file again when the
 * diagnostic is printed.
 */
#include "io/LpReader.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <unordered_set>
#include <utility>
#include <vector>

#include "HConfig.h"  // for ZLIB_FOUND
#ifdef ZLIB_FOUND
#include "../extern/zstr/zstr.hpp"
#endif

namespace lp_reader {

enum class TokenKind {
  kIdentifier,
  kNumber,
  kPlus,
  kMinus,
  kTimes,
  kDivide,
  kPower,
  kOpenBracket,
  kCloseBracket,
  kColon,
  kLess,
  kGreater,
  kEqual,
  kImplies,
  kEndOfFile
};

struct Token {
  TokenKind kind = TokenKind::kEndOfFile;
  // Location of the token as a byte offset into the file, its length in
  // bytes, and its line number
  size_t pos = 0;
  size_t len = 0;
  size_t line = 1;
  // Whether this is the first token on its line
  bool line_start = false;
  // The value of a kNumber token
  double value = 0;
};

// An error or warning about a location in the file
struct Diagnostic {
  std::string message;
  size_t pos;
  size_t len;
  size_t line;
  std::string label;
  std::string help;
};

// Thrown when the file is not valid
struct ParseError {
  Diagnostic diagnostic;
};

// Only this many warnings are printed
const size_t kMaxWarnings = 10;

// The size of the buffer in zstr::ifstream. 1 << 16 = 64kB.
const size_t kZlibBufferSize = 1 << 16;

// The size of the initial buffer (in bytes) when reading the lexer.
// 1 << 14 = 16kB. If the file is larger than this, the buffer is expanded, up
// to a maximum of kMaximumLexerBufferSize.
const size_t kInitialLexerBufferSize = 1 << 14;

// The maximum number of bytes that the lexer reads at a time. 1 << 20 = 1MiB.
const size_t kMaximumLexerBufferSize = 1 << 20;

// The maximum number of tokens that the parser can look ahead. The parser
// never peeks more than four tokens ahead, so this is generous. It must be a
// power of two.
const size_t kMaxLookahead = 8;

// Open a file that may be compressed. If the file cannot be opened, the
// stream is in a failed state.
static std::unique_ptr<std::istream> openFile(const std::string& filename) {
#ifdef ZLIB_FOUND
  try {
    return std::unique_ptr<std::istream>(
        new zstr::ifstream(filename, std::ios::in, kZlibBufferSize));
  } catch (const strict_fstream::Exception&) {
    std::unique_ptr<std::istream> failed(new std::istringstream());
    failed->setstate(std::ios::failbit);
    return failed;
  }
#else
  return std::unique_ptr<std::istream>(
      new std::ifstream(filename, std::ios::in | std::ios::binary));
#endif
}

static bool isUtf8Continuation(const char c) {
  return (static_cast<unsigned char>(c) & 0xC0) == 0x80;
}

// The number of UTF-8 characters in text[begin, end)
static size_t countChars(const std::string& text, size_t begin, size_t end) {
  size_t count = 0;
  for (size_t i = begin; i < end; i++)
    if (!isUtf8Continuation(text[i])) count++;
  return count;
}

// Render a diagnostic in the style of rustc:
//
//   <message>
//    --> <filename>:<line>:<column>
//     |
//   2 |  obj: x y
//     |         ^ <label>
//     |
//     = help: <help>
//
// `line` is the text of the line that the diagnostic refers to, and
// `line_pos` is the location of the start of the line in the file.
static std::string renderDiagnostic(const std::string& filename,
                                    const std::string& line, size_t line_pos,
                                    const Diagnostic& d) {
  const size_t line_end = line.size();
  const size_t pos = std::min(d.pos - line_pos, line_end);
  const size_t column = 1 + countChars(line, 0, pos);
  // Show at most kContext bytes either side of the token
  const size_t kContext = 40;
  size_t caret_end = std::min(pos + std::min(d.len, kContext), line_end);
  while (caret_end < line_end && isUtf8Continuation(line[caret_end]))
    caret_end++;
  size_t begin = 0;
  if (pos > kContext) {
    begin = pos - kContext;
    while (begin < pos && isUtf8Continuation(line[begin])) begin++;
  }
  size_t end = line_end;
  if (line_end - caret_end > kContext) {
    end = caret_end + kContext;
    while (end > caret_end && isUtf8Continuation(line[end])) end--;
  }
  std::string excerpt = line.substr(begin, end - begin);
  for (char& c : excerpt)
    if (static_cast<unsigned char>(c) < 0x20) c = ' ';
  size_t caret_offset = countChars(line, begin, pos);
  if (begin > 0) {
    excerpt = "..." + excerpt;
    caret_offset += 3;
  }
  if (end < line_end) excerpt += "...";
  const size_t num_caret =
      std::max(size_t{1}, countChars(line, pos, caret_end));
  const std::string line_str = std::to_string(d.line);
  const std::string gutter(line_str.size(), ' ');
  std::string out = d.message + "\n";
  out += gutter + "--> " + filename + ":" + line_str + ":" +
         std::to_string(column) + "\n";
  out += gutter + " |\n";
  out += line_str + " | " + excerpt + "\n";
  out += gutter + " | " + std::string(caret_offset, ' ') +
         std::string(num_caret, '^') + " " + d.label + "\n";
  if (!d.help.empty()) {
    out += gutter + " |\n";
    out += gutter + " = help: " + d.help + "\n";
  }
  return out;
}

// Render diagnostics that are sorted by their location, reading the lines
// that they refer to from the file
static std::vector<std::string> renderDiagnostics(
    const std::string& filename, const std::vector<Diagnostic>& diagnostics) {
  std::vector<std::string> rendered;
  if (diagnostics.empty()) return rendered;
  std::unique_ptr<std::istream> file = openFile(filename);
  std::string line;
  size_t line_number = 0;
  size_t line_pos = 0;
  size_t next_line_pos = 0;
  for (const Diagnostic& d : diagnostics) {
    while (line_number < d.line) {
      line.clear();
      try {
        std::getline(*file, line);
      } catch (const std::ios_base::failure&) {
        // For example, zstr throws if a compressed file is corrupt. The
        // diagnostic is rendered with the part of the line that was read.
      }
      line_pos = next_line_pos;
      next_line_pos += line.size() + 1;
      line_number++;
    }
    std::string text = line;
    text.erase(text.find_last_not_of('\r') + 1);
    rendered.push_back(renderDiagnostic(filename, text, line_pos, d));
  }
  return rendered;
}

static bool isDigit(const int c) { return c >= '0' && c <= '9'; }

static bool isIdentifierChar(const int c) {
  return c > 0x20 && c != 0x7F && std::strchr("\\:+-*^<>=[]", c) == nullptr;
}

// Lexer::lex() checks for a number before an identifier, so this doesn't need
// to exclude digits
static bool isIdentifierStart(const int c) {
  return isIdentifierChar(c) && c != '/';
}

class Lexer {
 public:
  explicit Lexer(std::istream& in) : in_(in) {}

  // Look at the token n places ahead, without consuming it. References stay
  // valid until the token is consumed.
  const Token& peek(size_t n = 0) {
    assert(n < kMaxLookahead);
    while (num_tokens_ <= n) {
      tokens_[(first_token_ + num_tokens_) % kMaxLookahead] = lex();
      num_tokens_++;
    }
    return tokens_[(first_token_ + n) % kMaxLookahead];
  }

  Token next() {
    Token token = peek();
    first_token_ = (first_token_ + 1) % kMaxLookahead;
    num_tokens_--;
    // The text of the most recent token stays available
    keep_ = token.pos;
    return token;
  }

  // The text of a token, which must be the most recent token returned by
  // next() or a token that has been peeked at. The pointer is valid until the
  // next call to peek() or next().
  const char* data(const Token& token) const {
    return buf_.data() + (token.pos - buf_start_);
  }

  std::string text(const Token& token) const {
    return std::string(data(token), token.len);
  }

 private:
  // The byte at `pos`, or -1 at the end of the file
  int at(size_t pos) {
    const size_t offset = pos - buf_start_;
    if (offset < buf_.size()) return static_cast<unsigned char>(buf_[offset]);
    return fill(pos);
  }
  int fill(size_t pos);
  Token lex();
  void lexNumber(Token& token);

  std::istream& in_;
  // The part of the file from buf_start_ that has been read and may still be
  // needed
  std::string buf_;
  size_t buf_start_ = 0;
  bool eof_ = false;
  size_t chunk_ = kInitialLexerBufferSize;
  // The bytes before keep_ are no longer needed
  size_t keep_ = 0;
  // The location of the lexer
  size_t pos_ = 0;
  size_t line_ = 1;
  bool line_start_ = true;
  // The end of the most recent token, where the end of the file is reported
  size_t last_end_ = 0;
  size_t last_line_ = 1;
  // The tokens that have been lexed but not yet consumed, as a ring buffer of
  // num_tokens_ tokens starting at first_token_
  Token tokens_[kMaxLookahead];
  size_t first_token_ = 0;
  size_t num_tokens_ = 0;
};

// Read the file until it contains `pos`, discarding bytes that are no longer
// needed. Returns the byte at `pos`, or -1 at the end of the file.
int Lexer::fill(size_t pos) {
  buf_.erase(0, keep_ - buf_start_);
  buf_start_ = keep_;
  // While we have not read the byte in position `pos`...
  while (pos - buf_start_ >= buf_.size()) {
    // The byte at position `pos` is past the end of the file.
    if (eof_) return -1;
    // Read in a new chunk to the buffer, and count how many bytes we read.
    const size_t size = buf_.size();
    buf_.resize(size + chunk_);
    in_.read(&buf_[size], chunk_);
    const size_t got = static_cast<size_t>(in_.gcount());
    if (got < chunk_) {
      // We hit the end of the file before we filled the chunk.
      buf_.resize(size + got);
      eof_ = true;
    } else if (chunk_ < kMaximumLexerBufferSize) {
      // Double the size of the chunk we read in the next round.
      chunk_ = std::min(2 * chunk_, kMaximumLexerBufferSize);
    }
  }
  return static_cast<unsigned char>(buf_[pos - buf_start_]);
}

Token Lexer::lex() {
  // Skip whitespace and comments
  int c = at(pos_);
  while (c >= 0) {
    if (c == '\\') {
      while (c >= 0 && c != '\n') c = at(++pos_);
      continue;
    }
    if (c == '\n') {
      line_++;
      line_start_ = true;
    } else if (c == ' ' || c == '\t') {
      // Performance optimization: special case the most common ASCII whitespace
      // to avoid calling std::isspace(c)
    } else if (0x21 <= c && c <= 0x7E) {
      // Performance optimization: special case ASCII characters that are not
      // whitespace to avoid calling std::isspace(c)
      break;
    } else if (!std::isspace(c)) {
      break;
    }
    c = at(++pos_);
  }
  Token token;
  if (c < 0) {
    // kEndOfFile
    token.pos = last_end_;
    token.line = last_line_;
    return token;
  }
  token.pos = pos_;
  token.line = line_;
  token.line_start = line_start_;
  line_start_ = false;
  token.len = 1;
  const int next = at(pos_ + 1);
  if (isDigit(c) || (c == '.' && isDigit(next))) {
    lexNumber(token);
  } else if (isIdentifierStart(c)) {
    token.kind = TokenKind::kIdentifier;
    while (isIdentifierChar(at(pos_ + token.len))) token.len++;
  } else if (c == '+') {
    token.kind = TokenKind::kPlus;
  } else if (c == '-') {
    token.kind = TokenKind::kMinus;
    if (next == '>') {
      token.kind = TokenKind::kImplies;
      token.len = 2;
    }
  } else if (c == '*') {
    token.kind = TokenKind::kTimes;
  } else if (c == '/') {
    token.kind = TokenKind::kDivide;
  } else if (c == '^') {
    token.kind = TokenKind::kPower;
  } else if (c == '[') {
    token.kind = TokenKind::kOpenBracket;
  } else if (c == ']') {
    token.kind = TokenKind::kCloseBracket;
  } else if (c == ':') {
    token.kind = TokenKind::kColon;
  } else if (c == '<' || c == '>') {
    // <, <=, >, >=
    token.kind = c == '<' ? TokenKind::kLess : TokenKind::kGreater;
    if (next == '=') token.len = 2;
  } else if (c == '=') {
    // =, ==, =<, =>
    token.kind = TokenKind::kEqual;
    if (next == '<') token.kind = TokenKind::kLess;
    if (next == '>') token.kind = TokenKind::kGreater;
    if (next == '=' || next == '<' || next == '>') token.len = 2;
  } else {
    // Only control characters can get here
    char shown[16];
    snprintf(shown, sizeof(shown), "0x%02X", c);
    throw ParseError{{std::string("unexpected control character ") + shown,
                      pos_, 1, line_,
                      "this character is not valid in an LP file", ""}};
  }
  pos_ += token.len;
  last_end_ = pos_;
  last_line_ = line_;
  return token;
}

// <constant> := <digits> ["." [<digits>]] [<exponent>]
//             | "." <digits> [<exponent>]
// <exponent> := ("e" | "E") ["+" | "-"] <digits>
void Lexer::lexNumber(Token& token) {
  size_t end = pos_;
  while (isDigit(at(end))) end++;
  if (at(end) == '.') {
    end++;
    while (isDigit(at(end))) end++;
  }
  const int e = at(end);
  if (e == 'e' || e == 'E') {
    size_t exponent = end + 1;
    const int sign = at(exponent);
    if (sign == '+' || sign == '-') exponent++;
    if (isDigit(at(exponent))) {
      end = exponent;
      while (isDigit(at(end))) end++;
    }
  }
  token.kind = TokenKind::kNumber;
  token.len = end - pos_;
  token.value = std::strtod(text(token).c_str(), nullptr);
}

enum class Section {
  kNone,
  kMinimize,
  kMaximize,
  kConstraints,
  kBounds,
  kGeneral,
  kBinary,
  kSemi,
  kSos,
  kEnd
};

// Parses an LP file, building the model as it goes
class Parser {
 public:
  Parser(std::istream& in, HighsModel& model)
      : lexer_(in),
        lp_(model.lp_),
        hessian_(model.hessian_),
        col_index_(16, NameHash{this}, NameEqual{this}) {
    lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    lp_.a_matrix_.start_.assign(1, 0);
  }

  // Parse the file. Throws ParseError if the file is not valid, in which
  // case the model is incomplete.
  void parse();

  // Complete the model after parse() has succeeded
  void finish(const HighsLogOptions& log_options);

  // The first kMaxWarnings warnings, in order of their location
  const std::vector<Diagnostic>& warnings() const { return warnings_; }
  size_t numWarnings() const { return num_warnings_; }

 private:
  // An entry of the current row
  struct RowEntry {
    HighsInt col;
    double value;
    bool repeated;
  };

  // An entry of the Hessian, which may be repeated
  struct HessianEntry {
    HighsInt col;
    HighsInt row;
    double value;
    bool operator<(const HessianEntry& other) const {
      return col < other.col || (col == other.col && row < other.row);
    }
  };

  // col_index_ stores column indices, and hashes and compares them by the
  // column names in lp_, so the names are not stored twice. The index -1
  // refers to the name that is being looked up.
  struct NameHash {
    const Parser* parser;
    size_t operator()(HighsInt col) const {
      return std::hash<std::string>()(parser->columnName(col));
    }
  };
  struct NameEqual {
    const Parser* parser;
    bool operator()(HighsInt col1, HighsInt col2) const {
      return parser->columnName(col1) == parser->columnName(col2);
    }
  };
  const std::string& columnName(HighsInt col) const {
    return col < 0 ? query_ : lp_.col_names_[col];
  }

  void parseObjective();
  void parseConstraint();
  void parseBound();
  void parseVariableList(Section section);
  void parseExpression(bool objective);
  void parseTerm(double sign, bool objective, const Token& start);
  void parseQuadratic(double sign);
  void parseQuadraticTerm(double sign, double scale);
  bool parseName(std::string& name);
  double parseNumber(const std::string& label);
  HighsInt parseVariable(const std::string& label, Token& var);
  HighsInt parseVariable(const std::string& label) {
    Token var;
    return parseVariable(label, var);
  }
  void addVariableTerm(HighsInt col, const Token& var, double coef,
                       bool objective);
  void addHessianEntry(HighsInt col1, HighsInt col2, double coef);
  void finishRow(std::string& name, double lower, double upper);
  void checkNewLine(const std::string& label, const std::string& help);
  void checkUnsupportedConstraint();

  Section keywordAt(size_t i, size_t& num_tokens);
  bool atKeyword() {
    size_t num_tokens;
    return keywordAt(0, num_tokens) != Section::kNone;
  }
  bool atSectionEnd() {
    return lexer_.peek().kind == TokenKind::kEndOfFile || atKeyword();
  }
  bool equalsIgnoreCase(const Token& token, const char* lower) const;
  std::string describe(const Token& token) const {
    if (token.kind == TokenKind::kEndOfFile) return "end of file";
    return quote(lexer_.data(token), token.len);
  }
  static std::string quote(const char* data, size_t len);
  std::string keywordHelp(const Token& token);
  bool isKeywordWord(const Token& token) const;
  [[noreturn]] void error(const Token& token, const std::string& message,
                          const std::string& label,
                          const std::string& help = "") const {
    throw ParseError{{message, token.pos, token.len, token.line, label, help}};
  }
  void warning(const Token& token, const std::string& message,
               const std::string& label, const std::string& help) {
    if (num_warnings_++ < kMaxWarnings)
      warnings_.push_back(
          {message, token.pos, token.len, token.line, label, help});
  }
  HighsInt getColumn(const Token& token);

  Lexer lexer_;
  HighsLp& lp_;
  HighsHessian& hessian_;
  std::vector<Diagnostic> warnings_;
  size_t num_warnings_ = 0;

  // The name that is being looked up in col_index_
  std::string query_;
  std::unordered_set<HighsInt, NameHash, NameEqual> col_index_;

  bool has_integrality_ = false;
  std::vector<HighsInt> binary_cols_;
  std::vector<HessianEntry> hessian_entries_;

  // Whether a row has been given a name beginning kRowPrefix, and whether a
  // row has a name beginning kRowPrefix in the file
  bool used_prefix_ = false;
  bool prefix_ok_ = true;

  // The entries of the current row
  std::vector<RowEntry> row_entries_;
  bool row_has_constant_ = false;
  // For each column, its position in row_entries_, or -1
  std::vector<HighsInt> col_slot_;
};

const std::string kRowPrefix = "HiGHS_R";

bool Parser::equalsIgnoreCase(const Token& token, const char* lower) const {
  if (token.len != std::strlen(lower)) return false;
  const char* data = lexer_.data(token);
  for (size_t i = 0; i < token.len; i++)
    if (std::tolower(static_cast<unsigned char>(data[i])) != lower[i])
      return false;
  return true;
}

// The section keywords that are one token long
static const std::pair<const char*, Section> kKeywords[] = {
    {"minimize", Section::kMinimize},
    {"minimise", Section::kMinimize},
    {"minimum", Section::kMinimize},
    {"min", Section::kMinimize},
    {"maximize", Section::kMaximize},
    {"maximise", Section::kMaximize},
    {"maximum", Section::kMaximize},
    {"max", Section::kMaximize},
    {"st", Section::kConstraints},
    {"s.t.", Section::kConstraints},
    {"st.", Section::kConstraints},
    {"bounds", Section::kBounds},
    {"bound", Section::kBounds},
    {"general", Section::kGeneral},
    {"generals", Section::kGeneral},
    {"gen", Section::kGeneral},
    {"integer", Section::kGeneral},
    {"integers", Section::kGeneral},
    {"binary", Section::kBinary},
    {"binaries", Section::kBinary},
    {"bin", Section::kBinary},
    {"semi", Section::kSemi},
    {"semis", Section::kSemi},
    {"sos", Section::kSos},
    {"end", Section::kEnd}};

static bool isComparison(const Token& token) {
  return token.kind == TokenKind::kLess || token.kind == TokenKind::kGreater ||
         token.kind == TokenKind::kEqual;
}

static bool isSign(const Token& token) {
  return token.kind == TokenKind::kPlus || token.kind == TokenKind::kMinus;
}

// If the tokens from position i are a section keyword, return the section and
// set num_tokens to the number of tokens in the keyword
Section Parser::keywordAt(size_t i, size_t& num_tokens) {
  const Token& token = lexer_.peek(i);
  if (token.kind != TokenKind::kIdentifier || !token.line_start) {
    return Section::kNone;
  }
  // A keyword followed by ":" is a name
  const Token& second = lexer_.peek(i + 1);
  if (second.kind == TokenKind::kColon) return Section::kNone;
  num_tokens = 2;
  if ((equalsIgnoreCase(token, "subject") && equalsIgnoreCase(second, "to")) ||
      (equalsIgnoreCase(token, "such") && equalsIgnoreCase(second, "that"))) {
    return second.line_start ? Section::kNone : Section::kConstraints;
  }
  num_tokens = 3;
  // "semi-continuous", possibly with spaces, but on one line
  if (equalsIgnoreCase(token, "semi") && second.kind == TokenKind::kMinus &&
      !second.line_start) {
    const Token& third = lexer_.peek(i + 2);
    if (third.kind == TokenKind::kIdentifier && !third.line_start &&
        equalsIgnoreCase(third, "continuous"))
      return Section::kSemi;
  }
  num_tokens = 1;
  for (const auto& keyword : kKeywords)
    if (equalsIgnoreCase(token, keyword.first)) return keyword.second;
  return Section::kNone;
}

std::string Parser::quote(const char* data, size_t len) {
  const size_t kMaxLength = 30;
  if (len <= kMaxLength) return "`" + std::string(data, len) + "`";
  len = kMaxLength;
  // Back up to the start of a UTF-8 character, unless there isn't one (the
  // token is invalid UTF-8)
  while (len > 0 && isUtf8Continuation(data[len])) len--;
  if (len == 0) len = kMaxLength;
  return "`" + std::string(data, len) + "...`";
}

// If the token would be a section keyword at the start of a line, explain
// that keywords must start a line
std::string Parser::keywordHelp(const Token& token) {
  if (token.kind != TokenKind::kIdentifier || token.line_start ||
      !isKeywordWord(token))
    return "";
  return "section keywords like " + describe(token) +
         " must be at the start of a line";
}

// Whether the token is a section keyword, or the first word of one
bool Parser::isKeywordWord(const Token& token) const {
  if (equalsIgnoreCase(token, "subject") || equalsIgnoreCase(token, "such"))
    return true;
  for (const auto& keyword : kKeywords)
    if (equalsIgnoreCase(token, keyword.first)) return true;
  return false;
}

// The index of the column named by the token, adding the column if it is new
HighsInt Parser::getColumn(const Token& token) {
  query_.assign(lexer_.data(token), token.len);
  auto it = col_index_.find(-1);
  if (it != col_index_.end()) return *it;
  const HighsInt col = static_cast<HighsInt>(lp_.col_names_.size());
  lp_.col_names_.push_back(query_);
  col_index_.insert(col);
  lp_.col_cost_.push_back(0);
  lp_.col_lower_.push_back(0);
  lp_.col_upper_.push_back(kHighsInf);
  lp_.integrality_.push_back(HighsVarType::kContinuous);
  col_slot_.push_back(-1);
  return col;
}

void Parser::parse() {
  // An empty file (or one with only comments) is an empty model
  if (lexer_.peek().kind == TokenKind::kEndOfFile) return;
  bool have_objective = false;
  size_t objective_line = 0;
  while (true) {
    const Token token = lexer_.peek();
    size_t num_tokens = 0;
    const Section section = keywordAt(0, num_tokens);
    if (section == Section::kNone) {
      // Each section consumes everything up to the next keyword, so we can
      // only get here at the start or the end of the file.
      if (token.kind == TokenKind::kEndOfFile)
        error(token, "expected `end`, found end of file",
              "expected `end` after this",
              "an LP file must finish with the keyword `end`");
      error(token, "expected a section keyword, found " + describe(token),
            "expected a keyword like `minimize` or `maximize`",
            "an LP file must start with a section keyword");
    }
    for (size_t i = 0; i < num_tokens; i++) lexer_.next();
    switch (section) {
      case Section::kMinimize:
      case Section::kMaximize:
        if (have_objective)
          error(token, "the objective is defined more than once",
                "second objective section",
                "the first objective section is on line " +
                    std::to_string(objective_line));
        have_objective = true;
        objective_line = token.line;
        lp_.sense_ = section == Section::kMinimize ? ObjSense::kMinimize
                                                   : ObjSense::kMaximize;
        parseObjective();
        break;
      case Section::kConstraints:
        while (!atSectionEnd()) parseConstraint();
        break;
      case Section::kBounds:
        while (!atSectionEnd()) parseBound();
        break;
      case Section::kGeneral:
      case Section::kBinary:
      case Section::kSemi:
        parseVariableList(section);
        break;
      case Section::kSos:
        error(token, "SOS constraints are not supported by HiGHS",
              "SOS section");
      default: {
        // Section::kEnd
        const Token& after = lexer_.peek();
        if (after.kind != TokenKind::kEndOfFile)
          error(after, "expected end of file, found " + describe(after),
                "unexpected content after `end`");
        return;
      }
    }
  }
}

// <name> := (<identifier> | <number>) ":"
bool Parser::parseName(std::string& name) {
  const TokenKind kind = lexer_.peek(0).kind;
  if ((kind != TokenKind::kIdentifier && kind != TokenKind::kNumber) ||
      lexer_.peek(1).kind != TokenKind::kColon)
    return false;
  name = lexer_.text(lexer_.next());
  lexer_.next();  // :
  return true;
}

// <number> := [<signs>] (<constant> | "inf" | "infinity")
double Parser::parseNumber(const std::string& label) {
  double sign = 1;
  while (isSign(lexer_.peek()))
    if (lexer_.next().kind == TokenKind::kMinus) sign = -sign;
  const Token token = lexer_.peek();
  if (token.kind == TokenKind::kNumber) {
    lexer_.next();
    return sign * token.value;
  }
  if (token.kind == TokenKind::kIdentifier &&
      (equalsIgnoreCase(token, "inf") || equalsIgnoreCase(token, "infinity"))) {
    lexer_.next();
    return sign * kHighsInf;
  }
  error(token, "expected a number, found " + describe(token), label,
        keywordHelp(token));
}

HighsInt Parser::parseVariable(const std::string& label, Token& var) {
  var = lexer_.peek();
  if (var.kind != TokenKind::kIdentifier || atKeyword())
    error(var, "expected a variable, found " + describe(var), label);
  return getColumn(lexer_.next());
}

// The next token must be on a new line
void Parser::checkNewLine(const std::string& label, const std::string& help) {
  const Token& token = lexer_.peek();
  if (token.kind != TokenKind::kEndOfFile && !token.line_start)
    error(token, "expected a new line, found " + describe(token), label, help);
}

// <objective> := ("minimize" | "maximize") [<name>] <expression>
void Parser::parseObjective() {
  parseName(lp_.objective_name_);
  parseExpression(true);
}

// <expression> := [<signs>] <term> { <signs> <term> }
//
// In the objective, the expression ends at the next section. In a
// constraint, it ends at the comparison operator.
void Parser::parseExpression(bool objective) {
  bool first = true;
  while (true) {
    const Token token = lexer_.peek();
    if (token.kind == TokenKind::kEndOfFile || atKeyword() ||
        (!objective && isComparison(token)))
      return;
    bool has_sign = false;
    double sign = 1;
    while (isSign(lexer_.peek())) {
      if (lexer_.next().kind == TokenKind::kMinus) sign = -sign;
      has_sign = true;
    }
    if (!first && !has_sign) {
      std::string help = keywordHelp(token);
      const Token& next = lexer_.peek(1);
      if (token.kind == TokenKind::kNumber &&
          next.kind == TokenKind::kIdentifier &&
          next.pos == token.pos + token.len) {
        Token name = token;
        name.len += next.len;
        help = describe(name) +
               " looks like a name, but names cannot start with a digit";
      }
      if (help.empty())
        help = "terms in an expression must be separated by `+` or `-`";
      error(token, "expected `+` or `-`, found " + describe(token),
            "expected `+` or `-` before this", help);
    }
    parseTerm(sign, objective, token);
    first = false;
  }
}

// <term> := <constant> [["*"] <identifier>] | <identifier> | <quadratic>
//
// `start` is the first token of the term, including any signs
void Parser::parseTerm(double sign, bool objective, const Token& start) {
  const Token token = lexer_.peek();
  const bool has_sign = start.pos != token.pos;
  if (token.kind == TokenKind::kNumber) {
    lexer_.next();
    const double coef = sign * token.value;
    const Token& next = lexer_.peek();
    if (next.kind == TokenKind::kTimes) {
      lexer_.next();
      Token var;
      const HighsInt col = parseVariable("expected a variable after `*`", var);
      addVariableTerm(col, var, coef, objective);
    } else if (next.kind == TokenKind::kIdentifier && !atKeyword()) {
      const Token var = lexer_.next();
      addVariableTerm(getColumn(var), var, coef, objective);
    } else if (next.kind == TokenKind::kOpenBracket) {
      error(next, "expected `+` or `-`, found `[`",
            "expected `+` or `-` before this",
            "a coefficient cannot multiply `[ ]`; move it inside the brackets");
    } else if (objective) {
      lp_.offset_ += coef;
    } else if (coef != 0 && !row_has_constant_) {
      // Ignore the constant, for backwards compatibility
      row_has_constant_ = true;
      Token constant = start;
      constant.len = token.pos + token.len - start.pos;
      warning(constant,
              "a constant on the left-hand side of a constraint is ignored",
              "this constant is ignored",
              "move the constant to the right-hand side");
    }
  } else if (token.kind == TokenKind::kIdentifier && !atKeyword()) {
    const Token var = lexer_.next();
    addVariableTerm(getColumn(var), var, sign, objective);
  } else if (token.kind == TokenKind::kOpenBracket) {
    if (!objective)
      error(token, "quadratic constraints are not supported by HiGHS",
            "quadratic term in a constraint");
    parseQuadratic(sign);
  } else {
    error(token, "expected a term, found " + describe(token),
          has_sign ? "expected a number, a variable, or `[` after the sign"
                   : "expected a number, a variable, or `[`",
          token.kind == TokenKind::kIdentifier
              ? describe(token) +
                    " is a section keyword because it is at the start of a "
                    "line"
              : "");
  }
}

void Parser::addVariableTerm(HighsInt col, const Token& var, double coef,
                             bool objective) {
  if (objective) {
    lp_.col_cost_[col] += coef;
  } else {
    // Sum repeated variables, warning once for each repeated variable
    HighsInt& slot = col_slot_[col];
    if (slot < 0) {
      slot = static_cast<HighsInt>(row_entries_.size());
      row_entries_.push_back({col, coef, false});
    } else {
      RowEntry& entry = row_entries_[slot];
      if (!entry.repeated)
        warning(
            var,
            "variable " +
                quote(lp_.col_names_[col].data(), lp_.col_names_[col].size()) +
                " appears more than once in this constraint",
            "repeated here", "the coefficients are summed");
      entry.repeated = true;
      entry.value += coef;
    }
  }
  const Token& next = lexer_.peek();
  if (next.kind == TokenKind::kTimes || next.kind == TokenKind::kPower)
    error(next, "quadratic terms must be inside `[` and `]`",
          "unexpected " + describe(next) + " outside `[ ]`");
}

// <quadratic> := "[" [[<signs>] <quad-term> { <signs> <quad-term> }] "]"
//                ["/" "2"]
void Parser::parseQuadratic(double sign) {
  const Token open = lexer_.next();
  // The bracketed terms are q(x). With "/ 2" they are q(x) / 2 = 0.5 x'Qx,
  // so x'Qx = q(x), otherwise x'Qx = 2 q(x). We don't know which until the
  // end, so add the terms of 2 q(x), and halve them if there is a "/ 2".
  const size_t first_entry = hessian_entries_.size();
  bool first = true;
  while (lexer_.peek().kind != TokenKind::kCloseBracket) {
    const Token token = lexer_.peek();
    if (token.kind == TokenKind::kEndOfFile || atKeyword())
      error(
          token, "expected `]`, found " + describe(token),
          "expected `]` before this",
          "the `[` on line " + std::to_string(open.line) + " is never closed");
    bool has_sign = false;
    double term_sign = sign;
    while (isSign(lexer_.peek())) {
      if (lexer_.next().kind == TokenKind::kMinus) term_sign = -term_sign;
      has_sign = true;
    }
    if (!first && !has_sign)
      error(token, "expected `+`, `-`, or `]`, found " + describe(token),
            "expected `+` or `-` before this",
            "terms in an expression must be separated by `+` or `-`");
    parseQuadraticTerm(term_sign, 2);
    first = false;
  }
  lexer_.next();  // ]
  if (lexer_.peek().kind == TokenKind::kDivide) {
    lexer_.next();
    const Token two = lexer_.peek();
    if (two.kind != TokenKind::kNumber || two.value != 2)
      error(two, "expected `2`, found " + describe(two),
            "expected `2` after `/`",
            "the quadratic part of the objective can only be divided by 2");
    lexer_.next();
    for (size_t i = first_entry; i < hessian_entries_.size(); i++)
      hessian_entries_[i].value /= 2;
  }
}

// <quad-term> := [<constant> ["*"]] <identifier>
//                ("^" "2" | "*" <identifier>)
void Parser::parseQuadraticTerm(double sign, double scale) {
  double coef = sign * scale;
  if (lexer_.peek().kind == TokenKind::kNumber) {
    coef *= lexer_.next().value;
    if (lexer_.peek().kind == TokenKind::kTimes) lexer_.next();
  }
  const HighsInt col1 = parseVariable("expected a variable in this term");
  const Token op = lexer_.peek();
  HighsInt col2 = col1;
  if (op.kind == TokenKind::kPower) {
    lexer_.next();
    const Token two = lexer_.peek();
    if (two.kind != TokenKind::kNumber || two.value != 2)
      error(two, "expected `2`, found " + describe(two),
            "expected `2` after `^`",
            "only squared terms like `x ^ 2` are supported");
    lexer_.next();
  } else if (op.kind == TokenKind::kTimes) {
    lexer_.next();
    col2 = parseVariable("expected a variable after `*`");
  } else {
    error(op, "expected `^` or `*`, found " + describe(op),
          "expected `^ 2` or `* <variable>`",
          "linear terms must be outside `[` and `]`");
  }
  addHessianEntry(col1, col2, coef);
}

// Add the term coef * x[col1] * x[col2] to x'Qx
void Parser::addHessianEntry(HighsInt col1, HighsInt col2, double coef) {
  if (col1 == col2) {
    hessian_entries_.push_back({col1, col1, coef});
  } else {
    hessian_entries_.push_back({col1, col2, coef / 2});
    hessian_entries_.push_back({col2, col1, coef / 2});
  }
}

// Report constraints that are valid LP syntax, but not supported by HiGHS
void Parser::checkUnsupportedConstraint() {
  const Token& token = lexer_.peek(0);
  // <name> S1:: x1:1 x2:2
  if ((token.kind == TokenKind::kIdentifier ||
       token.kind == TokenKind::kNumber) &&
      lexer_.peek(1).kind == TokenKind::kColon &&
      lexer_.peek(2).kind == TokenKind::kColon)
    error(token, "SOS constraints are not supported by HiGHS",
          "SOS constraint");
  // z = 1 -> x + y <= 1
  if (token.kind == TokenKind::kIdentifier &&
      lexer_.peek(1).kind == TokenKind::kEqual &&
      lexer_.peek(2).kind == TokenKind::kNumber &&
      lexer_.peek(3).kind == TokenKind::kImplies)
    error(lexer_.peek(3), "indicator constraints are not supported by HiGHS",
          "indicator constraint");
}

// <constraint> := [<name>] <expression> <comparison> <number> <newline>
void Parser::parseConstraint() {
  checkUnsupportedConstraint();
  std::string name;
  if (parseName(name)) checkUnsupportedConstraint();
  row_has_constant_ = false;
  parseExpression(false);
  const Token comparison = lexer_.peek();
  if (!isComparison(comparison)) {
    std::string help;
    if (comparison.kind == TokenKind::kIdentifier)
      help = describe(comparison) +
             " is a section keyword because it is at the start of a line";
    error(comparison, "expected a comparison, found " + describe(comparison),
          "expected `<=`, `>=`, or `=` before this", help);
  }
  lexer_.next();
  const double rhs = parseNumber("expected the right-hand side");
  const Token& after = lexer_.peek();
  if (after.kind == TokenKind::kIdentifier && !after.line_start)
    error(after, "expected a new line, found " + describe(after),
          "unexpected variable on the right-hand side",
          "move the variables to the left-hand side");
  checkNewLine("expected the constraint to end before this",
               "a constraint ends with a comparison and a number; the next "
               "constraint must start on a new line");
  double lower = -kHighsInf;
  double upper = kHighsInf;
  if (comparison.kind != TokenKind::kLess) lower = rhs;
  if (comparison.kind != TokenKind::kGreater) upper = rhs;
  finishRow(name, lower, upper);
}

void Parser::finishRow(std::string& name, double lower, double upper) {
  HighsSparseMatrix& matrix = lp_.a_matrix_;
  for (const RowEntry& entry : row_entries_) {
    col_slot_[entry.col] = -1;
    if (entry.value == 0) continue;
    matrix.index_.push_back(entry.col);
    matrix.value_.push_back(entry.value);
  }
  row_entries_.clear();
  matrix.start_.push_back(static_cast<HighsInt>(matrix.index_.size()));
  // Give an unnamed row a name, but remember if it could clash with a name
  // in the file
  if (name.compare(0, kRowPrefix.size(), kRowPrefix) == 0) {
    prefix_ok_ = false;
  } else if (name.empty()) {
    name = kRowPrefix + std::to_string(lp_.row_names_.size());
    used_prefix_ = true;
  }
  lp_.row_names_.push_back(std::move(name));
  lp_.row_lower_.push_back(lower);
  lp_.row_upper_.push_back(upper);
}

// <bound> := <identifier> "free" <newline>
//          | <identifier> <comparison> <number> <newline>
//          | <number> <comparison> <identifier> [<comparison> <number>]
//            <newline>
void Parser::parseBound() {
  const Token token = lexer_.peek();
  const bool infinity =
      equalsIgnoreCase(token, "inf") || equalsIgnoreCase(token, "infinity");
  if (token.kind == TokenKind::kIdentifier && !infinity) {
    // x free | x <= 1
    const HighsInt col = getColumn(lexer_.next());
    const Token op = lexer_.peek();
    if (op.kind == TokenKind::kIdentifier && equalsIgnoreCase(op, "free")) {
      lexer_.next();
      lp_.col_lower_[col] = -kHighsInf;
      lp_.col_upper_[col] = kHighsInf;
    } else if (isComparison(op)) {
      lexer_.next();
      const double value = parseNumber("expected the value of the bound");
      if (op.kind != TokenKind::kLess) lp_.col_lower_[col] = value;
      if (op.kind != TokenKind::kGreater) lp_.col_upper_[col] = value;
    } else {
      error(op, "expected `free` or a comparison, found " + describe(op),
            "expected `free`, `<=`, `>=`, or `=`");
    }
  } else {
    // 1 <= x | 1 <= x <= 2
    if (!isSign(token) && token.kind != TokenKind::kNumber && !infinity)
      error(token, "expected a bound, found " + describe(token),
            "expected a variable or a number");
    const double value = parseNumber("expected a number");
    const Token op = lexer_.peek();
    if (!isComparison(op))
      error(op, "expected a comparison, found " + describe(op),
            "expected `<=`, `>=`, or `=`");
    lexer_.next();
    const HighsInt col = parseVariable("expected a variable");
    const Token op2 = lexer_.peek();
    if (isComparison(op2)) {
      // 1 <= x <= 2 | 2 >= x >= 1
      if (op.kind == TokenKind::kEqual)
        error(op, "a bound on both sides of a variable cannot use `=`",
              "expected `<=` or `>=`");
      if (op2.kind != op.kind)
        error(op2,
              "the comparisons in a bound on both sides of a variable must "
              "match",
              op.kind == TokenKind::kLess ? "expected `<=`" : "expected `>=`");
      lexer_.next();
      const double value2 = parseNumber("expected the value of the bound");
      lp_.col_lower_[col] = op.kind == TokenKind::kLess ? value : value2;
      lp_.col_upper_[col] = op.kind == TokenKind::kLess ? value2 : value;
    } else {
      if (op.kind != TokenKind::kGreater) lp_.col_lower_[col] = value;
      if (op.kind != TokenKind::kLess) lp_.col_upper_[col] = value;
    }
  }
  checkNewLine("expected the bound to end before this",
               "each bound must be on its own line");
}

// { <identifier> } in the general, binary and semi-continuous sections
void Parser::parseVariableList(Section section) {
  while (!atSectionEnd()) {
    const HighsInt col = parseVariable("expected a variable");
    HighsVarType& type = lp_.integrality_[col];
    const bool semi = type == HighsVarType::kSemiContinuous ||
                      type == HighsVarType::kSemiInteger;
    const bool integer =
        type == HighsVarType::kInteger || type == HighsVarType::kSemiInteger;
    if (section == Section::kSemi) {
      type =
          integer ? HighsVarType::kSemiInteger : HighsVarType::kSemiContinuous;
    } else {
      type = semi ? HighsVarType::kSemiInteger : HighsVarType::kInteger;
      // A binary variable has an upper bound of 1, unless the bounds section
      // gives it one
      if (section == Section::kBinary) binary_cols_.push_back(col);
    }
    has_integrality_ = true;
  }
}

void Parser::finish(const HighsLogOptions& log_options) {
  for (const HighsInt col : binary_cols_)
    if (lp_.col_upper_[col] == kHighsInf) lp_.col_upper_[col] = 1;
  if (!has_integrality_) lp_.integrality_.clear();

  if (used_prefix_ && !prefix_ok_) {
    lp_.row_names_.clear();
    highsLogUser(log_options, HighsLogType::kWarning,
                 "Cannot create row name beginning \"HiGHS_R\" due to others "
                 "with same prefix: row names cleared\n");
  }

  lp_.num_col_ = static_cast<HighsInt>(lp_.col_names_.size());
  lp_.num_row_ = static_cast<HighsInt>(lp_.row_lower_.size());
  HighsSparseMatrix& matrix = lp_.a_matrix_;
  matrix.num_col_ = lp_.num_col_;
  matrix.num_row_ = lp_.num_row_;
  matrix.ensureColwise();

  // Sum repeated entries of the Hessian, and store it column-wise
  std::sort(hessian_entries_.begin(), hessian_entries_.end());
  hessian_.clear();
  hessian_.start_.clear();
  size_t i = 0;
  while (i < hessian_entries_.size()) {
    const HessianEntry& entry = hessian_entries_[i];
    double value = 0;
    for (;
         i < hessian_entries_.size() && hessian_entries_[i].col == entry.col &&
         hessian_entries_[i].row == entry.row;
         i++)
      value += hessian_entries_[i].value;
    if (value == 0) continue;
    while (static_cast<HighsInt>(hessian_.start_.size()) <= entry.col)
      hessian_.start_.push_back(static_cast<HighsInt>(hessian_.index_.size()));
    hessian_.index_.push_back(entry.row);
    hessian_.value_.push_back(value);
  }
  if (hessian_.index_.empty()) {
    hessian_.clear();
  } else {
    hessian_.dim_ = lp_.num_col_;
    while (static_cast<HighsInt>(hessian_.start_.size()) <= lp_.num_col_)
      hessian_.start_.push_back(static_cast<HighsInt>(hessian_.index_.size()));
    hessian_.format_ = HessianFormat::kSquare;
  }
}

}  // namespace lp_reader

FilereaderRetcode readLpFile(const HighsLogOptions& log_options,
                             const std::string& filename, HighsModel& model) {
  std::unique_ptr<std::istream> file = lp_reader::openFile(filename);
  if (file->fail()) return FilereaderRetcode::kFileNotFound;
  lp_reader::Parser parser(*file, model);
  try {
    parser.parse();
  } catch (const lp_reader::ParseError& e) {
    file.reset();
    highsLogUser(
        log_options, HighsLogType::kError, "%s",
        lp_reader::renderDiagnostics(filename, {e.diagnostic})[0].c_str());
    return FilereaderRetcode::kParserError;
  } catch (const std::ios_base::failure& e) {
    // For example, zstr throws if a compressed file is corrupt
    highsLogUser(log_options, HighsLogType::kError, "Error reading %s: %s\n",
                 filename.c_str(), e.what());
    return FilereaderRetcode::kParserError;
  }
  file.reset();
  parser.finish(log_options);
  for (const std::string& warning :
       lp_reader::renderDiagnostics(filename, parser.warnings()))
    highsLogUser(log_options, HighsLogType::kWarning, "%s", warning.c_str());
  if (parser.numWarnings() > lp_reader::kMaxWarnings) {
    const size_t num_hidden = parser.numWarnings() - lp_reader::kMaxWarnings;
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%d more warning%s not shown\n", int(num_hidden),
                 num_hidden == 1 ? " was" : "s were");
  }
  return parser.numWarnings() == 0 ? FilereaderRetcode::kOk
                                   : FilereaderRetcode::kWarning;
}
