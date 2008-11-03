/***************************************************************************
 * blitz/generate/bzfstream.h    Definition of the bzofstream class
 *
 * $Id: bzfstream.h 1413 2005-11-01 22:04:15Z cookedm $
 *
 * Copyright (C) 1997-2001 Todd Veldhuizen <tveldhui@oonumerics.org>
 *
 * This code was relicensed under the modified BSD license for use in SciPy
 * by Todd Veldhuizen (see LICENSE.txt in the weave directory).
 *
 *
 * Suggestions:          blitz-dev@oonumerics.org
 * Bugs:                 blitz-bugs@oonumerics.org
 *
 * For more information, please see the Blitz++ Home Page:
 *    http://oonumerics.org/blitz/
 *
 ***************************************************************************/


#include <fstream>
#include <iomanip>
#include <iostream>

class bzofstream : public std::ofstream {

public:
    bzofstream(const char* filename, const char* description,
        const char* sourceFile, const char* mnemonic)
        : std::ofstream(filename)
    {
        (*this) << 
"/***************************************************************************\n"
" * blitz/" << filename << "\t" << description << std::endl <<
" *\n"
" * This code was relicensed under the modified BSD license for use in SciPy\n"
" * by Todd Veldhuizen (see LICENSE.txt in the weave directory).\n"
" *\n"
" *\n"
" * Suggestions:          blitz-suggest@cybervision.com\n"
" * Bugs:                 blitz-bugs@cybervision.com\n"
" *\n"
" * For more information, please see the Blitz++ Home Page:\n"
" *    http://seurat.uwaterloo.ca/blitz/\n"
" *\n"
" ***************************************************************************\n"
" *\n"
" */ " 
       << std::endl << std::endl
       << "// Generated source file.  Do not edit. " << std::endl
       << "// " << sourceFile << " " << __DATE__ << " " << __TIME__ 
       << std::endl << std::endl
       << "#ifndef " << mnemonic << std::endl
       << "#define " << mnemonic << std::endl << std::endl;
    }

    void include(const char* filename)
    {
        (*this) << "#include <blitz/" << filename << ">" << std::endl;
    }

    void beginNamespace()
    {
        (*this) << "BZ_NAMESPACE(blitz)" << std::endl << std::endl;
    }

    ~bzofstream()
    {
        (*this) << "BZ_NAMESPACE_END" << std::endl << std::endl
                << "#endif" << std::endl;
    }

};

