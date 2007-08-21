
#include <iostream>
#include <fstream>
#include "bzfstream.h"
#include "optuple.h"

#ifdef NO_BOOL
 #define bool int
 #define true 1
 #define false 0
#endif

int main()
{
    std::cout << "Generating <vecbops.cc>" << std::endl;

    OperandTuple operands(2);

    bzofstream ofs("../vecbops.cc", "Vector expression templates (2 operands)",
        __FILE__, "BZ_VECBOPS_CC");

    ofs << "#ifndef BZ_VECEXPR_H" << std::endl
        << " #error <blitz/vecbops.cc> must be included via <blitz/vecexpr.h>" 
        << std::endl << "#endif" << std::endl << std::endl;

    ofs.beginNamespace();

    struct {
        const char* opSymbol;
        bool        nonIntOperands;
        bool        complexOperands;
        const char* opApplicName;
        const char* comment;
    } ops[] = {
     { "+",  true,  true,  "_bz_Add",            "Addition Operators" },
     { "-",  true,  true,  "_bz_Subtract",       "Subtraction Operators" },
     { "*",  true,  true,  "_bz_Multiply",       "Multiplication Operators" },
     { "/",  true,  true,  "_bz_Divide",         "Division Operators" },
     { "%",  false, false, "_bz_Mod",            "Modulus Operators" },
     { "^",  false, false, "_bz_BitwiseXOR",     "Bitwise XOR Operators" },
     { "&",  false, false, "_bz_BitwiseAnd",     "Bitwise And Operators" },
     { "|",  false, false, "_bz_BitwiseOr",      "Bitwise Or Operators" },
     { ">>", false, false, "_bz_ShiftRight",     "Shift right Operators" },
     { "<<", false, false, "_bz_ShiftLeft",      "Shift left Operators" },
     { ">",  true,  false, "_bz_Greater",        "Greater-than Operators" },
     { "<",  true,  false, "_bz_Less",           "Less-than Operators" },
     { ">=", true,  false, "_bz_GreaterOrEqual", "Greater or equal (>=) operators" },
     { "<=", true,  false, "_bz_LessOrEqual",    "Less or equal (<=) operators" },
     { "==", true,  true,  "_bz_Equal",          "Equality operators" },
     { "!=", true,  true,  "_bz_NotEqual",       "Not-equal operators" },
     { "&&", false, false, "_bz_LogicalAnd",     "Logical AND operators" },
     { "||", false, false, "_bz_LogicalOr",      "Logical OR operators" }
    };

    for (int i=0; i < 18; ++i)
    {
    ofs << "/****************************************************************************" << std::endl
        << " * " << ops[i].comment << std::endl
        << " ****************************************************************************/" << std::endl;

    operands.reset();

    do {
        // Can't declare operator+(int,Range) or operator+(Range,int)
        // because these would conflict with the versions defined
        // in range.h.  Also, the versions in range.h will be
        // much faster.
        if (operands[0].isScalar() && operands[0].isInteger()
            && operands[1].isRange())
                continue;
        if (operands[1].isScalar() && operands[1].isInteger()
            && operands[0].isRange())
                continue;

        if (ops[i].nonIntOperands == false)
        {
            if ((operands[0].isScalar() && !operands[0].isInteger())
             || (operands[1].isScalar() && !operands[1].isInteger()))
                continue;
        }

        // Vector<P_numtype1> + _bz_VecExpr<P_expr2>

        if (operands.anyComplex())
            ofs << "#ifdef BZ_HAVE_COMPLEX" << std::endl;

        ofs << std::endl << "// ";
        operands[0].printName(ofs);
        ofs << " " << ops[i].opSymbol << " ";
        operands[1].printName(ofs);
        ofs << std::endl;

        operands.printTemplates(ofs);
        ofs << std::endl << "inline" << std::endl;

        // _bz_VecExpr<_bz_VecExprOp<VectorIterConst<T_numtype1>,
        //     _bz_VecExpr<T_expr2>, _bz_Add<T_numtype1,typename T_expr2::T_numtype> > >
        ofs << "_bz_VecExpr<_bz_VecExprOp<";
        operands.printIterators(ofs, 1);
        ofs << "," << std::endl << "      " << ops[i].opApplicName << "<";
        operands[0].printNumtype(ofs);
        ofs << ", ";    
        operands[1].printNumtype(ofs);
        ofs << " > > >" << std::endl;
     
        // operator+(const Vector<T_numtype1>& d1, _bz_VecExpr<T_expr2> d2)
				if (ops[i].opSymbol[0] == 'm')
					ofs << ops[i].opSymbol << "(";
				else
        	ofs << "operator" << ops[i].opSymbol << "(";
        operands.printArgumentList(ofs, 1);
        ofs << ")" << std::endl << "{" << std::endl;

        // typedef _bz_VecExprOp<VectorIterConst<T_numtype1>,
        // _bz_VecExpr<T_expr2>, _bz_Add<T_numtype1,typename T_expr2::T_numtype> > T_expr;
        ofs << "    typedef _bz_VecExprOp<";
        operands.printIterators(ofs, 1);
        ofs << ", " << std::endl << "      " << ops[i].opApplicName << "<";
        operands[0].printNumtype(ofs);
        ofs << ", ";
        operands[1].printNumtype(ofs);
        ofs << "> > T_expr;" << std::endl << std::endl;

        // return _bz_VecExpr<T_expr>(T_expr(a.begin(), b));
        ofs << "    return _bz_VecExpr<T_expr>(T_expr(";
        operands.printInitializationList(ofs,1);
        ofs << "));" << std::endl << "}" << std::endl;

        if (operands.anyComplex())
            ofs << "#endif // BZ_HAVE_COMPLEX" << std::endl << std::endl;

    } while (++operands);

   }

   std::cout << operands.numSpecializations() << " operators written." << std::endl;
}

