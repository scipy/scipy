#include "head.c"

double f1_pypy_g_bar(long l_n_2) {
	long l_n_3; bool_t l_v190; bool_t l_v192; bool_t l_v196;
	double l_v191; double l_v194; double l_v195; double l_v199;
	double l_v200; long l_v187; long l_v188; long l_v189; long l_v193;
	long l_v197; long l_v198; long l_v201;

    block0:
	OP_INT_LT(l_n_2, 10L, l_v190);
	if (l_v190) {
		l_v201 = 0L;
		l_n_3 = l_n_2;
		l_v188 = 1L;
		l_v189 = 10L;
		goto block6;
	}
	goto block1;

    block1:
	switch (l_n_2) {
    case 10L:
		l_v200 = 42.0;
		goto block2;
    case 11L:
		l_v200 = 3.1400000000000001;
		goto block2;
    case 12L:
		goto block3;
    default:
		goto block4;
	}

    block2:
	RPY_DEBUG_RETURN();
	return l_v200;

    block3:
	l_v191 = cos(3.1415926535897931);
	l_v200 = l_v191;
	goto block2;

    block4:
	OP_INT_GT(l_n_2, 12L, l_v192);
	if (l_v192) {
		goto block5;
	}
	l_v200 = 5.0;
	goto block2;

    block5:
	OP_INT_MUL(l_n_2, l_n_2, l_v193);
	OP_CAST_INT_TO_FLOAT(l_v193, l_v194);
	OP_FLOAT_TRUEDIV(l_v194, 1.2345600000000001, l_v195);
	l_v200 = l_v195;
	goto block2;

    block6:
	OP_INT_ADD(l_n_3, l_v201, l_v187);
	OP_INT_GE(l_v188, l_v189, l_v196);
	while (!l_v196) {
		goto block7;
    	  block6_back:
		OP_INT_ADD(l_n_3, l_v201, l_v187);
		OP_INT_GE(l_v188, l_v189, l_v196);
	}
	goto block8;

    block7:
	OP_INT_ADD(l_v188, 1L, l_v197);
	OP_INT_MUL(l_v188, l_v188, l_v198);
	l_v201 = l_v198;
	l_n_3 = l_v187;
	l_v188 = l_v197;
	goto block6_back;

    block8:
	OP_CAST_INT_TO_FLOAT(l_v187, l_v199);
	l_v200 = l_v199;
	goto block2;
}

double f2_pypy_g_bar(double l_n_6) {
	double l_n_7; bool_t l_v546; bool_t l_v547; bool_t l_v548;
	bool_t l_v549; bool_t l_v550; bool_t l_v554; double l_v545;
	double l_v551; double l_v552; double l_v553; double l_v557;
	double l_v558; double l_v559; long l_v543; long l_v544; long l_v555;
	long l_v556;

    block0:
	OP_FLOAT_LT(l_n_6, 10.0, l_v546);
	if (l_v546) {
		l_v559 = 0.0;
		l_v543 = 10L;
		l_n_7 = l_n_6;
		l_v544 = 1L;
		goto block8;
	}
	goto block1;

    block1:
	OP_FLOAT_EQ(l_n_6, 10.0, l_v547);
	if (l_v547) {
		l_v558 = 42.0;
		goto block5;
	}
	goto block2;

    block2:
	OP_FLOAT_EQ(l_n_6, 11.0, l_v548);
	if (l_v548) {
		l_v558 = 3.1400000000000001;
		goto block5;
	}
	goto block3;

    block3:
	OP_FLOAT_EQ(l_n_6, 12.0, l_v549);
	if (l_v549) {
		goto block7;
	}
	goto block4;

    block4:
	OP_FLOAT_GT(l_n_6, 12.0, l_v550);
	if (l_v550) {
		goto block6;
	}
	l_v558 = 5.0;
	goto block5;

    block5:
	RPY_DEBUG_RETURN();
	return l_v558;

    block6:
	OP_FLOAT_MUL(l_n_6, l_n_6, l_v551);
	OP_FLOAT_TRUEDIV(l_v551, 1.2345600000000001, l_v552);
	l_v558 = l_v552;
	goto block5;

    block7:
	l_v553 = cos(3.1415926535897931);
	l_v558 = l_v553;
	goto block5;

    block8:
	OP_FLOAT_ADD(l_n_7, l_v559, l_v545);
	OP_INT_GE(l_v544, l_v543, l_v554);
	while (!l_v554) {
		goto block9;
    	  block8_back:
		OP_FLOAT_ADD(l_n_7, l_v559, l_v545);
		OP_INT_GE(l_v544, l_v543, l_v554);
	}
	l_v558 = l_v545;
	goto block5;

    block9:
	OP_INT_ADD(l_v544, 1L, l_v555);
	OP_INT_MUL(l_v544, l_v544, l_v556);
	OP_CAST_INT_TO_FLOAT(l_v556, l_v557);
	l_v559 = l_v557;
	l_n_7 = l_v545;
	l_v544 = l_v555;
	goto block8_back;
}

