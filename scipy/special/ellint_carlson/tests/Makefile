.PHONY: test test_complex test_real clean

TEST_COMPLEX_EXECS = test_complex_for_rf.t test_complex_for_rd.t test_complex_for_rc.t test_complex_for_rj.t test_complex_for_rg.t
TEST_REAL_EXECS = test_real_for_rf.t test_real_for_rd.t test_real_for_rc.t test_real_for_rj.t test_real_for_rg.t

test: test_complex test_real

test_complex: $(TEST_COMPLEX_EXECS)
	@for testexec in $^; do echo "$${testexec}"; "./$${testexec}" || exit "$$?"; done

test_real: $(TEST_REAL_EXECS)
	@for testexec in $^; do echo "$${testexec}"; "./$${testexec}" || exit "$$?"; done

test_common_complex.o: test_common.c test_common.h
	$(CC) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) -c -fPIC -DELLINT_POLY_COMPLEX $< -o $@ -I. -I..

test_common_real.o: test_common.c test_common.h
	$(CC) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) -c -fPIC -DELLINT_POLY_REAL $< -o $@ -I. -I..

libellint_carlson.so:
	ln -s ../$@

test_complex_for_%.t: test_complex_for_%.c test_common_complex.o libellint_carlson.so
	$(CC) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $+ -o $@ -lm -I. -I.. -L..

test_real_for_%.t: test_real_for_%.c test_common_real.o libellint_carlson.so
	$(CC) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $+ -o $@ -lm -I. -I.. -L..

clean:
	rm -rf test_*_for_*.t *.dSYM *.o *.so
