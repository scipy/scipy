# Increment this whenever this file is changed.
#serial 2

dnl GLIB_GSETTINGS
dnl Defines GSETTINGS_SCHEMAS_INSTALL which controls whether
dnl the schema should be compiled
dnl

AC_DEFUN([GLIB_GSETTINGS],
[
  dnl We can't use PKG_PREREQ because that needs 0.29.
  m4_ifndef([PKG_PROG_PKG_CONFIG],
            [pkg.m4 version 0.28 or later is required])

  m4_pattern_allow([AM_V_GEN])
  AC_ARG_ENABLE(schemas-compile,
                AS_HELP_STRING([--disable-schemas-compile],
                               [Disable regeneration of gschemas.compiled on install]),
                [case ${enableval} in
                  yes) GSETTINGS_DISABLE_SCHEMAS_COMPILE=""  ;;
                  no)  GSETTINGS_DISABLE_SCHEMAS_COMPILE="1" ;;
                  *) AC_MSG_ERROR([bad value ${enableval} for --enable-schemas-compile]) ;;
                 esac])
  AC_SUBST([GSETTINGS_DISABLE_SCHEMAS_COMPILE])
  PKG_PROG_PKG_CONFIG([0.16])
  AC_SUBST(gsettingsschemadir, [${datadir}/glib-2.0/schemas])
  AS_IF([test x$cross_compiling != xyes],
        [PKG_CHECK_VAR([GLIB_COMPILE_SCHEMAS], [gio-2.0], [glib_compile_schemas])],
        [AC_PATH_PROG([GLIB_COMPILE_SCHEMAS], [glib-compile-schemas])])
  AC_SUBST(GLIB_COMPILE_SCHEMAS)
  if test "x$GLIB_COMPILE_SCHEMAS" = "x"; then
    ifelse([$2],,[AC_MSG_ERROR([glib-compile-schemas not found.])],[$2])
  else
    ifelse([$1],,[:],[$1])
  fi

  GSETTINGS_RULES='
.PHONY : uninstall-gsettings-schemas install-gsettings-schemas clean-gsettings-schemas

mostlyclean-am: clean-gsettings-schemas

gsettings__enum_file = $(addsuffix .enums.xml,$(gsettings_ENUM_NAMESPACE))

%.gschema.valid: %.gschema.xml $(gsettings__enum_file)
	$(AM_V_GEN) $(GLIB_COMPILE_SCHEMAS) --strict --dry-run $(addprefix --schema-file=,$(gsettings__enum_file)) --schema-file=$< && mkdir -p [$](@D) && touch [$]@

all-am: $(gsettings_SCHEMAS:.xml=.valid)
uninstall-am: uninstall-gsettings-schemas
install-data-am: install-gsettings-schemas

.SECONDARY: $(gsettings_SCHEMAS)

install-gsettings-schemas: $(gsettings_SCHEMAS) $(gsettings__enum_file)
	@$(NORMAL_INSTALL)
	if test -n "$^"; then \
		test -z "$(gsettingsschemadir)" || $(MKDIR_P) "$(DESTDIR)$(gsettingsschemadir)"; \
		$(INSTALL_DATA) $^ "$(DESTDIR)$(gsettingsschemadir)"; \
		test -n "$(GSETTINGS_DISABLE_SCHEMAS_COMPILE)$(DESTDIR)" || $(GLIB_COMPILE_SCHEMAS) $(gsettingsschemadir); \
	fi

uninstall-gsettings-schemas:
	@$(NORMAL_UNINSTALL)
	@list='\''$(gsettings_SCHEMAS) $(gsettings__enum_file)'\''; test -n "$(gsettingsschemadir)" || list=; \
	files=`for p in $$list; do echo $$p; done | sed -e '\''s|^.*/||'\''`; \
	test -n "$$files" || exit 0; \
	echo " ( cd '\''$(DESTDIR)$(gsettingsschemadir)'\'' && rm -f" $$files ")"; \
	cd "$(DESTDIR)$(gsettingsschemadir)" && rm -f $$files
	test -n "$(GSETTINGS_DISABLE_SCHEMAS_COMPILE)$(DESTDIR)" || $(GLIB_COMPILE_SCHEMAS) $(gsettingsschemadir)

clean-gsettings-schemas:
	rm -f $(gsettings_SCHEMAS:.xml=.valid) $(gsettings__enum_file)

ifdef gsettings_ENUM_NAMESPACE
$(gsettings__enum_file): $(gsettings_ENUM_FILES)
	$(AM_V_GEN) glib-mkenums --comments '\''<!-- @comment@ -->'\'' --fhead "<schemalist>" --vhead "  <@type@ id='\''$(gsettings_ENUM_NAMESPACE).@EnumName@'\''>" --vprod "    <value nick='\''@valuenick@'\'' value='\''@valuenum@'\''/>" --vtail "  </@type@>" --ftail "</schemalist>" [$]^ > [$]@.tmp && mv [$]@.tmp [$]@
endif
'
  _GSETTINGS_SUBST(GSETTINGS_RULES)
])

dnl _GSETTINGS_SUBST(VARIABLE)
dnl Abstract macro to do either _AM_SUBST_NOTMAKE or AC_SUBST
AC_DEFUN([_GSETTINGS_SUBST],
[
AC_SUBST([$1])
m4_ifdef([_AM_SUBST_NOTMAKE], [_AM_SUBST_NOTMAKE([$1])])
]
)
