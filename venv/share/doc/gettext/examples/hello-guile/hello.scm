#!@GUILE@ -s
!#
;;; Example for use of GNU gettext.
;;; This file is in the public domain.

;;; Source code of the GNU guile program.

(use-modules (ice-9 format))

(catch #t (lambda () (setlocale LC_ALL "")) (lambda args #f))
(textdomain "hello-guile")
(bindtextdomain "hello-guile" "@localedir@")
(define _ gettext)

(display (_ "Hello, world!"))
(newline)
(format #t (_ "This program is running as process number ~D.") (getpid))
(newline)
