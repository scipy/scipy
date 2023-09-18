;;;; libgpg-error-package.lisp

;;; Copyright (C) 2006 g10 Code GmbH
;;;
;;; This file is part of libgpg-error.
;;;
;;; libgpg-error is free software; you can redistribute it and/or
;;; modify it under the terms of the GNU Lesser General Public License
;;; as published by the Free Software Foundation; either version 2.1 of
;;; the License, or (at your option) any later version.
;;;
;;; libgpg-error is distributed in the hope that it will be useful, but
;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;; Lesser General Public License for more details.
;;;
;;; You should have received a copy of the GNU Lesser General Public
;;; License along with libgpg-error; if not, write to the Free
;;; Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
;;; 02111-1307, USA.

;;; Conventions
;;;
;;; Error sources and codes are represented as keywords like
;;; :gpg-err-source-gpg and :gpg-err-unknown-packet.
;;;
;;; Errors are represented as lists '(SOURCE CODE).  Other
;;; representations are also accepted in some places.
;;;
;;; The following functions are defined which are not defined in the C API:
;;; gpg-err-source-as-key, gpg-err-source-as-value
;;; gpg-err-code-as-key, gpg-err-code-as-value
;;; gpg-err-canonicalize, gpg-err-as-value
;;; Conversion between keywords and values for error sources and codes.
;;; 
;;; The following functions from the C API are omitted:
;;; gpg-strerror-r
;;;
;;; The following features work slightly differently:
;;; *gpg-err-source-default* is a dynamic variable that can be set to
;;; change the default for gpg-error.

(defpackage #:gpg-error
  (:use #:common-lisp #:cffi)

  (:export :gpg-err-code-as-key
	   :gpg-err-code-as-value
	   :gpg-err-source-as-key
	   :gpg-err-source-as-value
	   :gpg-err-canonicalize
	   :gpg-err-as-value
	   :gpg-err-make
	   :*gpg-err-source-default*
	   :gpg-error
	   :gpg-err-code
	   :gpg-err-source
	   :gpg-strerror
	   :gpg-strsource
	   :gpg-err-code-from-errno
	   :gpg-err-code-to-errno
           :gpg-err-code-from-syserror
	   :gpg-err-make-from-errno
	   :gpg-error-from-errno
           :gpg-error-from-syserror))
