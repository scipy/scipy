;;;; libgpg-error.lisp

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

;;; Set up the library.

(in-package :gpg-error)

(define-foreign-library libgpg-error
  (:unix "libgpg-error.so")
  (t (:default "libgpg-error")))
   
(use-foreign-library libgpg-error)

;;; System dependencies.

(defctype size-t :unsigned-int "The system size_t type.")

;;; Error sources.

(defcenum gpg-err-source-t
  "The GPG error source type."
  (:gpg-err-source-unknown 0)
  (:gpg-err-source-gcrypt 1)
  (:gpg-err-source-gpg 2)
  (:gpg-err-source-gpgsm 3)
  (:gpg-err-source-gpgagent 4)
  (:gpg-err-source-pinentry 5)
  (:gpg-err-source-scd 6)
  (:gpg-err-source-gpgme 7)
  (:gpg-err-source-keybox 8)
  (:gpg-err-source-ksba 9)
  (:gpg-err-source-dirmngr 10)
  (:gpg-err-source-gsti 11)
  (:gpg-err-source-any 31)
  (:gpg-err-source-user-1 32)
  (:gpg-err-source-user-2 33)
  (:gpg-err-source-user-3 34)
  (:gpg-err-source-user-4 35))

(defconstant +gpg-err-source-dim+ 256)

;;; The error code type gpg-err-code-t.

;;; libgpg-error-codes.lisp is loaded by ASDF.

(defctype gpg-error-t :unsigned-int "The GPG error code type.")

;;; Bit mask manipulation constants.

(defconstant +gpg-err-code-mask+ (- +gpg-err-code-dim+ 1))

(defconstant +gpg-err-source-mask+ (- +gpg-err-source-dim+ 1))
(defconstant +gpg-err-source-shift+ 24)

;;; Constructor and accessor functions.

;;; If we had in-library versions of our static inlines, we wouldn't
;;; need to replicate them here.  Oh well.

(defun c-gpg-err-make (source code)
  "Construct an error value from an error code and source.
   Within a subsystem, use gpg-error instead."
  (logior
   (ash (logand source +gpg-err-source-mask+)
	+gpg-err-source-shift+)
   (logand code +gpg-err-code-mask+)))

(defun c-gpg-err-code (err)
  "retrieve the error code from an error value." 
  (logand err +gpg-err-code-mask+))

(defun c-gpg-err-source (err)
  "retrieve the error source from an error value." 
  (logand (ash err (- +gpg-err-source-shift+))
	  +gpg-err-source-mask+))

;;; String functions.

(defcfun ("gpg_strerror" c-gpg-strerror) :string
  (err gpg-error-t))

(defcfun ("gpg_strsource" c-gpg-strsource) :string
  (err gpg-error-t))

;;; Mapping of system errors (errno).

(defcfun ("gpg_err_code_from_errno" c-gpg-err-code-from-errno) gpg-err-code-t
  (err :int))

(defcfun ("gpg_err_code_to_errno" c-gpg-err-code-to-errno) :int
  (code gpg-err-code-t))

(defcfun ("gpg_err_code_from_syserror"
           c-gpg-err-code-from-syserror) gpg-err-code-t)

;;; Self-documenting convenience functions.

;;; See below.

;;;
;;;
;;; Lispy interface.
;;;
;;;

;;; Low-level support functions.

(defun gpg-err-code-as-value (code-key)
  (foreign-enum-value 'gpg-err-code-t code-key))

(defun gpg-err-code-as-key (code)
  (foreign-enum-keyword 'gpg-err-code-t code))

(defun gpg-err-source-as-value (source-key)
  (foreign-enum-value 'gpg-err-source-t source-key))

(defun gpg-err-source-as-key (source)
  (foreign-enum-keyword 'gpg-err-source-t source))

(defun gpg-err-canonicalize (err)
  "Canonicalize the error value err."
  (gpg-err-make (gpg-err-source err) (gpg-err-code err)))

(defun gpg-err-as-value (err)
  "Get the integer representation of the error value ERR."
  (let ((error (gpg-err-canonicalize err)))
    (c-gpg-err-make (gpg-err-source-as-value (gpg-err-source error))
		    (gpg-err-code-as-value (gpg-err-code error)))))

;;; Constructor and accessor functions.

(defun gpg-err-make (source code)
  "Construct an error value from an error code and source.
   Within a subsystem, use gpg-error instead."
  ;; As an exception to the rule, the function gpg-err-make will use
  ;; the error source value as is when provided as integer, instead of
  ;; parsing it as an error value.
  (list (if (integerp source)
	    (gpg-err-source-as-key source)
	    (gpg-err-source source))
	(gpg-err-code code)))

(defvar *gpg-err-source-default* :gpg-err-source-unknown
  "define this to specify a default source for gpg-error.")

(defun gpg-error (code)
  "Construct an error value from an error code, using the default source."
  (gpg-err-make *gpg-err-source-default* code))

(defun gpg-err-code (err)
    "Retrieve an error code from the error value ERR."
    (cond ((listp err) (second err))
	  ((keywordp err) err) ; FIXME
	  (t (gpg-err-code-as-key (c-gpg-err-code err)))))

(defun gpg-err-source (err)
    "Retrieve an error source from the error value ERR."
    (cond ((listp err) (first err))
	  ((keywordp err) err) ; FIXME
	  (t (gpg-err-source-as-key (c-gpg-err-source err)))))

;;; String functions.

(defun gpg-strerror (err)
  "Return a string containig a description of the error code."
  (c-gpg-strerror (gpg-err-as-value err)))

;;; FIXME: maybe we should use this as the actual implementation for
;;; gpg-strerror.

;; (defcfun ("gpg_strerror_r" c-gpg-strerror-r) :int
;;   (err gpg-error-t)
;;   (buf :string)
;;   (buflen size-t))

;; (defun gpg-strerror-r (err)
;;   "Return a string containig a description of the error code."
;;   (with-foreign-pointer-as-string (errmsg 256 errmsg-size)
;;     (c-gpg-strerror-r (gpg-err-code-as-value (gpg-err-code err))
;; 		      errmsg errmsg-size)))

(defun gpg-strsource (err)
  "Return a string containig a description of the error source."
  (c-gpg-strsource (gpg-err-as-value err)))

;;; Mapping of system errors (errno).

(defun gpg-err-code-from-errno (err)
  "Retrieve the error code for the system error.  If the system error
   is not mapped, :gpg-err-unknown-errno is returned."
  (gpg-err-code-as-key (c-gpg-err-code-from-errno err)))

(defun gpg-err-code-to-errno (code)
  "Retrieve the system error for the error code.  If this is not a
   system error, 0 is returned."
  (c-gpg-err-code-to-errno (gpg-err-code code)))

(defun gpg-err-code-from-syserror ()
  "Retrieve the error code directly from the system ERRNO.  If the system error
   is not mapped, :gpg-err-unknown-errno is returned and 
   :gpg-err-missing-errno if ERRNO has the value 0."
  (gpg-err-code-as-key (c-gpg-err-code-from-syserror)))


;;; Self-documenting convenience functions.

(defun gpg-err-make-from-errno (source err)
  (gpg-err-make source (gpg-err-code-from-errno err)))

(defun gpg-error-from-errno (err)
  (gpg-error (gpg-err-code-from-errno err)))

(defun gpg-error-from-syserror ()
  (gpg-error (gpg-err-code-from-syserror)))

