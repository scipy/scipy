;;; -*- Mode: lisp -*-

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

(defpackage #:gpg-error-system
  (:use #:common-lisp #:asdf))

(in-package #:gpg-error-system)

(defsystem gpg-error
    :description "Common error values for all GnuPG components."
    :author "g10 Code GmbH"
    :version "1.47-unknown"
    :licence "LGPL"
    :depends-on ("cffi")
    :components ((:file "gpg-error-package")
		 (:file "gpg-error-codes"
			:depends-on ("gpg-error-package"))
		 (:file "gpg-error" :depends-on ("gpg-error-codes"))))
