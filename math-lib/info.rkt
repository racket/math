#lang info

(define collection 'multi)

(define deps '(["base" #:version "8.16.0.4"]
               ["typed-racket-lib" #:version "1.14"]
               ("math-i386-macosx" #:platform "i386-macosx")
               ("math-x86_64-macosx" #:platform "x86_64-macosx")
               ("math-ppc-macosx" #:platform "ppc-macosx")
               ("math-aarch64-macosx" #:platform "aarch64-macosx")
               ("math-win32-i386" #:platform "win32\\i386")
               ("math-win32-x86_64" #:platform "win32\\x86_64")
               ("math-win32-arm64" #:platform "win32\\arm64")
               ("math-x86_64-linux-natipkg" #:platform "x86_64-linux-natipkg")))

(define build-deps '())

(define pkg-desc "Math library")

(define pkg-authors '(ntoronto))

(define version "1.2")

(define license
  '(Apache-2.0 OR MIT))
