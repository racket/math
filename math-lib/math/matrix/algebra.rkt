#lang racket/base

(require typed/untyped-utils)

(require "../private/matrix/matrix-types.rkt"
         (except-in "../private/matrix/algebra/matrix-solve.rkt"
                    make-matrix-determinant))

(require/untyped-contract
 (begin (require "../private/matrix/matrix-types.rkt"))
 "../private/matrix/algebra/matrix-solve.rkt"
 [make-matrix-determinant
  ((Any * -> Any) (Any -> Any)
   (Any * -> Any) (Any -> Any)
   (Any Any -> Boolean)
   (Any Any -> Boolean)
   ->
   ((Matrix Any) -> Any))])

(provide (all-from-out "../private/matrix/matrix-types.rkt"
                       "../private/matrix/algebra/matrix-solve.rkt")
         ;; matrix-solve.rkt
         make-matrix-determinant)
