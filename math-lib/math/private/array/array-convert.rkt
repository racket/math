#lang racket/base

(require typed/untyped-utils
         (except-in "typed-array-convert.rkt"
                    list*->array
                    vector*->array)
         (prefix-in typed: (only-in "typed-array-convert.rkt"
                                    list*->array
                                    vector*->array))
         (prefix-in untyped: (only-in "untyped-array-convert.rkt"
                                      list*->array
                                      vector*->array)))

(provide list*->array
         vector*->array
         array->list*
         array->vector*
         array->list
         array->vector)

(module shallow-defs typed/racket/shallow
  (require typed/racket/unsafe "array-struct.rkt" "mutable-array.rkt" "utils.rkt")
  (provide shallow:list*->array shallow:vector*->array)
  (unsafe-require/typed "untyped-array-convert.rkt"
    [(list*->array shallow:list*->array)
     (All (A) ((Listof* A) ((Listof* A) -> Any : A) -> (Array A)))]
    [(vector*->array shallow:vector*->array)
     (All (A) ((Vectorof* A) ((Vectorof* A) -> Any : A) -> (Mutable-Array A)))]))
(require 'shallow-defs)

(define-typed/untyped-identifier list*->array
  typed:list*->array
  untyped:list*->array
  shallow:list*->array
  shallow:list*->array)

(define-typed/untyped-identifier vector*->array
  typed:vector*->array
  untyped:vector*->array
  shallow:vector*->array
  shallow:vector*->array)
