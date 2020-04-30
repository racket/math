#lang typed/racket/base

(require "matrix-types.rkt"
         "matrix-comprehension.rkt"
         "matrix-basic.rkt"
         "utils.rkt")

(provide matrix-kronecker)

;; Based on Tim Brown's code. Used with permission.
(: matrix-kronecker/2 (All (A) ((Matrix A) (Matrix A) (A A -> A) -> (Matrix A))))
(define (matrix-kronecker/2 a b ×)
  (define-values (row-a col-a) (matrix-shape a))
  (define-values (row-b col-b) (matrix-shape b))
  (define row-out (* row-a row-b))
  (define col-out (* col-a col-b))
  (for*/matrix: row-out col-out ([r (in-range row-out)] [c (in-range col-out)]) : A
    (define-values (index-row-a index-row-b) (quotient/remainder r row-b))
    (define-values (index-col-a index-col-b) (quotient/remainder c col-b))
    (× (matrix-ref a index-row-a index-col-a) (matrix-ref b index-row-b index-col-b))))

(: matrix-kronecker/ns
   (case-> ((Matrix Flonum) (Listof (Matrix Flonum)) -> (Matrix Flonum))
           ((Matrix Real) (Listof (Matrix Real)) -> (Matrix Real))
           ((Matrix Float-Complex) (Listof (Matrix Float-Complex)) -> (Matrix Float-Complex))
           ((Matrix Number) (Listof (Matrix Number)) -> (Matrix Number))))
(define (matrix-kronecker/ns a as)
  (for/fold ([result a]) ([x (in-list as)])
    (matrix-kronecker/2 result x *)))

(: matrix-kronecker
   (case-> ((Matrix Flonum) (Matrix Flonum) * -> (Matrix Flonum))
           ((Matrix Real) (Matrix Real) * -> (Matrix Real))
           ((Matrix Float-Complex) (Matrix Float-Complex) * -> (Matrix Float-Complex))
           ((Matrix Number) (Matrix Number) * -> (Matrix Number))))
(define (matrix-kronecker a . as) (call/ns (λ () (matrix-kronecker/ns a as))))
