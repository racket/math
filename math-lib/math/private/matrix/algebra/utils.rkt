#lang typed/racket/base

(require racket/fixnum
         "../base/utils.rkt"
         "../../unsafe.rkt"
         "../../vector/vector-mutate.rkt")

(provide (all-from-out "../base/utils.rkt")
         (all-defined-out))

(: generic-find-partial-pivot
   (All (A)
        ((Vectorof (Vectorof A))
         Index
         Index
         Index
         (A A -> Boolean)
         ->
         (Values Index A))))
(define (generic-find-partial-pivot rows m i j comp)
  (for/fold ([p : Index i] [pivot (unsafe-vector2d-ref rows i j)])
            ([l (in-range (add1 i) m)]
             #:when (index? l))
    (define new-pivot (unsafe-vector2d-ref rows l j))
    (if (comp new-pivot pivot)
        (values l new-pivot)
        (values p pivot))))

(: make-elim-rows!
   (All (A)
        ((A * -> A) (A -> A)
         (A * -> A) (A -> A)
         (A A -> Boolean)
         ->
         ((Vectorof (Vectorof A)) Index Index Index A Nonnegative-Fixnum -> Void))))
(define (make-elim-rows! +F -F *F /F =F)
  (define 0F (+F))
  (: elim-rows! ((Vectorof (Vectorof A)) Index Index Index A Nonnegative-Fixnum -> Void))
  (define (elim-rows! rows m i j pivot start)
    (define row_i (unsafe-vector-ref rows i))
    (for ([l : Nonnegative-Fixnum (in-range start m)]
          #:unless (fx= l i)
          #:do [(define row_l (unsafe-vector-ref rows l))
                (define x_lj (unsafe-vector-ref row_l j))]
          #:unless (=F x_lj 0F))
      (vector-generic-scaled-add! row_l row_i (-F (*F x_lj (/F pivot))) j +F *F)
      ;; Make sure the element below the pivot is zero
      (unsafe-vector-set! row_l j 0F)))
  elim-rows!)
