#lang typed/racket/base

(require typed/racket/unsafe
         racket/fixnum
         "matrix-types.rkt"
         "base/utils.rkt"
         "algebra/utils.rkt"
         "../../base.rkt")

(provide (all-from-out "base/utils.rkt")
         (all-defined-out))

(: find-partial-pivot
   (case-> ((Vectorof (Vectorof Flonum)) Index Index Index -> (Values Index Flonum))
           ((Vectorof (Vectorof Real)) Index Index Index -> (Values Index Real))
           ((Vectorof (Vectorof Float-Complex)) Index Index Index -> (Values Index Float-Complex))
           ((Vectorof (Vectorof Number)) Index Index Index -> (Values Index Number))))
;; Find the element with maximum magnitude in a column
(define (find-partial-pivot rows m i j)
  (generic-find-partial-pivot rows m i j magnitude-compare))

(: find-first-pivot
   (case-> ((Vectorof (Vectorof Flonum)) Index Index Index -> (Values Index Flonum))
           ((Vectorof (Vectorof Real)) Index Index Index -> (Values Index Real))
           ((Vectorof (Vectorof Float-Complex)) Index Index Index -> (Values Index Float-Complex))
           ((Vectorof (Vectorof Number)) Index Index Index -> (Values Index Number))))
;; Find the first nonzero element in a column
(define (find-first-pivot rows m i j)
  (define pivot (unsafe-vector2d-ref rows i j))
  (if ((magnitude pivot) . > . 0)
      (values i pivot)
      (let loop ([#{l : Nonnegative-Fixnum} (fx+ i 1)])
        (cond [(l . fx< . m)
               (define pivot (unsafe-vector2d-ref rows l j))
               (if ((magnitude pivot) . > . 0) (values l pivot) (loop (fx+ l 1)))]
              [else
               (values i pivot)]))))

(: elim-rows!
   (case-> ((Vectorof (Vectorof Flonum)) Index Index Index Flonum Nonnegative-Fixnum -> Void)
           ((Vectorof (Vectorof Real)) Index Index Index Real Nonnegative-Fixnum -> Void)
           ((Vectorof (Vectorof Float-Complex)) Index Index Index Float-Complex Nonnegative-Fixnum
                                                -> Void)
           ((Vectorof (Vectorof Number)) Index Index Index Number Nonnegative-Fixnum -> Void)))
(define (elim-rows! rows m i j pivot start)
  (define x00 (unsafe-vector2d-ref rows 0 0))
  (define +F (add x00))
  (define -F (sub x00))
  (define *F (mul x00))
  (define /F (div x00))
  (define =F (eqv x00))
  ((make-elim-rows! +F -F *F /F =F) rows m i j pivot start))


(: magnitude-compare
   (case-> (Flonum Flonum -> Boolean)
           (Real Real -> Boolean)
           (Float-Complex Float-Complex -> Boolean)
           (Number Number -> Boolean)))
(define (magnitude-compare a b) (> (magnitude a) (magnitude b)))

(: one (case-> (Flonum -> Nonnegative-Flonum)
               (Real -> (U 1 Nonnegative-Flonum))
               (Float-Complex -> Nonnegative-Flonum)
               (Number -> (U 1 Nonnegative-Flonum))))
(define (one x)
  (cond [(flonum? x)  1.0]
        [(real? x)  1]
        [(float-complex? x)  1.0]
        [else  1]))

(: zero (case-> (Flonum -> Flonum-Positive-Zero)
                (Real -> (U 0 Flonum-Positive-Zero))
                (Float-Complex -> Flonum-Positive-Zero)
                (Number -> (U 0 Flonum-Positive-Zero))))
(define (zero x)
  (cond [(flonum? x)  0.0]
        [(real? x)  0]
        [(float-complex? x)  0.0]
        [else  0]))

(: one* (case-> (Flonum -> Nonnegative-Flonum)
                (Real -> (U 1 Nonnegative-Flonum))
                (Float-Complex -> Float-Complex)
                (Number -> (U 1 Nonnegative-Flonum Float-Complex))))
(define (one* x)
  (cond [(flonum? x)  1.0]
        [(real? x)  1]
        [(float-complex? x)  1.0+0.0i]
        [else  1]))

(: zero* (case-> (Flonum -> Flonum-Positive-Zero)
                 (Real -> (U 0 Flonum-Positive-Zero))
                 (Float-Complex -> Float-Complex)
                 (Number -> (U 0 Flonum-Positive-Zero Float-Complex))))
(define (zero* x)
  (cond [(flonum? x)  0.0]
        [(real? x)  0]
        [(float-complex? x)  0.0+0.0i]
        [else  0]))

(unsafe-require/typed/provide "untyped-utils.rkt"
  [eqv (case-> (Flonum -> (Flonum Flonum -> Boolean))
               (Real -> (Real Real -> Boolean))
               (Float-Complex -> (Float-Complex Float-Complex -> Boolean))
               (Number -> (Number Number -> Boolean)))]
  [sub (case-> (Flonum -> (Flonum * -> Flonum))
               (Real -> (Real * -> Real))
               (Float-Complex -> (Float-Complex * -> Float-Complex))
               (Number -> (Number * -> Number)))]
  [div (case-> (Flonum -> (Flonum * -> Flonum))
               (Real -> (Real * -> Real))
               (Float-Complex -> (Float-Complex * -> Float-Complex))
               (Number -> (Number * -> Number)))]
  [add (case-> (Flonum -> (Flonum * -> Flonum))
               (Real -> (Real * -> Real))
               (Float-Complex -> (Float-Complex * -> Float-Complex))
               (Number -> (Number * -> Number)))]
  [mul (case-> (Flonum -> (Flonum * -> Flonum))
               (Real -> (Real * -> Real))
               (Float-Complex -> (Float-Complex * -> Float-Complex))
               (Number -> (Number * -> Number)))]
  [add* (case-> (Flonum -> (Flonum * -> Flonum))
                (Real -> (Real * -> Real))
                (Float-Complex -> (Float-Complex * -> Float-Complex))
                (Number -> (Number * -> Number)))]
  [mul* (case-> (Flonum -> (Flonum * -> Flonum))
                (Real -> (Real * -> Real))
                (Float-Complex -> (Float-Complex * -> Float-Complex))
                (Number -> (Number * -> Number)))])
