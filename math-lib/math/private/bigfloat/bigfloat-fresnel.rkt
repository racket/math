#lang typed/racket/base

;Bigfloat implementation:
;  adaptation of powerseries as shown on wikipedia
;  this implementation is potentially exact
;  but for large x (> 5!) it is really slow and needs a lot of bits in bf-precision (~2x²)!


(require "bigfloat-struct.rkt")

(provide bfFresnel-S bfFresnel-RS bfFresnel-C bfFresnel-RC)

(define (precision-check [a : Bigfloat][maxp : Integer]) : Integer
  (define p (bf-precision))
  (define a2 (expt (bigfloat->real a) 2))
  (define min-precision-needed (* 2 a2))
  (define expected-precision-loss (* 4/5 a2))
  (define mp (+ 5 (round (inexact->exact (max min-precision-needed (+ p expected-precision-loss))))))
  (cond
    [(<= mp maxp) mp]
    [else
     (error (format "bfFresnel: calculation aborted
 Minimum precision needed for calculating ~a... is ~a
 This is more than the maximum allowed calculating precision ~a->~a."
                    (bigfloat->flonum a) mp p maxp))]))
(define (bfFresnel-RS [x : Bigfloat][maxp : Integer (* 2 (bf-precision))])
  (define p (bf-precision))
  (bfcopy
   (parameterize ([bf-precision (precision-check x maxp)])
     (define X3 (bfexpt x (bf 3)))
     (define X4 (bfexpt x (bf 4)))
     (define prsn (bfexpt (bf 1/2) (bf p)))
     (define-values (s l)
       (for/fold : (Values Bigfloat Bigfloat)
         ([s : Bigfloat (bf/ X3 (bf 3))]
          [l : Bigfloat (bf/ X3 (bf 3))])
         ([n (in-naturals 1)]
          #:break (bf< (bfabs (bf/ l s)) prsn))
         (define l+ (bf* (bf -1) X4 (bf* l (bf (/ (- (* 4 n) 1)
                                                  (* 2 n)(+ (* 2 n) 1)(+ (* 4 n) 3))))))
         (values (bf+ s l+) l+)))
     s)))
(define (bfFresnel-S [x : Bigfloat][maxp : Integer (* 2 (bf-precision))])
  (bf* (bfsqrt (bf/ (bf 2) pi.bf))
       (bfFresnel-RS (bf* (bfsqrt (bf/ pi.bf (bf 2))) x) maxp)))
  
(define (bfFresnel-RC [x : Bigfloat][maxp : Integer (* 2 (bf-precision))])
  (define p (bf-precision))
  (bfcopy
   (parameterize ([bf-precision (precision-check x maxp)])
     (define X4 (bfexpt x (bf 4)))
     (define prsn (bfexpt (bf 1/2) (bf (bf-precision))))
     (define-values (s l)
       (for/fold : (Values Bigfloat Bigfloat)
         ([s : Bigfloat x]
          [l : Bigfloat x])
         ([n (in-naturals 1)]
          #:break (bf< (bfabs (bf/ l s)) prsn))
         (define l+ (bf* (bf -1) X4 (bf* l (bf (/ (- (* 4 n) 3)
                                                  (* 2 n)(- (* 2 n) 1)(+ (* 4 n) 1))))))
         (values (bf+ s l+) l+)))
     s)))
(define (bfFresnel-C [x : Bigfloat][maxp : Integer (* 2 (bf-precision))])
  (bf* (bfsqrt (bf/ (bf 2) pi.bf))
       (bfFresnel-RC (bf* (bfsqrt (bf/ pi.bf (bf 2))) x) maxp)))
