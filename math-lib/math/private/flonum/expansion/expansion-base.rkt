#lang typed/racket/base

#|
Arithmetic based on:

Jonathan Richard Shewchuk
Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates
Discrete & Computational Geometry 18(3):305–363, October 1997

and

Tight and rigourous error bounds for basic building blocks of double-word arithmetic
Mioara Joldes, Jean-Michel Muller, Valentina Popescu
ACM Transactions on Mathematical Software, Association for Computing Machinery, 2017, 44 (2)
https://hal.archives-ouvertes.fr/hal-01351529v3

Note: Joldes et al provides algorithms that do not consider
over/underflow, and thus some adaptations have been made for that. In
most cases, that simply checks for overflow in intermediate
computations and then produces a "simple" answer which is the
one-double approximation using just the high-components and 0.0 as the
second component.

Note: better algorithms exist in crlibm but are under the LGPL so
they're not used here.

|#

(require racket/math
         "../flonum-functions.rkt"
         "../flonum-bits.rkt"
         "../flonum-error.rkt"
         (only-in "../flonum-error.rkt"
                  ;; These names are used in Joldes et al.
                  ;; One would hope that using the `fast-` variants
                  ;; would be possible here, but it results in
                  ;; inaccuracy/overflow.
                  [fl+/error 2Sum]
                  [fl+/error Fast2Sum]
                  [fl*/error 2Prod])
         "../flonum-constants.rkt"
         "../utils.rkt")

(provide fl2? fl2zero? fl2rational? fl2positive? fl2negative? fl2infinite? fl2nan?
         fl2 fl2->real
         fl2ulp fl2ulp-error fl2step fl2next fl2prev
         +max.hi +max.lo -max.hi -max.lo
         +max-subnormal.hi -max-subnormal.hi
         fl2abs fl2+ fl2- fl2= fl2> fl2< fl2>= fl2<=
         fl2*split-fl fl2* fl2sqr fl2/
         fl2sqrt flsqrt/error)

(: floverlapping? (Flonum Flonum -> Boolean))
(define (floverlapping? x2 x1)
  (define-values (s2 e2) (flonum->sig+exp (flabs x2)))
  (define-values (s1 e1) (flonum->sig+exp (flabs x1)))
  (define-values (n1 n2)
    (if (e2 . > . e1)
        (values s1 (arithmetic-shift s2 (- e2 e1)))
        (values (arithmetic-shift s1 (- e1 e2)) s2)))
  (not (= (bitwise-ior n1 n2)
          (bitwise-xor n1 n2))))

(: fl2? (Flonum Flonum -> Boolean))
(define (fl2? x2 x1)
  (cond [(flrational? x2)
         (cond [(flrational? x1)
                (cond [((flabs x2) . < . (flabs x1))  #f]
                      [else (not (floverlapping? x2 x1))])]
               [else  #f])]
        [else
         (fl= x1 0.0)]))

(define-syntax-rule (define-simple-fl2-predicate fl2pred? flpred?)
  (begin
    (: fl2pred? (Flonum Flonum -> Boolean))
    (define (fl2pred? x2 x1)
      (flpred? (fl+ x2 x1)))))

(define-simple-fl2-predicate fl2zero? (λ (x) (fl= x 0.0)))
(define-simple-fl2-predicate fl2positive? (λ (x) (fl> x 0.0)))
(define-simple-fl2-predicate fl2negative? (λ (x) (fl< x 0.0)))
(define-simple-fl2-predicate fl2rational? flrational?)
(define-simple-fl2-predicate fl2nan? flnan?)
(define-simple-fl2-predicate fl2infinite? flinfinite?)

;; ===================================================================================================
;; Conversion

(: fl2 (case-> (Real -> (Values Flonum Flonum))
               (Flonum Flonum -> (Values Flonum Flonum))))
(define fl2
  (case-lambda
    [(x)
     (cond [(flonum? x)  (values x 0.0)]
           [(single-flonum? x)  (values (fl x) 0.0)]
           [else
            (define x2 (fl x))
            (if (flinfinite? x2)
                (values x2 0.0)
                (let* ([x  (- x (inexact->exact x2))]
                       [x1  (fl x)]
                       [x  (- x (inexact->exact x1))])
                  (let-values ([(x2 x1)  (fl+/error x2 x1)])
                    (values x2 (fl+ x1 (fl x))))))])]
    [(x2 x1)
     (if (and (fl= x2 0.0) (fl= x1 0.0))
         (values x2 0.0)
         (fl+/error x2 x1))]))

(: fl2eqv? (case-> (Flonum Flonum Flonum -> Boolean)
                   (Flonum Flonum Flonum Flonum -> Boolean)))
(define (fl2eqv? x2 x1 y2 [y1 0.0])
  (and (eqv? x2 y2) (fl= x1 y1)))

(: fl2->real (Flonum Flonum -> Real))
(define (fl2->real x2 x1)
  (if (flrational? x2)
      (+ (inexact->exact x2) (inexact->exact x1))
      x2))

(: fl4->fl2 (Flonum Flonum Flonum Flonum -> (Values Flonum Flonum)))
(define (fl4->fl2 e4 e3 e2 e1)
  (values e4 (fl+ e3 (fl+ e2 e1))))

;; ===================================================================================================
;; Error

(: fl2ulp (Flonum Flonum -> Flonum))
(define (fl2ulp x2 x1)
  (cond [(fl= x2 0.0)  0.0]
        [else  (flmax +min.0 (fl* (flulp x2) epsilon.0))]))

(: fl2ulp-error (Flonum Flonum Real -> Flonum))
(define (fl2ulp-error x2 x1 r)
  (define x (fl2->real x2 x1))
  (define-values (r2 r1) (fl2 r))
  (cond [(eqv? x r)  0.0]
        [(and (fl= x2 0.0) (fl= r2 0.0))  0.0]
        [(and (fl= x2 +inf.0) (fl= r2 +inf.0))  0.0]
        [(and (fl= x2 -inf.0) (fl= r2 -inf.0))  0.0]
        [(zero? r)  +inf.0]
        [(and (rational? x) (flrational? r2))
         (flabs (fl (/ (- (inexact->exact x) (inexact->exact r))
                       (inexact->exact (flmax +min.0 (fl2ulp r2 r1))))))]
        [else  +inf.0]))

(define-values (+max.hi +max.lo)
  (values +max.0 (flprev (* 0.5 (flulp +max.0)))))

(define-values (-max.hi -max.lo)
  (values (- +max.hi) (- +max.lo)))

(: fl2step (Flonum Flonum Integer -> (Values Flonum Flonum)))
(define (fl2step x2 x1 n)
  (let-values ([(x2 x1)  (fast-fl+/error x2 x1)])
    (cond [(flnan? x2)  (values +nan.0 0.0)]
          [(fl= x2 +inf.0)  (fl+/error +max.hi (flstep +max.lo (+ n 1)))]
          [(fl= x2 -inf.0)  (fl+/error -max.hi (flstep -max.lo (- n 1)))]
          [else  (fl+/error x2 (flstep x1 n))])))

(: fl2next (Flonum Flonum -> (Values Flonum Flonum)))
(define (fl2next x2 x1) (fl2step x2 x1 1))

(: fl2prev (Flonum Flonum -> (Values Flonum Flonum)))
(define (fl2prev x2 x1) (fl2step x2 x1 -1))

(define +min-normal.hi (fl/ (flnext +max-subnormal.0) epsilon.0))

(define-values (+max-subnormal.hi +max-subnormal.lo)
  (fl2prev +min-normal.hi 0.0))

(define-values (-max-subnormal.hi -max-subnormal.lo)
  (values (- +max-subnormal.hi) (- +max-subnormal.lo)))

;; ===================================================================================================
;; Absolute value

(: fl2abs (case-> (Flonum -> (Values Flonum Flonum))
                  (Flonum Flonum -> (Values Flonum Flonum))))
(define fl2abs
  (case-lambda
    [(x)  (values (flabs x) 0.0)]
    [(x2 x1)
     (cond [(flnan? x2)  (values +nan.0 0.0)]
           [(fl= x2 0.0)  (values 0.0 0.0)]
           [(fl> x2 0.0)  (values x2 x1)]
           [else  (values (- x2) (- x1))])]))

;; ===================================================================================================
;; Addition and subtraction

(: fl3->fl2 (Flonum Flonum Flonum -> (Values Flonum Flonum)))
(define (fl3->fl2 e3 e2 e1)
  (values e3 (fl+ e2 e1)))

(: raw-fl2+fl (Flonum Flonum Flonum -> (Values Flonum Flonum Flonum)))
(define (raw-fl2+fl e2 e1 b)
  (let*-values ([(Q h1)  (fast-fl+/error b e1)]
                [(h3 h2)  (fast-fl+/error Q e2)])
    (values h3 h2 h1)))

(: raw-fl2+ (Flonum Flonum Flonum Flonum -> (Values Flonum Flonum Flonum Flonum)))
(define (raw-fl2+ e2 e1 f2 f1)
  (let*-values ([(h3 h2 h1)  (raw-fl2+fl e2 e1 f1)]
                [(h4 h3 h2)  (raw-fl2+fl h3 h2 f2)])
    (values h4 h3 h2 h1)))

(: fl2+ (case-> (Flonum Flonum Flonum -> (Values Flonum Flonum))
                (Flonum Flonum Flonum Flonum -> (Values Flonum Flonum))))
;; Algorithm 6 from Joldes et al
(define (fl2+ xh xl yh [yl 0.0])
  (define r (fl+ xh yh))
  (cond
    ;; bail out in weird cases:
    [(not (flrational? r)) (values r 0.0)]
    [(and (fl= xh 0.0) (fl= yh 0.0)) (values r 0.0)]
    [else
     (let*-values ([(sh sl) (2Sum xh yh)]
                   [(th tl) (2Sum xl yl)]
                   [(c) (fl+ sl th)]
                   [(vh vl) (Fast2Sum sh c)]
                   [(w) (fl+ tl vl)]
                   [(zh zl) (Fast2Sum vh w)])
       (values zh zl))]))

(: fl2- (case-> (Flonum Flonum Flonum -> (Values Flonum Flonum))
                (Flonum Flonum Flonum Flonum -> (Values Flonum Flonum))))
(define (fl2- x2 x1 y2 [y1 0.0])
  (fl2+ x2 x1 (- y2) (- y1)))

;; ===================================================================================================
;; Comparison

(define-syntax-rule (define-fl2-comparison name flcomp)
  (begin
    (: name (Flonum Flonum Flonum Flonum -> Boolean))
    (define (name x2 x1 y2 y1)
      (let-values ([(z2 z1)  (fl2- x2 x1 y2 y1)])
        ((fl+ z2 z1) . flcomp . 0.0)))))

(define-fl2-comparison fl2= fl=)
(define-fl2-comparison fl2> fl>)
(define-fl2-comparison fl2< fl<)
(define-fl2-comparison fl2>= fl>=)
(define-fl2-comparison fl2<= fl<=)

;; ===================================================================================================
;; Multiplication and division

(: raw-split-fl2*split-fl (Flonum Flonum Flonum Flonum Flonum Flonum
                                  -> (Values Flonum Flonum Flonum Flonum)))
(define (raw-split-fl2*split-fl e2-hi e2-lo e1-hi e1-lo b-hi b-lo)
  (let*-values ([(b)   (fl+ b-lo b-hi)]
                [(Q1)  (fl* (fl+ e1-lo e1-hi) b)]
                [(h1)  (- (- Q1
                             (fl* e1-hi b-hi)
                             (fl* e1-lo b-hi)
                             (fl* e1-hi b-lo)
                             (fl* e1-lo b-lo)))]
                [(T)  (fl* (fl+ e2-lo e2-hi) b)]
                [(t)  (- (- T
                            (fl* e2-hi b-hi)
                            (fl* e2-lo b-hi)
                            (fl* e2-hi b-lo)
                            (fl* e2-lo b-lo)))]
                [(Q2 h2)  (fast-fl+/error Q1 t)]
                [(h4 h3)  (fast-mono-fl+/error T Q2)])
    (values h4 h3 h2 h1)))

(: split-fl2*split-fl (Flonum Flonum Flonum Flonum Flonum Flonum -> (Values Flonum Flonum)))
(define (split-fl2*split-fl e2-hi e2-lo e1-hi e1-lo b-hi b-lo)
  (let-values ([(h4 h3 h2 h1)  (raw-split-fl2*split-fl e2-hi e2-lo e1-hi e1-lo b-hi b-lo)])
    (fl4->fl2 h4 h3 h2 h1)))

(: fl2*split-fl (Flonum Flonum Flonum Flonum -> (Values Flonum Flonum)))
(define (fl2*split-fl e2 e1 b-hi b-lo)
  (let*-values ([(e2-hi e2-lo)  (flsplit e2)]
                [(e1-hi e1-lo)  (flsplit e1)]
                [(h4 h3 h2 h1)  (raw-split-fl2*split-fl e2-hi e2-lo e1-hi e1-lo b-hi b-lo)])
    (fl4->fl2 h4 h3 h2 h1)))

(: fl2* (case-> (Flonum Flonum Flonum -> (Values Flonum Flonum))
                (Flonum Flonum Flonum Flonum -> (Values Flonum Flonum))))
(define (fl2* xh xl yh [yl 0.0])
  (define z (fl* xh yh))
  (cond [(fl= 0.0 z) (values z 0.0)]
        [(flsubnormal? z) (values z 0.0)]
        [(and (flrational? xh) (flrational? yh) (z . fl> . -inf.0) (z . fl< . +inf.0))
         ;; Algorithm 10 from Joldes et al
         (let*-values ([(ch cl1) (fl*/error xh yh)]
                       [(tl1) (fl* xh yl)]
                       [(tl2) (fl* xl yh)]
                       [(cl2) (fl+ tl1 tl2)]
                       [(cl3) (fl+ cl1 cl2)])
           (if (or (flinfinite? tl1) (flinfinite? tl2))
               (values z 0.0)
               (Fast2Sum ch cl3)))]
        [else (values z 0.0)]))

(: fl2sqr (case-> (Flonum -> (Values Flonum Flonum))
                  (Flonum Flonum -> (Values Flonum Flonum))))
(define fl2sqr
  (case-lambda
    [(x)  (flsqr/error x)]
    [(x2 x1) (fl2* x2 x1 x2 x1)]))

(: fl2/ (case-> (Flonum Flonum Flonum -> (Values Flonum Flonum))
                (Flonum Flonum Flonum Flonum -> (Values Flonum Flonum))))
(define fl2/
  (lambda
    (xh xl yh [yl 0.0])
     ;; Algorithm 17 from Joldes et al
    (define z (fl/ xh yh))
    (cond [(and (flrational? z) (not (fl= z 0.0)) (flrational? yh))
           (let*-values ([(th) (fl/ xh yh)]
                         [(rh rl) (fl2* yh yl th)]
                         [(pih) (fl- xh rh)]
                         [(dl) (fl- xl rl)]
                         [(d) (fl+ pih dl)]
                         [(tl) (fl/ d yh)])
             (if (and (flrational? rh) (flrational? rl) (flrational? tl))
                 (Fast2Sum th tl)
                 (values th 0.0)))]
          [else (values z 0.0)])))

;; ===================================================================================================
;; Square roots

;; One-flonum estimate followed by one Newton's method iteration
;; This could be a little faster if `y' were split only once
(: flsqrt/error (Flonum -> (Values Flonum Flonum)))
(define (flsqrt/error x)
  (let*-values ([(y)  (flsqrt x)]
                [(z2 z1)  (flsqr/error y)]
                [(dy2 dy1)  (fl2+ (- z2) (- z1) x)]
                [(dy2 dy1)  (fl2/ dy2 dy1 y)])
    (if (and (flrational? dy2) (flrational? dy1))
        (fl2+ (* 0.5 dy2) (* 0.5 dy1) y)
        (values y 0.0))))

(: fl2sqrt (case-> (Flonum -> (Values Flonum Flonum))
                   (Flonum Flonum -> (Values Flonum Flonum))))
(define (fl2sqrt x2 [x1 0.0])
  (let*-values ([(y)  (flsqrt (fl+ x1 x2))]
                [(z2 z1)  (flsqr/error y)]
                [(dy2 dy1)  (fl2- x2 x1 z2 z1)]
                [(dy2 dy1)  (fl2/ dy2 dy1 y)]
                [(r2 r1) (fl2+ (* 0.5 dy2) (* 0.5 dy1) y)])
    (if (and (flrational? dy2) (flrational? dy1))
        (values r2 r1)
        (values (flsqrt x2) 0.0))))
