#lang typed/racket
(require "../../functions/beta.rkt"
         "../../functions/continued-fraction.rkt"
         "../../../flonum.rkt")

;; *************************************************************************************************
;; ** compute (Beta 1/2 a/2) but more accurate if a is integer?
;; *************************************************************************************************
(provide beta1/2 logbeta1/2)

;; inner valid for POSITIVE integers
(define (beta-inner [n : Integer]) : (values Real Real)
  (if (= n 1) (values 1 1)
      ;; n will be odd, but round to satisfy type-checker
      (let ([n (round (/ (- n 1) 2))])
        (define n0 (floor (/ n 2)))
        (define j0 (if (even? n) (+ n 1) (+ n 2)))
        (define e (expt 2 (+ n n0))) ;; flexp2
        (define stp 1000)
        (let lp : (values Real Real) ([A : Integer 1][B : Integer 1][i : Integer 1][j : Integer j0])
          (if (<= i n0)
              (if (< i stp)
                  (lp (* A i) (* B j) (+ i 1) (+ j 2))
                  (let ([C (/ (* A i) (* B j))])
                    (set! stp (+ stp 1000))
                    (lp (numerator C) (denominator C) (+ i 1) (+ j 2))))
              (values (* A e) B))))))

(: beta1/2 (-> Real Flonum))
(define (beta1/2 a)
  (if (and (integer? a) (< 1 a) (or (exact? a) (< a 111111)))
      ;; if a > 1e6 this starts to take too long
      ;; a ~ 2×a/4 multiplications of exact integers
      (let ([a (exact-round a)])
        (if (even? a)
            (let*-values ([(A) (- a 1)]
                          [(b a) (beta-inner A)])
              (fl (/ (* 2 b) A a)))
            (let-values ([(a b) (beta-inner a)])
              (* pi (fl (/ b a))))))
      (beta 0.5 (fl/ (fl a) 2.))))

(: logbeta1/2 (-> Real Flonum))
(define (logbeta1/2 a)
  (if (and (integer? a) (< 1 a) (or (exact? a) (< a 111111)))
      (let ([a (exact-round a)])
        (if (even? a)
            (let*-values ([(A) (- a 1)]
                          [(b a) (beta-inner A)])
              (fllog (fl (/ (* 2 b) A a))))
            (let-values ([(a b) (beta-inner (exact-round a))])
              (fllog (* pi (fl (/ b a)))))))
      (cond
        [(= a +min.0) 745.1332191019412]
        ;; log of beta works better in following region,
        ;; but for 3 to 15 both are bad, and around 6.7 they are terrible
        [(< 3 a 500) (fllog (beta1/2 a))]
        [else        (fllog-beta 0.5 (fl/ (fl a) 2.))])))

;; *************************************************************************************************
;; ** beta1/2-inc : calculates (* 1/2 (beta-inc ν/2 1/2 x)) REGURALIZED 
;; *************************************************************************************************
(provide flbeta1/2-regularized)

;; copy of "math/private/functions/incomplete-beta.rkt"
(: in-bounds? (Flonum Flonum -> Boolean))
(define (in-bounds? A x)
  (and (x . fl> . 0.0) (x . fl< . 1.0) 
       (A . fl> . 0.0) (A . fl< . +inf.0)))

(: maybe1- (Flonum Any Any -> Flonum))
(define (maybe1- z log? 1-?)
  (cond [1-?  (if log? (lg- (fllog 0.5) z) (fl- 0.5 z))]
        [else  z]))

(: flbeta1/2-regularized-limits (Flonum Flonum Any -> Flonum))
(define (flbeta1/2-regularized-limits A x r?)
  (cond [(or (x . fl< . 0.0) (x . fl> . 1.0)
             (A . fl< . 0.0))
         +nan.0]
        [(fl= x 1.0)     0.5]
        [(fl= A +inf.0)  (if r? 0.5 0.0)]
        [(fl= A 0.0)     (if r? 0.5 0.0)]
        [(fl= x 0.0)     0.0]
        [else           +nan.0]))

;; ------------- flbeta-const -------------------------
(: flbeta-regularized-constA (Flonum Flonum Flonum Flonum Flonum Any -> Flonum))
(define (flbeta-regularized-constA A x y log-x log-y log?); a = 0.5A, b = 0.5
  ;; lanczos never needed b < 1.
  (cond
    [log?      (flsum (list (* 0.5 A log-x) (* 0.5 log-y) (- (logbeta1/2 A))))]
    [(< 2e2 A) (fl/ (flexp (flsum (list (* 0.5 A log-x) (* 0.5 log-y)))) (beta1/2 A))]
    [else      (/ (* (flexpt x (* 0.5 A)) (flsqrt y)) (beta1/2 A))])
  #;(if (or log? (< 3e2 A))
      ((if log? values flexp)
       (flsum (list (* 0.5 A log-x) (* 0.5 log-y) (- (logbeta1/2 A)))))
      (/ (* (flexpt x (* 0.5 A)) (flsqrt y)) (beta1/2 A)))
  #;((if log? values flexp)
   (flsum (list (* 0.5 A log-x) (* 0.5 log-y) (- (logbeta1/2 A))))))

(: flbeta-regularized-constB (Flonum Flonum Flonum Flonum Flonum Any -> Flonum))
(define (flbeta-regularized-constB A x y log-x log-y log?); a = 0.5, b = 0.5A
  ;; lanczos never needed a < 1.
  #;(define log-t (flsum (list (* 0.5 log-x) (* 0.5 A log-y) (- (logbeta1/2 A)))))
  #;(if log? log-t (flexp log-t))
  (cond
    [log?      (flsum (list (* 0.5 log-x) (* 0.5 A log-y) (- (logbeta1/2 A))))]
    [(< 2e2 A) (fl/ (flexp (flsum (list (* 0.5 log-x) (* 0.5 A log-y)))) (beta1/2 A))]
    [else      (/ (* (flsqrt x) (flexpt y (* 0.5 A))) (beta1/2 A))]))

;; ------------- hypergeom -------------------------
(: hypergeom-fac (Flonum Flonum Flonum -> Flonum))
(define (hypergeom-fac a b x)
  (define a+b (fl+ a b))
  (define a+1 (fl+ a 1.0))
  (fl* (fl/ (fl+ a+b 20.0) (fl+ a+1 20.0)) x))

(: hypergeom-facA (Flonum Flonum -> Flonum))
;; Returns the adjustment to the hypergeometric series coefficient at n = 20
;; If this is < 0.5, the series would converge in < 50 iterations or so
(define (hypergeom-facA A x) ; a b
  (fl* (fl/ (fl+ A 41.0) (fl+ A 42.0)) x))
(: hypergeom-facB (Flonum Flonum -> Flonum))
(define (hypergeom-facB A x); b a
  (fl* (fl/ (fl+ A 41.0) 43.0) x))

(: get-hypergeom-params1/2
   (Flonum Flonum Flonum Flonum Flonum Any Any
           -> (Values Flonum Flonum Flonum Flonum Flonum Any Any)))
(define (get-hypergeom-params1/2 A x y log-x log-y 1-? r?)
  (if (if r?
          ((hypergeom-fac (* 0.5 A) 0.5 y) . fl< . (hypergeom-fac 0.5 (* 0.5 A) x))
          ;((hypergeom-facA A y) . fl< . (hypergeom-facB A x))
          
          ;((hypergeom-fac 0.5 (* 0.5 A) y) . fl< . (hypergeom-fac (* 0.5 A) 0.5 x))
          ((hypergeom-facB A y) . fl< . (hypergeom-facA A x)))
      (values A y x log-y log-x (not 1-?) (not r?)) ; b a
      (values A x y log-x log-y 1-? r?))) ; a b

(: flbeta-regularized-hypergeom1/2 (Flonum Flonum Flonum Flonum Flonum Any Any -> Flonum))
;; Computes lower incomplete beta using the hypergeometric series
(define (flbeta-regularized-hypergeom1/2 A x y log-x log-y log? r?) ; r?=#f : a b => a = 0.5A, b=0.5
;(println `(flbeta-regularized-hypergeom1/2 ,A ,x ,y ,log-x ,log-y ,log? ,r?))
  (define a+b/2 (fl+ A 1.0))                                        ; r?=#t : b a => a = 0.5, b= 0.5A
  (define a+2/2 (if r? 3. (fl+ A 2.0)))
  (define: z : Flonum
    (let loop ([z 0.0] [dz 1.0] [2n 0.0] [i -1.0])
      ;(printf "z = ~v  dz = ~v  i = ~v~n" z dz i)
      (define new-dz (fl* (fl* dz (fl/ (fl+ a+b/2 2n)
                                       (fl+ a+2/2 2n)))
                          x)
                     ;no real influence
                     #;(fl* (fl* dz (fl- 1. (fl/ 1. (fl+ a+2/2 2n)))) x)
                     #;(fl* (fl- dz (fl/ dz (fl+ a+2/2 2n))) x))
      (cond [(zero? i)  (fl+ z new-dz)]
            [else
             (let ([i  (if (and (i . fl< . 0.0)
                                ((flabs new-dz) . fl<= . (flabs dz))
                                ((flabs new-dz) . fl<= . (fl* (fl* 0.5 epsilon.0) (flabs z))))
                           3.0
                           (fl- i 1.0))])
               (loop (fl+ z new-dz) new-dz (fl+ 2n 2.0) i))])))
  (define c ((if r? flbeta-regularized-constB flbeta-regularized-constA) A x y log-x log-y log?))
;(println (list a+b/2 a+2/2 (/ a+b/2 a+2/2) z c))
  (cond [log?  (flsum (list c (if r? 0. (- (fllog A))) (fllog1p z)))]
        [r?    (fl+ c (fl* z c))]
        [else  ; (fl/ (fl+ c (fl* z c)) a)
               ; but we want
               ; (fl* 0.5 (fl/ (fl+ c (fl* z c)) a))
               ; = (fl* 0.5 (fl/ (fl+ c (fl* z c)) 0.5 A))
               (fl/ (fl+ c (fl* z c)) A)]))

;; ------------- Frac & largeParam -------------------------
(: get-large-params1/2
   (Flonum Flonum Flonum Flonum Flonum Any Any
           -> (Values Flonum Flonum Flonum Flonum Flonum Flonum Any Any)))
(define (get-large-params1/2 A x y log-x log-y 1-? r?)
  (define l
    (let ([a  (inexact->exact A)]
          [x  (inexact->exact x)])
      (exact->inexact (/ (- (if r? 1 a) (* (+ a 1) x)) 2))))
  (if (fl< l 0.)
      (values A y x log-y log-x (- l) (not 1-?) (not r?))
      (values A x y log-x log-y l 1-? r?)))

(: flbeta-regularized-frac1/2A (Flonum Flonum Flonum Flonum Flonum Flonum Any -> Flonum))
;; Didonato and Morris's continued fraction
(define (flbeta-regularized-frac1/2A A x y log-x log-y l log?); a b => a = 0.5A, b=0.5
  (define a (* 0.5 A))
;(println 'pingA)
  (define-values (s t)
    (continued-fraction-parts
     1.0
     (λ: ([n : Flonum] [s : Flonum])
       (let ([a+n-1  (fl+ a n -1.)]
             [a+2n-1  (fl+ a (fl* 2.0 n) -1.0)]
             [a+b+n-1  (fl+ a n -.5)])
         (fl/ (fl* (fl* (fl* (fl* (fl* a+n-1 a+b+n-1) n) (fl- 0.5 n)) x) x)
              (fl* a+2n-1 a+2n-1))))
     (fl* (fl/ a (fl+ a 1.0)) (fl+ l 1.0))
     (λ: ([n : Flonum] [t : Flonum])
       (let ([a+2n  (fl+ a (fl* 2.0 n))])
         (fl+ (fl+ (fl/ (fl* (fl* n (fl- 0.5 n)) x)
                        (fl+ a+2n -1.0))
                   (fl* (fl/ (fl+ a n) (fl+ a+2n 1.0))
                        (fl+ (fl+ l 1.0) (fl* n (fl+ 1.0 y)))))
              n)))
     (fl* 0.5 epsilon.0)))
  (define c (flbeta-regularized-constA A x y log-x log-y log?))
  (cond [log?  (fl+ (fllog 0.5) c (fllog-quotient s t))]
        [else  (fl* 0.5 c (fl/ s t))]))

(: flbeta-regularized-frac1/2B (Flonum Flonum Flonum Flonum Flonum Flonum Any -> Flonum))
;; Didonato and Morris's continued fraction
(define (flbeta-regularized-frac1/2B A x y log-x log-y l log?) ; r?=#f : a b => a = 0.5A, b=0.5
  (define b (* 0.5 A))
;(println 'pingB)
  (define-values (s t)
    (continued-fraction-parts
     1.0
     (λ: ([n : Flonum] [s : Flonum])
       (let ([a+n-1  (fl+ n -.5)]
             [a+2n-1  (fl+ (fl* 2.0 n) -.5)]
             [a+b+n-1  (fl+ b n -.5)])
         (fl/ (fl* (fl* (fl* (fl* (fl* a+n-1 a+b+n-1) n) (fl- b n)) x) x)
              (fl* a+2n-1 a+2n-1))))
     (fl* #i1/3 (fl+ l 1.0))
     (λ: ([n : Flonum] [t : Flonum])
       (let ([a+2n  (fl+ 0.5 (fl* 2.0 n))])
         (fl+ (fl+ (fl/ (fl* (fl* n (fl- b n)) x)
                        (fl+ a+2n -1.0))
                   (fl* (fl/ (fl+ 0.5 n) (fl+ a+2n 1.0))
                        (fl+ (fl+ l 1.0) (fl* n (fl+ 1.0 y)))))
              n)))
     (fl* 0.5 epsilon.0)))
  (define c (flbeta-regularized-constA A x y log-x log-y log?))
  (cond [log?  (fl+ (fllog 0.5) c (fllog-quotient s t))]
        [else  (fl* 0.5 c (fl/ s t))]))

;; ------------- all together now -------------------------
(: flbeta1/2-regularized (Flonum Flonum Flonum Flonum Any Any Any -> Flonum))
(define (flbeta1/2-regularized A X x y log? 1-? r?)
;(println `(flbeta1/2-regularized ,A ,x ,y ,log? ,1-? ,r?))
  ; r?=#f : a b => a = 0.5A, b=0.5
  ; r?=#t : b a => a = 0.5, b= 0.5A
  (define log-x (if (< x 0.5) (fllog x) (fllog1p (- y))))
  (define log-y (if (< y 0.5) (fllog y) (fllog1p (- x))))
  (cond
    [(not (in-bounds? A x))
     (define z (flbeta1/2-regularized-limits A x r?))
     (maybe1- (if log? (fllog z) z)
              log? 1-?)]
    [(A . fl< . 2.0) ;; (and (a < 1.) (b < 1.))
     (let-values ([(A x y log-x log-y 1-? r?)
                   (get-hypergeom-params1/2 A x y log-x log-y 1-? r?)])
       (maybe1- (flbeta-regularized-hypergeom1/2 A x y log-x log-y log? r?)
                log? 1-?))]
    [(and (< -1.6 X) ;; works well for cfd4 in range x -1e3 .. -1e0 for ν 1e0 .. 1e20
          (((if r? hypergeom-facB hypergeom-facA) A x) . fl< . 0.75))
     (maybe1- (flbeta-regularized-hypergeom1/2 A x y log-x log-y log? r?)
              log? 1-?)]
    [(and (< X -1e2) ;; works well for cfd4 in range x -1e3 .. -1e0 for ν 1e0 .. 1e20
          (((if r? hypergeom-facA hypergeom-facB) A y) . fl< . 0.75))
     (maybe1- (flbeta-regularized-hypergeom1/2 A y x log-y log-x log? (not r?))
              log? (not 1-?))]
    [else
     (let-values ([(A x y log-x log-y l 1-? r?)
                   (get-large-params1/2 A x y log-x log-y 1-? r?)])
       ; assym branch can never be #t since only valid if both a && b > 100
       (define z ((if r? flbeta-regularized-frac1/2B flbeta-regularized-frac1/2A) A x y log-x log-y l log?))
       (maybe1- z log? 1-?))]))

