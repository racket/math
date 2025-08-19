#lang typed/racket
(require racket/performance-hint
         (only-in racket/math nan?)
         "../../base/base-constants.rkt"
         "../../../flonum.rkt"
         "../dist-struct.rkt"
         "delta-dist.rkt")

(provide make-weibull-pdf
         make-weibull-cdf
         make-weibull-inverse-cdf
         make-weibull-sampler)
                     
;;; Weibull distribution

;; Parameter               - Domain
;   k   - shape            - any positive real number
;   d   - location         - any real number
;   s   - scale parameter  - any real number


;;;
;;; PDF (->* (Real) (Any) Flonum))
;;;

(: make-std-nrm-weibull-pdf (-> Flonum (-> Flonum Flonum)))
(define (make-std-nrm-weibull-pdf k)
  (cond
    [(fl< 3e19 k)
     (define k/e (/ k euler.0))
     (λ (x) (if (= x 1.) k/e 0.))]
    [(fl< 0. k)
     (begin-encourage-inline
       (define (pdf33 [x : Flonum]) : Flonum
         ;; x^k => 0 == (exp (- x^k)) => 1
         (define x^k/2 (flexpt x (fl/ (fl- k 1.) 2.)))
         (fl* x^k/2 k x^k/2))
       (define (pdf34 [x : Flonum]) : Flonum
         ;; x^k => 0 == (exp (- x^k)) => 1
         ;; but don't use -1 for large exponents
         (define x^k/2 (flexpt x (fl/ k 2.)))
         (fl* x^k/2 (fl/ k x) x^k/2))
       (define (pdf411 [x : Flonum]) : Flonum
         (define x^k (flexpt x k))
         (fl* (fl/ k x) x^k (flexp (fl- x^k))))
       (define (pdf412 [x : Flonum]) : Flonum
         (define x^k (flexpt x k))
         (fl* k (flexpt x (fl- k 1.)) (flexp (fl- x^k))))
       (define (pdf413 [x : Flonum]) : Flonum
         (define x^k (flexpt x k))
         (/ (fl* k x^k (flexp (fl- x^k))) x))
       (define (pdf43 [x : Flonum]) : Flonum
         (define-values (lh ll) (fl2log x 0.))
         (define-values (xh xl) (fl2* lh ll k))
         (define-values (x^k x^kl) (fl2exp xh xl))
         (define-values (e^h e^l) (fl2exp (fl- x^k) (fl- x^kl)))
         (fl* x^k (fl/ k x) e^h))
       (define (pdf437 [x : Flonum]) : Flonum
         (define-values (lh ll) (fl2log x 0.))
         (define-values (xh xl) (fl2* lh ll k))
         (define-values (x^k x^kl) (fl2exp xh xl))
         (define-values (e^h e^l) (fl2exp (fl* -0.5 x^k) (fl* -0.5 x^kl)))
         (define T (fl* (fl/ k x) x^k))
         (fl+ (fl* T e^h e^h)
              (fl* 2. T e^h e^l)
              (fl* T e^l e^l))))

     (define smallX : (-> Flonum Flonum) ;; x < 1
       (cond
         [(<= k (* 1e3 +max-subnormal.0))
          pdf411]
         [(< k 1.)
          (λ (x)
            (cond
              [(< x +max-subnormal.0)
               (pdf413 x)]
              [else
               (pdf411 x)]))]
         [(< k 3.)
          pdf412]
         [else
          (define LIMIT00 (flexpt +max-subnormal.0 (fl/ 0.99 (fl- k 1.))))
          (define LIMIT1 (flexpt +min.0 (fl/ (max 1. (fl+ 1. (fl/ (fllog k) 200.))) (fl- k 1.))))
          (cond
            [(< k 8.e15)
             (define LIMIT0 (if (< k 50.)
                                (fl* 1.5 (flexpt +max-subnormal.0 (fl/ 1.00 (fl- k 1.))))
                                LIMIT00))
             (λ (x)
               (cond
                 [(< LIMIT0 x) (pdf412 x)]
                 [(< LIMIT1 x) (pdf33 x)]
                 [else          0.]))]
            [else
             (λ (x)
               (cond
                 [(< LIMIT00 x) (pdf411 x)]
                 [(< LIMIT1 x)  (pdf34 x)]
                 [else           0.]))])]))
     (define bigX : (-> Flonum Flonum) ;; 1 < x
       (cond
         [(< k 1e-4) pdf411]
         [else
          (define LIMIT0 (flexpt 3.5 (/ k)))
          (define LIMIT1 (flexpt 708. (/ k)))
          (define LIMIT2 (flexpt (fl* 750. (max 1. (fl+ 1. (fl/ (fllog k) 100.)))) (/ k)))
          (λ (x)
            (cond
              [(< x LIMIT0) (pdf411 x)]
              [(< x LIMIT1) (pdf43 x)]
              [(< x LIMIT2) (pdf437 x)]
              [else         0.]))]))
     (λ (x)
       (cond
         [(< x 0.) 0.]
         [(= x 1.) (/ k euler.0)]
         [(< x 1.) (smallX x)]
         [else     (bigX x)]))]
    [else (λ (x) +nan.0)]))

(: make-std-log-weibull-pdf (-> Flonum (-> Flonum Flonum)))
(define (make-std-log-weibull-pdf k)
  (cond
    [(< 0. k)
     (define lnK (fllog k))
     (define nrm-pdf (make-std-nrm-weibull-pdf k))
     (begin-encourage-inline
       (define (lpdf0 [x : Flonum]) : Flonum
         (fllog (nrm-pdf x)))
       (define (lpdf1 [x : Flonum]) : Flonum
         (define X (fllog x))
         (flsum (list lnK (fl* X (fl- k 1.)) (fl- (flexpt x k)))))
       (define (lpdf6 [x : Flonum]) : Flonum
         (define-values (Xh Xl) (fl2log x 0.))
         (define-values (Kh Kl) (fl2log k 0.))
         (define-values (kXh kXl) (fl2* Xh Xl k))
         (define-values (xkh xkl) (fl2exp kXh kXl))
         (flsum (list Kh Kl (fl- Xh) (fl- Xl) kXh kXl (fl- xkh) (fl- xkl)))))

     ; (fllog +max.0) -> 710.
     (define xMax (flexp (/ (fllog +max.0) k)))

     (define bigX : (-> Flonum Flonum) ;; 1 < x
       (cond
         [(<= k 1.15)
          (λ (x)
            (cond
              [(< x xMax) (lpdf1 x)]
              [else        -inf.0]))]
         [else
          (define xx (flexpt 1e280 (fl/ -1. (- k 1.))))
          (define XX (min 1. (- (* .28 k) .312)))

          (cond
            [(<= k 1.5)
             (λ (x)
               (cond
                 [(< xx x XX) (lpdf0 x)]
                 [(< x xMax)  (lpdf1 x)]
                 [else         -inf.0]))]
            [else
             ;; lower band
             (define x00 (max .20 (flexpt (fl* 2.75 k) (fl/ -1. (- k 1.)))))
             (define x01 (if (< k 3.8) 1. (min 1. (max (+ 0.778 (* 0.006 k))
                                                       (flexpt (min (- k 3.) 8.) (fl/ -1. k))))))

             ;; upper band
             (define x10 (if (< k 3.8) 1. (flexpt (min (- k 3.) 3.) (fl/ 1. 2. k))))
             (define x11 (flexpt (- k .5) (fl/  1. k)))
             (λ (x)
               (cond
                 [(or (< x00 x x01) (< x10 x x11)) (lpdf6 x)]
                 [(< xx x XX)                      (lpdf0 x)]
                 [(< x xMax)                       (lpdf1 x)]
                 [else                              -inf.0]))])]))
     (define smallX : (-> Flonum Flonum) ;; x < 1
       (cond
         [(< k 1.)
          (define x_min (min (/ k 8) (- 0.21 (* 0.2 k))))
          (define x_max (min k (- 1.1 k)))
          (cond
            [(<= k +max-subnormal.0)
             (λ (x) (if (or (<= x +max-subnormal.0) (< x_min x x_max))
                        (lpdf6 x)
                        (lpdf0 x)))]
            [(<= k .5)
             (λ (x) (if (< x_min x x_max) (lpdf6 x) (lpdf0 x)))]
            [else
             (λ (x) (if (< x_min x x_max) (lpdf6 x) (lpdf1 x)))])]
         [else
          bigX]))


     (λ (x)
       (cond
         [(< x 0.) -inf.0]
         [(= x 1.) (- lnK 1.)]
         [(< x 1.) (smallX x)]
         [else     (bigX x)]))]
    [else (λ (x) +nan.0)]))

(: make-std-weibull-pdf (-> Real (PDF Real)))
(define (make-std-weibull-pdf k)
  (define F_nrm (make-std-nrm-weibull-pdf (fl k)))
  (define F_log (make-std-log-weibull-pdf (fl k)))
  (λ (x [log? #f]) (if (nan? x) +nan.0 ((if log? F_log F_nrm) (fl x)))))

(: make-weibull-pdf : (case-> (Real           -> (PDF Real))
                              (Real Real Real -> (PDF Real))))
(define make-weibull-pdf
  (case-lambda
    [(k) (make-std-weibull-pdf k)]
    [(k d s)
     (cond
       [(= s 0) (λ (x [log? #f]) (fldelta-pdf (fl d) (fl x) log?))]
       [else
        (define F (make-std-weibull-pdf k))
        (define S (fl (abs s)))
        (λ (x [log? #f])
          (define q (F (/ (- x d) s) log?))
          (if log? (fl- q (fllog S)) (fl/ q S)))])]))

;;;
;;; CDF (->* (Real) (Any Any) Flonum))
;;;

(: make-std-weibull-cdf (-> Real (CDF Real)))
(define (make-std-weibull-cdf k)
  (let ([k (fl k)])
    (cond
      [(< 0. k)
       (define x_min (flexpt 2. (fl/ 1. k)))
       ; (fllog +max.0) ~> 714.
       (define x_max (flexpt 714. (fl/ 1. k)))

       (define (F_nrm [x : Flonum] [1-p? : Any]) : Flonum
         (cond
           [(<= x 0.)
            (if 1-p? 1. 0.)]
           [1-p?
            (if (< x_min x x_max)
                (let*-values ([(ah al) (fl2log x 0.)]
                              [(bh bl) (fl2* ah al k)]
                              [(ch cl) (fl2exp bh bl)]
                              [(dh dl) (fl2exp (- ch) (- cl))])
                  dh)
                (flexp (- (flexpt x k))))]
           [else
            (- (flexpm1 (- (flexpt x k))))]))
       (define (F_log [x : Flonum] [1-p? : Any]) : Flonum
         (cond
           [(<= x 0.)
            (if 1-p? 0. -inf.0)]
           [1-p?
            (- (flexpt x k))]
           [else
            (define -x^k (- (flexpt x k)))
            (cond
              [(< x_min x x_max)
               (let*-values ([(ah al) (fl2log x 0.)]
                             [(bh bl) (fl2* ah al k)]
                             [(ch cl) (fl2exp bh bl)]
                             [(dh dl) (fl2exp (- ch) (- cl))])
           
                 (fllog1p (- dh)))]
              [(<= -max-subnormal.hi -x^k)
               (fl* k (fllog x))]
              [else
               (lg1- -x^k)])]))
       
       (λ (x [log? #f] [1-p? #f]) ((if log? F_log F_nrm) (fl x) 1-p?))]
      [else
       (λ (x [log? #f] [1-p? #f]) +nan.0)])))

(: make-weibull-cdf : (case-> (Real           -> (CDF Real))
                              (Real Real Real -> (CDF Real))))
(define make-weibull-cdf
  (case-lambda
    [(k) (make-std-weibull-cdf (fl k))]
    [(k d s)
     (cond
       [(= s 0) (λ (x [log? #f] [1-p? #f]) (fldelta-cdf (fl d) (fl x) log? 1-p?))]
       [else
        (define F (make-std-weibull-cdf (fl k)))
        (λ (x [log? #f] [1-p? #f])
          (F (/ (- x d) s) log? 1-p?))])]))

;;;
;;; inverse-CDF (->* (Real) (Any Any) Flonum))
;;;

(: make-std-nrm-weibull-inverse-cdf (-> Flonum (-> Flonum Any Flonum)))
(define (make-std-nrm-weibull-inverse-cdf k)
  (cond
    [(< 0. k)
     (begin-encourage-inline
       (define (icdf-0 [p : Flonum]) : Flonum
         (flexpt (- (fllog p)) (/ k)))
       (define (icdf-3 [p : Flonum]) : Flonum
         (define-values (hi lo) (fl2log p 0.))
         (define-values (Kh kl) (fl//error 1. k))
         (fl* (flexpt+ (- hi) (- lo) Kh)
              (flexpt+ (- hi) (- lo) kl)))
       
       (define (icdf+0 [p : Flonum]) : Flonum
         (flexpt (- (fllog1p (- p))) (/ k)))
       (define (icdf+1 [p : Flonum]) : Flonum
         (define-values (hi lo) (fl2log1p (- p) 0.))
         (flexpt+ (- hi) (- lo) (/ k)))
       (define (icdf+3 [p : Flonum]) : Flonum
         (define-values (hi lo) (fl2log1p (- p) 0.))
         (define-values (Kh kl) (fl//error 1. k))
         (fl* (flexpt+ (- hi) (- lo) Kh)
              (flexpt+ (- hi) (- lo) kl))))

     (define inrange-p :(-> Flonum Flonum)
       (cond
         [(< 1e8 k) icdf+0]
         [(< 1e3 k) icdf+1]
         [(<= 1e-36 k)
          (define p_min (flprev (- (flexpm1 (- (flexpt +min.0 k))))))
          (define p_max (flnext (- (flexpm1 (- (flexpt +max.0 k))))))
          (λ (p) (if (<= p_min p p_max) (icdf+3 p) (icdf+0 p)))]
         [else  icdf+0]))
     (define inrange-1-p : (-> Flonum Flonum)
       (cond
         [(< 1e1 k) icdf-0]
         [(<= 1e-36 k)
          (define p_min (flprev (flexp (- (flexpt +max.0 k)))))
          (define p_max (flnext (flexp (- (flexpt +min.0 k)))))
          (λ (p) (if (<= p_min p p_max) (icdf-3 p) (icdf-0 p)))]
         [else icdf-0]))

     (define inrange : (-> Flonum Any Flonum)
       (cond
         [(= k +inf.0) (λ (p 1-p?) 1.0)]
         [else         (λ (p 1-p?) ((if 1-p? inrange-1-p inrange-p) p))]))
     
     (λ (p 1-p?)
       (if (<= 0. p 1.)
           (cond
             [(= p 0.) (if 1-p? +inf.0 0.)]
             [(= p 1.) (if 1-p? 0. +inf.0)]
             [else (inrange p 1-p?)])
           +nan.0))]
    [else
     (λ (p 1-p?) +nan.0)]))

(: make-std-log-weibull-inverse-cdf (-> Flonum (-> Flonum Any Flonum)))
(define (make-std-log-weibull-inverse-cdf k)
  (cond
    [(< 0. k)
     (begin-encourage-inline
       (define (ilcdf+0 [p : Flonum]) : Flonum
         (flexpt (- (lg1- p)) (/ k)))
       (define (ilcdf+1 [p : Flonum]) : Flonum
         (flexp (/ p k)))
       (define (ilcdf+2 [p : Flonum]) : Flonum
         (define-values (H L) (fl//error p k))
         (define-values (E e) (fl2exp H L))
         E)
       (define (ilcdf+3 [p : Flonum]) : Flonum
         (define-values (E e) (fl2exp p 0.))
         (define-values (L l) (fl2log1p (- E) (- e)))
         (define-values (I i) (fl2log (- L) (- l)))
         (define-values (D d) (fl2/ I i k))
         (define-values (R r) (fl2exp D d))
         R)
       (define (ilcdf+4 [p : Flonum]) : Flonum
         (define-values (I i) (fl2log (- p) 0.))
         (define-values (L l) (fl2log (- I) (- i)))
         (define-values (D d) (fl2/ L l k))
         (define-values (R r) (fl2exp D d))
         R))
     
     (define inrange-p : (-> Flonum Flonum)
       (cond
         [(< k 1.)
          (define LIMIT0 (lg1- (- (flexpt +min.0 k))))
          (define LIMIT1 (* (/ 8.3e1 (flexpt k .035)) (lg1- (- (flexpt +max.0 (* k 1e-19))))))
          (define LIMIT2 (lg1- (- (flexpt +max.0 k))))
          (define LIMIT3 (* (fllog +min.0) (* k 5e-18)))
          (λ (p)
            (cond
              [(<  p LIMIT0) (ilcdf+1 p)]
              [(<  p LIMIT1) (ilcdf+2 p)]
              [(<= p LIMIT2)
               (cond
                 [(<= p LIMIT3) (ilcdf+3 p)]
                 [else          (ilcdf+4 p)])]
              [else          (ilcdf+0 p)]))]
         [else
          (define LIMIT1 (* (/ 8.3e1 (flexpt k .035)) (lg1- (- (flexpt +max.0 (* k 1e-19))))))
          (define LIMIT3 (* (fllog +min.0) (* k 5e-18)))
          (define LIMIT4 (* (fllog +min.0) (* k 1e-3)))
          (define LIMIT5 (* (fllog +min.0) k))
          (define LIMIT6 (- (* (flexpt 2. 20.) (fllog +max.0))))
          (define LIMIT7 (* 8e1 (lg1- (- (flexpt +max.0 (* k 1e-19))))))
          (define (rst [p : Flonum]) : Flonum
            (cond
              [(< LIMIT1 p LIMIT4)
               (cond
                 [(<= p LIMIT3) (ilcdf+3 p)]
                 [else          (ilcdf+4 p)])]
              [(< p LIMIT5) (ilcdf+1 p)]
              [(< LIMIT6 p LIMIT7)
               (cond
                 [(< LIMIT5 p LIMIT4) (ilcdf+2 p)]
                 [else                (ilcdf+1 p)])]
              [else         (ilcdf+0 p)]))
          (cond
            [(< k 10.)
             (λ (p)
               (cond
                 [(< -1. p)
                  (cond
                    [(<= p LIMIT3) (ilcdf+3 p)]
                    [else          (ilcdf+4 p)])]
                 [else (rst p)]))]
            [else rst])]))
     (define inrange-1-p : (-> Flonum Flonum)
       (let ([p_min_low (flprev (- (flexpt +max.0 k)))]
             [p_max_low (flnext (- (flexpt +max.0 (/ k 710.))))]
             [p_min_hig (flprev (- (flexpt +min.0 (/ k 710.))))]
             [p_max_hig (flnext (- (flexpt +min.0 k)))])
         (λ (p)
           (cond
             [(if (< p -1)
                  (< p_min_low p p_max_low)
                  (< p_min_hig p p_max_hig))
              (define-values (L l) (fl2log (- p) 0.))
              (define-values (D d) (fl2/ L l k))
              (define-values (E e) (fl2exp D d))
              E]
             [else
              (flexpt (- p) (/ k))]))))
     
     (define inrange : (-> Flonum Any Flonum)
       (cond
         [(= k +inf.0) (λ (p 1-p?) 1.0)]
         [else (λ (p 1-p?) ((if 1-p? inrange-1-p inrange-p) p))]))
     
     (λ (p 1-p?)
       (if (<= -inf.0 p 0.)
           (cond
             [(= p -inf.0) (if 1-p? +inf.0 0.)]
             [(= p    0.0) (if 1-p? 0. +inf.0)]
             [else (inrange p 1-p?)])
           +nan.0))]
    [else
     (λ (p 1-p?) +nan.0)]))

(: make-std-weibull-inverse-cdf (-> Real (Inverse-CDF Flonum)))
(define (make-std-weibull-inverse-cdf k)
  (define F_nrm (make-std-nrm-weibull-inverse-cdf (fl k)))
  (define F_log (make-std-log-weibull-inverse-cdf (fl k)))
  (λ (p [log? #f] [1-p? #f]) ((if log? F_log F_nrm) (fl p) 1-p?)))

(: make-weibull-inverse-cdf : (case-> (Real           -> (Inverse-CDF Flonum))
                                      (Real Real Real -> (Inverse-CDF Flonum))))
(define make-weibull-inverse-cdf
  (case-lambda
    [(k) (make-std-weibull-inverse-cdf (fl k))]
    [(k d s)
     (cond
       [(= s 0) (λ (x [log? #f] [1-p? #f]) (fldelta-inv-cdf (fl d) (fl x) log? 1-p?))]
       [else
        (define F (make-std-weibull-inverse-cdf (fl k)))
        (λ (p [log? #f] [1-p? #f])
          (define x (F p log? 1-p?))
          (fl (+ d (* x s))))])]))

;;;
;;; sampler (-> Integer FlVector)
;;;

(: make-weibull-sampler : (case-> (Real           -> (Integer -> FlVector))
                                  (Real Real Real -> (Integer -> FlVector))))
(define make-weibull-sampler
  (case-lambda
    [(k)
     (define F (make-std-nrm-weibull-inverse-cdf (fl k)))
     (λ (n)
       (define PRNG (current-pseudo-random-generator))
       ;; from my test : (random) < 0.5 seems slighly faster than (random 2) = 1 (when using fl-operators)
       ;; this needs to improve: (random) only generates numbers between 2.3283...e-10 and 0.99...997671695
       (build-flvector n (λ (_) (F (fl* 0.5 (random PRNG)) (fl< (random PRNG) .5)))))]
    [(k d s)
     (define D (fl d))
     (cond
       [(= s 0) (λ (n) (build-flvector n (λ (_) D)))]
       [else
        (define F (make-std-nrm-weibull-inverse-cdf (fl k)))
        (λ (n)
          (define PRNG (current-pseudo-random-generator))
          (define S (fl s))
          (build-flvector n (λ (_) (+ D (* S (F (fl* 0.5 (random PRNG)) (fl< (random PRNG) .5)))))))])]))
