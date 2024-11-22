#lang typed/racket
(require "student-t-utils.rkt"
         "normal-cdf.rkt"
         "../dist-struct.rkt"
         "../../functions/incomplete-beta.rkt"
         "../../../flonum.rkt")

(provide make-student-t-cdf)
                     
;;; Student t distribution

;; Parameters
;   ν   - degrees of freedom
;   μ   - location parameter
;   σ   - scale parameter

;; Domains
;   ν   - any positive real number
;   μ   - any real number
;   σ   - any real number


;; The Cumulative distribution function (CDF)

; For t>0
;   F(t) = 1 - 0.5 I_x(t) (ν/2, 0.5),
;   where x(t) = ν/(t²+ν).

; Here I is the regularized incomplete beta function.



;; *************************************************************************************************
;; ** CDF IMPLEMENTAION
;; *************************************************************************************************
(: make-student-t-cdf : (case-> (Real           -> (CDF Real))
                                (Real Real Real -> (CDF Real))))
(define make-student-t-cdf
  (case-lambda
    ; Y ~ σX+μ
    [(ν μ σ)
     (define F (make-student-t-cdf ν))
     (let ([μ (fl μ)] [σ (fl σ)])
       (λ (y [log? #f] [1-p? #f])               
         (define x (/ (- (fl y) μ) σ))
         (F x log? 1-p?)))]
    
    ; X ~ t(ν)
    [(ν)
     (let ([ν (fl ν)])
       (define ν/2 (fl/ ν 2.))

       ;; *************************************************************************************************
       ;; ** CDF FUNCTIONS
       ;; *************************************************************************************************
       (: cdf4 (Flonum -> Flonum))
       (define (cdf4 x)
         (flbeta1/2-regularized ν x (/ (+ 1. (/ (/ ν x) x))) (/ (+ 1. (* x (/ x ν)))) #f #t #t))

       (: cdf4-hypergeom (Flonum -> Flonum))
       (define (cdf4-hypergeom x)
         (define X (/ (+ 1. (/ (/ ν x) x))))
         (define Y (/ (+ 1. (* x (/ x ν)))))
         (define log-X (if (< X 0.5) (fllog X) (fllog1p (- Y))))
         (define log-Y (if (< Y 0.5) (fllog Y) (fllog1p (- X))))
         (define log? #f) (define 1-? #t) (define r? #t)
  
         (if (or (< 1e4 ν) (< -0.5 x)) +nan.0
             (let ([A ν][x Y][y X][log-x log-Y][log-y log-X])
               (define a+2/2 (fl+ A 2.0))
               (define: z : Flonum
                 (let loop ([z 0.0] [dz 1.0] [2n 0.0] [i -1.0])
                   ;(printf "z = ~v  dz = ~v  i = ~v~n" z dz i)
                   (define new-dz (fl* (fl- dz (fl/ dz (fl+ a+2/2 2n))) x))
                   (cond [(zero? i)  (fl+ z new-dz)]
                         [else
                          (let ([i  (if (and (i . fl< . 0.0)
                                             ((flabs new-dz) . fl<= . (flabs dz))
                                             ((flabs new-dz) . fl<= . (fl* (fl* 0.5 epsilon.0) (flabs z))))
                                        3.0
                                        (fl- i 1.0))])
                            (loop (fl+ z new-dz) new-dz (fl+ 2n 2.0) i))])))
               (define c (/ (* (flexpt x (* 0.5 A)) (flsqrt y)) (beta1/2 A)))
               ;(println (list a+b/2 a+2/2 (/ a+b/2 a+2/2) z c))
               (fl/ (fl+ c (fl* z c)) A))))

       (: cdf-largex (Flonum -> Flonum))
       (define (cdf-largex x)
         ;;(B (ν/ν+x²) ; ν/2 1/2)
         ;; 0th therm of continuous fraction when ν+x²≈x²
         (define _x_ (flabs x))
         (fl/ (flexpt (fl/ (flsqrt ν) _x_) ν)
              (fl* ν (beta1/2 ν))))

       (: cdf-scalednorm (Flonum -> Flonum))
       (define (cdf-scalednorm x)
         (define V (/ 1. 4. ν))
         (standard-flnormal-cdf (/ (* x (- 1. V))
                                   (flsqrt (+ 1. (* x x 2. V))))))

       (: invert ((Flonum -> Flonum) -> (Flonum -> Flonum)))
       (define (invert f)
         (λ (x)
           (cond
             [(flnan? x) +nan.0]
             [else
              ((if (< x 0.) values (λ ([x : Flonum]) (fl- 1. x)))
               (let ([x (if (< x 0.) x (- x))])
                 (cond
                   [(= x -inf.0) 0.]
                   [(< -1e-17 x) 0.5]
                   [else         (f x)])))])))
       
       (: cdf : (Flonum -> Flonum))
       (define cdf
         (cond
           [(flnan? ν) (λ (x) +nan.0)]
           [(< ν 1e-19) (λ (x) 0.5)]
           [(< ν 1e3)
            (define large-lim (if (< ν 1.) (max -1e21 (* -1e8 ν)) -1e21))
            (define geo-lim (if (< ν 2e2) -1e0 -1e1))
            (invert
             (λ (x)
              (cond
                [(< x large-lim) (cdf-largex x)]
                [(< x geo-lim)   (cdf4-hypergeom x)]
                [else            (cdf4 x)])))]
           [(< ν 1e7)
            (invert
             (λ (x)
               (cond
                 [(< x -1e2) 0.]
                 [else       (cdf4 x)])))]
           [(< ν 1e20)
            (define norm-lim (* -1e-3 (flexpt ν #i1/3)))
            (invert
             (λ (x)
               (cond
                 [(< norm-lim x) (cdf-scalednorm x)]
                 [(< x -1e2)     0.]
                 [else           (cdf4 x)])))]
           [else
            (λ (x) (standard-flnormal-cdf x))]
           ))

       ;; *************************************************************************************************
       ;; ** CDF LOG FUNCTIONS
       ;; *************************************************************************************************
       (: lcdf1 (Flonum -> Flonum))
       (define (lcdf1 x)
         (fllog (cdf x)))

       (: lcdf1.1 (Flonum -> Flonum))
       (define (lcdf1.1 x)
         (fllog1p (- (cdf (- x)))))

       (: lcdf2 (Flonum -> Flonum))
       (define (lcdf2 x)
         (define z (fl/ ν (fl+ (fl* x x) ν)))
         (define ν/2 (/ ν 2.))
         (if (< x 0.)
             (+ (fllog 0.5) (log-beta-inc ν/2 0.5 z #f #t))
             (fllog1p (* -.5 (beta-inc ν/2 0.5 z #f #t)))))

       (: lcdf5 (Flonum -> Flonum))
       (define (lcdf5 x)
         (define z (fl/ (fl+ 1. (fl/ (fl/ ν x) x))))
         (define ν/2 (/ ν 2.))
         (if (< x 0.)
             (+ (fllog 0.5) (log-beta-inc 0.5 ν/2 z #t #t))
             (fllog1p (* -.5 (beta-inc 0.5 ν/2 z #t #t)))))

       (: lcdf9 (Flonum -> Flonum))
       (define (lcdf9 x)
         ;;(B (ν/ν+x²) ; ν/2 1/2)
         ;; 0th therm of continuous fraction when ν+x²≈x²
         (define _x_ (flabs x))
         (fl- (fl* ν (fl- (fl* 0.5 (fllog ν)) (fllog _x_)))
              (fl+ (fllog ν) (logbeta1/2 ν))))

       (define lg1/2 (fllog 0.5))
       (define 0.5sqrtmax (* .5 (flsqrt +max.0)))

       (: logcheck ((Flonum -> Flonum) -> (Flonum -> Flonum)))
       (define (logcheck f)
         (λ (x)
           (cond
             [(flnan? x) +nan.0]
             [(< (flabs x) 1e-17) lg1/2]
             [else (f x)])))
       
       (: log-cdf : (Flonum -> Flonum))
       (define log-cdf
         (cond
           [(nan? ν) (λ (x) +nan.0)]
           [(< ν 1e-20) (λ (x) lg1/2)]
           [(< ν 1.)
            (logcheck
             (λ (x)
               (cond
                 [(< x 0.) (lcdf1 x)]
                 [else     (lcdf1.1 x)])))]
           [(< ν 1e18)
            (define large-lim (* -1e7 (flexpt ν 0.5)))
            (define log-lim   (* -1e0 (flexpt ν 0.5)))
            (logcheck
             (λ (x)
               (cond
                 [(< x large-lim) (lcdf9 x)]
                 [(< x log-lim)   (lcdf2 x)]
                 [(< x -10.)      (lcdf5 x)]
                 [(< x 0.)        (lcdf1 x)]
                 [else            (lcdf1.1 x)])))]
           [(< ν 0.5sqrtmax)
            (define large-lim (if (< ν 1e20) (* -1e7 (flexpt ν 0.5)) (- 0.5sqrtmax)))
            (define log-lim (* -1e0 (flexpt ν 0.5)))
            (define mid-lim (* -1e-8 (flexpt ν 0.5)))
            (logcheck
             (λ (x)
               (cond
                 [(< x large-lim) (lcdf9 x)]
                 [(< x log-lim)   (lcdf2 x)]
                 [(< x mid-lim)   (lcdf5 x)]
                 [else
                  (standard-flnormal-log-cdf x)])))]
           [(< ν +inf.0)
            (define large-lim (- 0.5sqrtmax))
            (define log-lim (* -1e0 (flexpt ν 0.5)))
            (logcheck
             (λ (x)
               (cond
                 [(< x large-lim) (lcdf9 x)]
                 [(< x log-lim)   (lcdf2 x)]
                 [else
                  (standard-flnormal-log-cdf x)])))]
           [else
            (λ (x) (standard-flnormal-log-cdf x))]))
               
       (: result-cdf : (CDF Real))
       (define (result-cdf x [log? #f] [1-p? #f])
         (let* ([x    (fl (if 1-p? (- x) x))])
           (cond
             [log? (log-cdf x)]
             [else (cdf x)])))
       result-cdf)]
    ))

