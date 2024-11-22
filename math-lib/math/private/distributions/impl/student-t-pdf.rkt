#lang typed/racket
(require "student-t-utils.rkt"
         "../../../flonum.rkt"
         "../dist-struct.rkt"
         "normal-pdf.rkt")

(provide make-student-t-pdf)
                     
;;; Student t distribution

;; Parameters
;   ν   - degrees of freedom
;   μ   - location parameter
;   σ   - scale parameter

;; Domains
;   ν   - any positive real number
;   μ   - any real number
;   σ   - any real number


;; Probability Density Function (PDF)


; The density is proportional to
;             ν                           
; f(x) = ( -------- )^( (1+ν)/2 )         
;           (ν+x²)                        

; The proportionality constant is

;     Γ((ν+1)/2)                    
; -----------------                 
;  sqrt(πν) Γ(ν/2)                  

; Using the Β function this can also be written as

;          1                         
; --------------------               
;  sqrt(ν) B(1/2, ν/2)               

;; Reduction to one parameter

; If X ~ t(μ,σ,ν) then (X-μ)/σ ~ t(ν).
; This means we can concentrate on the one parameter version t(ν).
;     t(ν) = t(0, 1, ν)

;  t(ν)      - "the" Student t distribution
;  t(ν,μ,σ)  - the generalized Student t distribution
;              also called the location-scale-t-distribution



;; reimplentation of flexpt1p, but requiring the +1 to be done beforehand in fl2 (B & b)
;; and without checking of boundaries (-0.5 < B_b < +inf.0)
(: fl2expt (Flonum Flonum Flonum -> Flonum))
(define (fl2expt B b e) : Flonum
  (fl/ (flexpt B e)
       (flexp (fl* e (fllog1p (fl- (fl/ b B)))))))


;;;
;;; Implementation
;;;

(: make-student-t-pdf : (case-> (Real           -> (PDF Real))
                                (Real Real Real -> (PDF Real))))
(define make-student-t-pdf
  (case-lambda
    ; Y ~ σX+μ
    [(ν μ σ)
     (define f (make-student-t-pdf ν))
     (let ([μ (fl μ)] [σ (fl σ)])
       (λ (y [log? #f])
         (let ([y (fl y)])
           (cond
             [log? (define x     (fl/ (fl- y μ) σ))
                   (define log-p (fl- (f x #t) σ))
                   log-p]
             [else (define x     (fl/ (fl- y μ) σ))
                   (define p     (fl/ (f x) σ))
                   p]))))]
    
    ; X ~ t(ν)
    [(ν)
     (let ([BF (beta1/2 ν)]      ;; if (exact? ν) is provided, beta1/2 will be calculated exact,
           [LBF (logbeta1/2 ν)]  ;;  this can take a long time for big ν
           [ν (fl ν)])

       ;; *************************************************************************************************
       ;; ** PDF FUNCTIONS
       ;; *************************************************************************************************
       (define √ν (flsqrt ν))
       (define expo (* -.5 (fl+ 1. ν)))
       (define proportionality-constant (fl/ 1. (* √ν BF)))
       (define x-bnd (max 38.8 (flsqrt (fl* ν (fl- (flexp (fl/ (fl* -2. (fllog +min.0)) (fl+ 1. ν))) 1.)))))

       (: pdf1 (Flonum -> Flonum))
       (define (pdf1 x)
         (define base (flexpt+ ν (flexpt x 2.) ν))
         (fl/ (flexpt √ν ν)
              (flsqrt (fl+ (fl* base ν) (fl* x base x)))
              BF))

       (: pdf2 (Flonum -> Flonum))
       (define (pdf2 x)
         (fl* proportionality-constant (fl/ √ν x)))

       (: pdf3 (Flonum -> Flonum))
       (define (pdf3 x)
         (define base (fl/ √ν x))
         (fl* proportionality-constant base (flexpt base ν)))

       (: pdf4~ (Flonum -> Flonum))
       (define (pdf4~ x)
         (define x/ν (fl/ x ν))
         (if (< .1 x/ν 10.)
             (fl/ (flexpt x 2.) ν)
             (* x/ν x)))
       (: pdf4-e (Flonum -> Flonum))
       (define (pdf4-e x)
         (define base (pdf4~ x))
         (fl* proportionality-constant (flexpt1p base expo)))
       (: pdf4-o (Flonum -> Flonum))
       (define (pdf4-o x)
         (define base (pdf4~ x))
         (fl* proportionality-constant (flexpt1p base -.5) (flexpt1p base (* -.5 ν))))

       ;; limited f2 calc, same speed as 4, since fl2expt is similar to flext1p
       (: pdf5~ (Flonum -> (Values Flonum Flonum)))
       (define (pdf5~ x)
         (define-values (X² x²) (flsqr/error x))
         (define-values (X²/N x²/ν) (fl2/ X² x² ν))
         (fl2+ X²/N x²/ν 1.))
       (: pdf5-e (Flonum -> Flonum))
       (define (pdf5-e x)
         (define-values (A a) (pdf5~ x))
         (fl* proportionality-constant (fl2expt A a expo)))
       (: pdf5-o (Flonum -> Flonum))
       (define (pdf5-o x)
         (define-values (A a) (pdf5~ x))
         (fl* proportionality-constant (fl2expt A a -.5) (fl2expt A a (* -.5 ν))))

       ;; full f2 calc, 2.5 times slower than 4
       (: pdf6 (Flonum -> Flonum))
       (define (pdf6 x)
         (define-values (X² x²) (flsqr/error x))
         (define-values (X²/N x²/ν) (fl2/ X² x² ν))
         (define-values (LG lg) (fl2log1p X²/N x²/ν))
         (define-values (E e) (fl2* LG lg (* -.5 ν) -.5))
         (fl* proportionality-constant (flexp E) (flexp e)))
               
       (: pdf : (Flonum -> Flonum))
       (define pdf
         (cond
           [(or (<= ν 0.) (flnan? ν)) (λ (x) +nan.0)]
           [(< ν 1.12e-308)
            (λ (x)
              (cond
                [(< x 1e15)    (flexp (log-pdf x))]
                [(< x 1.1e77)  (pdf1 x)]
                [else          (pdf2 x)]))]
           [(< ν 1e-20)
            (λ (x)
              (cond
                [(= x 0.)       proportionality-constant]
                [(< x 1.1e77)  (pdf1 x)]
                [else          (pdf2 x)]))]
           [(< ν 1.)
            (λ (x)
              (cond
                [(= x 0.)       proportionality-constant]
                [(< x 1.1e77)  (pdf1 x)]
                [else          (pdf3 x)]))]
           [(< ν 3.5)
            (if (fleven? ν)
                (λ (x)
                  (cond
                    [(= x 0.)       proportionality-constant]
                    [(< x (min 1e30 x-bnd)) (pdf4-e x)]
                    [else         (pdf3 x)]))
                (λ (x)
                  (cond
                    [(= x 0.)       proportionality-constant]
                    [(< x (min 1e30 x-bnd)) (pdf4-o x)]
                    [else         (pdf3 x)])))]
           [(< ν 5e16)
            (if (and (< ν 5e15) (fleven? ν))
                (λ (x)
                  (cond
                    [(= x 0.)       proportionality-constant]
                    [(< x 2.5)     (pdf4-e x)]
                    [(< x x-bnd)   (pdf5-e x)]
                    [(flnan? x)    +nan.0]
                    [else          0.0]))
                (λ (x)
                  (cond
                    [(= x 0.)       proportionality-constant]
                    [(< x 2.5)     (pdf4-o x)]
                    [(< x x-bnd)   (pdf5-o x)]
                    [(flnan? x)    +nan.0]
                    [else          0.0])))]
           [(< ν 1e20)
            (λ (x)
              (cond
                [(= x 0.)       proportionality-constant]
                [(< x 2.5)     (pdf4-o x)]
                [(< x x-bnd)   (pdf6 x)]
                [(flnan? x)    +nan.0]
                [else          0.0]))]
           [else standard-flnormal-pdf]))

       ;; *************************************************************************************************
       ;; ** LOG PDF FUNCTIONS
       ;; *************************************************************************************************
       (define log-proportionality-constant
         (cond
           ;(fl/ +max-subnormal.0 2.)
           [(< 1.1125369292536007e-308 ν) (fllog proportionality-constant)]
           ;; below is not good enough in general, bad if ν -> +inf (300 ulp)
           ;; but both above and below are very inacurate in the range 3->15 (5000 ulp)
           ;; however this doesn't seem to matter for the final result
           [else (fl- (fl* -.5 (fllog ν)) LBF)]))
       (define lν (fllog ν))
       (define l√v (fl* 0.5 lν))
       (define νl√v (fl* 0.5 ν lν))
       (define lx-bnd (exp (+ 354.891356446692 l√v))) ;; x²/ν < +inf.0

       (: log1 (Flonum -> Flonum))
       (define (log1 x)
         (define base (flexpt+ ν (flexpt x 2.) ν))
         (flsum (list νl√v
                      (fl* -.5 (fllog (fl+ (fl* base ν) (fl* x base x))))
                      (fl- LBF))))
               
       (: log2 (Flonum -> Flonum))
       (define (log2 x)
         (define lx (fl- (fllog x)))
         (flsum (list log-proportionality-constant
                      l√v lx νl√v (fl* ν lx))))
               
       (: log3 (Flonum -> Flonum))
       (define (log3 x)
         (fl+ log-proportionality-constant
              (fl* -1. expo (fllog (fl/ ν (fl+ ν (flexpt x 2.)))))))
               
       (: log6 (Flonum -> Flonum))
       (define (log6 x)
         (fl+ log-proportionality-constant
              (fl* expo (fllog1p (fl* (fl/ x ν) x)))))
               
       (: log-pdf : (Flonum -> Flonum))
       (define log-pdf
         (cond
           [(or (<= ν 0.) (flnan? ν)) (λ (x) +nan.0)]
           [(< ν 1.)
            (λ (x)
              (cond
                [(< x 1.)   (log3 x)]
                [(< x 1e77) (log1 x)]
                [else       (log2 x)]))]
           [else
            (λ (x)
              (cond
                [(< x lx-bnd) (log6 x)]
                [else         (log2 x)]))]))
               
       (: result-pdf : (PDF Real))
       (define (result-pdf x [log? #f])
         (let ([x (fl x)])
           (if log? (log-pdf (flabs x)) (pdf (flabs x)))))               
       result-pdf)]
    ))
