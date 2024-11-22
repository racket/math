#lang typed/racket
(require "student-t-pdf.rkt"
         "student-t-cdf.rkt"
         "../gamma-dist.rkt"         
         "../../../flonum.rkt"
         "../dist-struct.rkt"
         "normal-inv-cdf.rkt"
         "normal-random.rkt")

(provide make-student-t-pdf
         make-student-t-cdf
         make-student-t-inverse-cdf
         flstudent-t-sample)
                     
;;; Student t distribution

;; Parameters
;   μ   - location parameter
;   σ   - scale parameter
;   ν   - degrees of freedom

;; Domains
;   μ   - any real number
;   σ   - any positive real number
;   ν   - any positive real number


;;;
;;; Implementation
;;;

;; (define-type (Inverse-CDF Out)
;;   (case-> (Real -> Out)
;;           (Real Any -> Out)
;;           (Real Any Any -> Out)))

(: make-student-t-inverse-cdf : (case-> (Real           -> (Inverse-CDF Flonum))
                              (Real Real Real -> (Inverse-CDF Flonum))))
(define make-student-t-inverse-cdf
  (case-lambda
    ; Y ~ σX+μ
    [(ν μ σ)
     (define inv-F (make-student-t-inverse-cdf ν))
     (let ([μ (fl μ)] [σ (fl σ)])
       (λ (p [log? #f] [1-p? #f])
         (define x (inv-F p log? 1-p?))
         (define y (+ (* σ x) μ))
         (fl y)))]
    
    ; X ~ t(ν)
    [(ν)
     (define lg1/2 (- (fllog 2.)))
     
     (: 1-?-swap (-> (Real Any -> Flonum) (->* (Real) (Any Any) Flonum)))
     (define (1-?-swap f)
       (λ (p* [log? #f] [1-? #f])
         (let ([p (fl p*)])
           (if log?
               (cond
                 [(= p -inf.0) -inf.0]
                 [(= p   0.)   +inf.0]
                 [(= p lg1/2)       0.0]
                 [(< 0. p)     +nan.0]
                 [else        (if 1-? (- (f p* log?))      (f p* log?))])
               (cond
                 [(= p 0.)     -inf.0]
                 [(= p 1.)     +inf.0]
                 [(< p 0.)     +nan.0]
                 [(< 1. p)     +nan.0]
                 [(= p 0.5)       0.0]
                 [(<= p 0.5)  (if 1-? (- (f p* log?))      (f p* log?))]
                 [else        (if 1-? (f (- 1 p*) log?) (- (f (- 1 p*) log?)))])))))
     (case ν
       ; special cases
       [(1 2 4)
        (: plain-inv-F : (Flonum -> Flonum))
        (define plain-inv-F
          (case ν
            [(1) (λ (p)
                   (fltan (* pi (fl- p 0.5))))]
            [(2) (λ (p)
                   (define α (fl* 4. p (fl- 1. p)))
                   (* 2. (fl- p 0.5) (flsqrt (fl/ 2. α))))]
            [(4) (λ (p)
                   (define α (fl* 2. (flsqrt p) (flsqrt (fl- 1. p))))
                   (define q (fl/ (flcos (fl/ (flacos α) 3.)) α))
                   (fl* (flsgn (fl- p 0.5)) 2. (flsqrt (fl- q 1.))))]
            [else (λ (p) 0.0)])) ; happy type checking

        (1-?-swap
         (λ (p log?)
           (let ([p (fl p)])
             (plain-inv-F (if log? (flexp p) p)))))]
       [(+nan.0)
        (λ (p [log? #f] [1-? #f]) +nan.0)]
       [(+inf.0)
        (1-?-swap
         (λ (p log?)
           (if log?
               (standard-flnormal-inv-log-cdf (fl p))
               (standard-flnormal-inv-cdf (fl p)))))]
       ; general
       [else
        (define F (make-student-t-cdf ν))
        (let ([ν (fl ν)])
          ;; lower limit for result bigger than 0.
          (define-values (p- x-)
            (cond
              [(< ν 1e-20) (values 0.5 -max.0)]
              [(< ν 2e0)   (values (F -max.0 #f) -max.0)]
              [(let ([V (map fllog '(2e0   1e1  1e2 1e3 1e4   1e6))]
                     [X (map fllog '(1e162 1e33 1e4 6e1 5.1e1 39.))]
                     [T (fllog ν)])
               (for/or : (Option (List Flonum Flonum))
                 ([v0 (in-list V)]
                  [v1 (in-list (cdr V))]
                  [y0 (in-list X)]
                  [y1 (in-list (cdr X))]
                  #:when (<= v0 T v1))
                 (list 0. (- (flexp (+ y0 (* (/ (- T v0) (- v1 v0)) (- y1 y0))))))))
               => (λ (l) (apply values l))]
              [else   (values 0. -39.)]))
          ;; limits for the log version
          (define lgp- (F -max.0 #t))
          (define lgp+ (F +max.0 #t))

          (: inv- (Flonum -> Flonum))
          (define inv-
            ;; keep a table to help finding the brackets for the root finding algorithm
            ;; table is only filled with p<0.5 and result<0. (other branch symetrical)
            (let ([TBL : (Listof (Pair Flonum Flonum)) '((0.5 . -1e-17))])
              (λ (p) ; input will be p in range [0 -> 0.5]
                (cond
                  [(< p p-)     -inf.0]
                  [else
                   ;; try finding the brackets for rootfinder in the table
                   (: a0 (Option Flonum))(: b0 Flonum)(: f0 Flonum)
                   (define-values (a0 b0 f0)
                     (apply
                      values
                      (or (for/or : (Option (List Flonum Flonum Flonum))
                            ([q0 (in-list TBL)]
                             [q1 (in-list (cdr TBL))]
                             #:when (<= (car q0) p (car q1)))
                            (list (cdr q0) (cdr q1) (car q1)))
                          (list #f (cdar TBL) (caar TBL)))))

                   ;; extend the table if necessary and find a correct bracket
                   (define-values (a b)
                     (if a0
                         (values a0 b0)
                         (let ([fx : Flonum f0][E : Flonum 2.])
                           (let find-bracket : (Values Flonum Flonum) ([a : Flonum (* b0 2.)] [b : Flonum b0] [fb : Flonum f0])
                             ; Since the function h is monotone, this strategy works.
                             (define fa (F a))
                             (when (< 0.075 (- fx fa))
                               (set! fx fa)
                               (set! TBL (cons (cons fa a) TBL)))
                             (when (< E (fllog (/ fx fa)))
                               (set! fx fa)(set! E (* 2. E))
                               (set! TBL (cons (cons fa a) TBL)))
                             (cond
                               [(<= fa p fb) (values a b)]
                               [else
                                (find-bracket (max x- (* 2. a)) a fa)])))))

                   (flbracketed-root (λ (x) (fl- (F x) p)) a b)]))))

          (: lginv (Flonum -> Flonum))
          (define lginv
            ;; keep two distinct tables for the positive (> lg1/2) and negative (< lg1/2) branch
            (let ([TBL+ : (Listof (Pair Flonum Flonum)) (list (cons lg1/2  1e-16))]
                  [TBL- : (Listof (Pair Flonum Flonum)) (list (cons lg1/2 -1e-16))])
              (λ (p) ; input will be p in range [-inf.0 -> 0.0]
                (cond
                  [(= p lg1/2)     0.0]
                  [(< p lgp-)   -inf.0]
                  [(< lgp+ p)   +inf.0]
                  [else
                   (define <? (< p lg1/2))
                   (define Q (if <? TBL- TBL+))
                   (define updateQ
                     (if <? (λ ([pair : (Pair Flonum Flonum)]) (set! TBL- (cons pair TBL-)))
                         (λ ([pair : (Pair Flonum Flonum)]) (set! TBL+ (cons pair TBL+)))))
                   (define in? (if <? <= >=))

                   ;; try finding the brackets for rootfinder in the table
                   (: a0 (Option Flonum))(: b0 Flonum)(: f0 Flonum)
                   (define-values (a0 b0 f0)
                     (apply
                      values
                      (or (for/or : (Option (List Flonum Flonum Flonum))
                            ([q0 (in-list Q)]
                             [q1 (in-list (cdr Q))]
                             #:when (in? (car q0) p (car q1)))
                            (list (cdr q0) (cdr q1) (car q1)))
                          (list #f (cdar Q) (caar Q)))))
                   
                   ;; extend the table if necessary and find a correct bracket
                   (define-values (a b)
                     (if a0
                         (values a0 b0)
                         (let ([fx : Flonum f0][E : Flonum 2.])
                           (let find-bracket : (Values Flonum Flonum)
                             ; if p<lg2 => x < 0, Q- is sorted from -inf.0 to +inf.0
                             ([a : Flonum (* b0 E)] [b : Flonum b0] [fb : Flonum f0])
                             ; Since the function h is monotone, this strategy works.
                             (define fa (F a #t))
                             (when (< E (if <? (/ fa fx) (/ fx fa)))
                               (set! fx fa)(set! E (* 2. E))
                               (updateQ (cons fa a)))
                             
                             (cond
                               [(in? fa p fb)
                                (when (if <? (< (* fx 2.) lgp-)
                                          (< lgp+ (/ fx 2.)))
                                  (updateQ (if <? (cons lgp- -max.0) (cons lgp+ +max.0))))
                                (values a b)]
                               [else
                                (find-bracket (* 2. a) a fa)])))))
                   
                   (flbracketed-root (λ (x) (fl- (F x #t) p)) a b)]))))
          
          (1-?-swap
           (λ (p log?)
             (if log?
                 (lginv (fl p))
                 (inv- (fl p))))))])]
    ))




(: flstudent-t-sample : (case-> (Flonum               Integer -> FlVector)
                                (Flonum Flonum Flonum Integer -> FlVector)))
(define flstudent-t-sample
  (let ([F (λ ([ν : Flonum] [n : Integer])
             ; Note: Our gamma distribution has a shape parameter.
             ;       A shape parameter of 2 corresponds to a a rate of 1/2.
             (define Xs  (flnormal-sample 0. 1. n))
             (define X²s (flgamma-sample (/ ν 2.) 2. n))
             (flvector/ Xs (flvector-sqrt (flvector-scale X²s (/ ν)))))])
    (case-lambda
    ; X ~ t(ν)    
    [(ν n)
     (cond
       [(n . < . 0)  (raise-argument-error 'sample-student-t "Natural" 1 n)]
       [else    
        (F ν n)])]

    ; Y ~ σX+μ
    [(μ σ ν n)
     (cond
       [(n . < . 0)  (raise-argument-error 'sample-student-t "Natural" 3 n)]
       [else    
        (define Xs  (flnormal-sample 0. 1. n))
        (define X²s (flgamma-sample (/ ν 2.) 2. n))
        (build-flvector n
                        (λ (i)
                          (define X  (flvector-ref Xs  i))
                          (define X² (flvector-ref X²s i))
                          (define x  (fl/ X (flsqrt (fl/ X² ν))))
                          (fl+ (fl* σ x) μ)))])])))

